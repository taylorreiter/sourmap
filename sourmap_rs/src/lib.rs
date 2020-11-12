use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::{Path, PathBuf};
use std::sync::atomic::{AtomicUsize, Ordering};

use log::info;
use nohash_hasher::BuildNoHashHasher;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use sourmash::signature::Signature;
use sourmash::sketch::minhash::KmerMinHash;
use sourmash::sketch::Sketch;

type HashToIdx = HashMap<u64, HashMap<usize, u64>, BuildNoHashHasher<u64>>;

#[derive(Serialize, Deserialize)]
pub struct RevIndex {
    hash_to_idx: HashToIdx,
    sig_files: Vec<PathBuf>,
    ref_sigs: Option<Vec<Signature>>,
    template: Sketch,
}

impl RevIndex {
    pub fn load<P: AsRef<Path>>(
        index_path: P,
        queries: Option<&[KmerMinHash]>,
    ) -> Result<RevIndex, Box<dyn std::error::Error>> {
        // TODO: avoid loading full revindex if query != None
        let (rdr, _) = niffler::from_path(index_path)?;
        let mut revindex: RevIndex = serde_json::from_reader(rdr)?;

        if let Some(qs) = queries {
            for q in qs {
                let hashes: HashSet<u64> = q.iter_mins().cloned().collect();
                revindex.hash_to_idx.retain(|hash, _| hashes.contains(hash));
            }
        }
        Ok(revindex)
    }

    pub fn new(
        search_sigs: &[PathBuf],
        template: &Sketch,
        threshold: usize,
        queries: Option<&[KmerMinHash]>,
        keep_sigs: bool,
    ) -> RevIndex {
        let processed_sigs = AtomicUsize::new(0);

        // If threshold is zero, let's merge all queries and save time later
        let merged_query = if let Some(qs) = queries {
            if threshold == 0 {
                let mut merged = qs[0].clone();
                for query in &qs[1..] {
                    merged.merge(query).unwrap();
                }
                Some(merged)
            } else {
                None
            }
        } else {
            None
        };

        let hash_to_idx = search_sigs
            .par_iter()
            .enumerate()
            .filter_map(|(dataset_id, filename)| {
                let i = processed_sigs.fetch_add(1, Ordering::SeqCst);
                if i % 1000 == 0 {
                    info!("Processed {} reference sigs", i);
                }

                let mut search_mh = None;
                let search_sig = Signature::from_path(&filename)
                    .unwrap_or_else(|_| panic!("Error processing {:?}", filename))
                    .swap_remove(0);
                if let Some(sketch) = search_sig.select_sketch(&template) {
                    if let Sketch::MinHash(mh) = sketch {
                        search_mh = Some(mh);
                    }
                }
                let search_mh = search_mh.unwrap();

                let mut hash_to_idx = HashToIdx::with_hasher(BuildNoHashHasher::default());
                let mut add_to = |matched_hashes: Vec<(u64, u64)>, intersection| {
                    if !matched_hashes.is_empty() || intersection > threshold as u64 {
                        matched_hashes.into_iter().for_each(|(hash, count)| {
                            let mut dataset_ids = HashMap::new();
                            dataset_ids.insert(dataset_id, count);
                            hash_to_idx.insert(hash, dataset_ids);
                        });
                    }
                };

                if let Some(qs) = queries {
                    if let Some(ref merged) = merged_query {
                        let (matched_hashes, intersection) =
                            merged.intersection(search_mh).unwrap();
                        let matched_hashes: HashSet<_> = matched_hashes.into_iter().collect();
                        let matched_counts = search_mh
                            .to_vec_abunds()
                            .into_iter()
                            .filter_map(|(h, c)| {
                                if matched_hashes.contains(&h) {
                                    Some((h, c))
                                } else {
                                    None
                                }
                            })
                            .collect();
                        add_to(matched_counts, intersection);
                    } else {
                        for query in qs {
                            let (matched_hashes, intersection) =
                                query.intersection(search_mh).unwrap();
                            let matched_hashes: HashSet<_> = matched_hashes.into_iter().collect();
                            let matched_counts = search_mh
                                .to_vec_abunds()
                                .into_iter()
                                .filter_map(|(h, c)| {
                                    if matched_hashes.contains(&h) {
                                        Some((h, c))
                                    } else {
                                        None
                                    }
                                })
                                .collect();
                            add_to(matched_counts, intersection);
                        }
                    }
                } else {
                    let matched = search_mh.to_vec_abunds();
                    let size = matched.len() as u64;
                    add_to(matched, size);
                };

                if hash_to_idx.is_empty() {
                    None
                } else {
                    Some(hash_to_idx)
                }
            })
            .reduce(
                || HashToIdx::with_hasher(BuildNoHashHasher::default()),
                |a, b| {
                    let (small, mut large) = if a.len() > b.len() { (b, a) } else { (a, b) };

                    small.into_iter().for_each(|(hash, ids)| {
                        let entry = large.entry(hash).or_insert_with(HashMap::new);
                        for (id, count) in ids {
                            entry.insert(id, count);
                        }
                    });

                    large
                },
            );

        // TODO: build this together with hash_to_idx?
        let ref_sigs = if keep_sigs {
            Some(
                search_sigs
                    .par_iter()
                    .map(|ref_path| {
                        Signature::from_path(&ref_path)
                            .unwrap_or_else(|_| panic!("Error processing {:?}", ref_path))
                            .swap_remove(0)
                    })
                    .collect(),
            )
        } else {
            None
        };

        RevIndex {
            hash_to_idx,
            sig_files: search_sigs.into(),
            ref_sigs,
            template: template.clone(),
        }
    }

    pub fn template(&self) -> Sketch {
        self.template.clone()
    }

    pub fn abundance_csv<P: AsRef<Path>>(
        &self,
        output: P,
    ) -> Result<(), Box<dyn std::error::Error>> {
        // FIXME: we actually want hashes as columns, datasets as rows
        //        (but it is much easier to do hashes as rows, datasets as columns =] )
        let mut out = BufWriter::new(File::create(output).unwrap());

        // TODO: write header: hash, <enumerate dataset_ids>

        let total_ids = self.sig_files.len();

        for (hash, dataset_ids) in &self.hash_to_idx {
            write!(out, "{}", hash)?;
            let mut current_id = 0;
            let mut sorted_ids: Vec<_> = dataset_ids.iter().collect();
            sorted_ids.sort();

            for (dataset_id, count) in sorted_ids {
                //dbg!(dataset_id, count, current_id);
                while current_id != *dataset_id {
                    current_id += 1;
                    if current_id >= total_ids {
                        break;
                    } else {
                        write!(out, ",0")?;
                    }
                }
                if current_id == *dataset_id {
                    write!(out, ",{}", count)?;
                }
            }
            for _ in current_id..total_ids {
                write!(out, ",0")?;
            }

            writeln!(out)?;
        }
        Ok(())
    }
}
