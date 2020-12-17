use std::fs::File;
use std::io::{BufRead, BufReader};
use std::path::{Path, PathBuf};

use log::info;
use structopt::StructOpt;

use sourmash::sketch::minhash::{max_hash_for_scaled, KmerMinHash};
use sourmash::sketch::Sketch;

use sourmap::RevIndex;

#[derive(StructOpt, Debug)]
enum Cli {
    Index {
        /// The path for output
        #[structopt(parse(from_os_str))]
        output: PathBuf,

        /// List of reference signatures
        #[structopt(parse(from_os_str))]
        siglist: PathBuf,

        /// ksize
        #[structopt(short = "k", long = "ksize", default_value = "31")]
        ksize: u8,

        /// scaled
        #[structopt(short = "s", long = "scaled", default_value = "1000")]
        scaled: usize,
    },
    AbundanceMatrix {
        /// Precomputed index or list of reference signatures
        #[structopt(parse(from_os_str))]
        siglist: PathBuf,

        /// ksize
        #[structopt(short = "k", long = "ksize", default_value = "31")]
        ksize: u8,

        /// scaled
        #[structopt(short = "s", long = "scaled", default_value = "1000")]
        scaled: usize,

        /// The path for output
        #[structopt(parse(from_os_str), short = "o", long = "output")]
        output: PathBuf,

        /// Is the index a list of signatures?
        #[structopt(long = "--from-file")]
        from_file: bool,
    },
}

fn read_paths<P: AsRef<Path>>(paths_file: P) -> Result<Vec<PathBuf>, Box<dyn std::error::Error>> {
    let paths = BufReader::new(File::open(paths_file)?);
    Ok(paths
        .lines()
        .map(|line| {
            let mut path = PathBuf::new();
            path.push(line.unwrap());
            path
        })
        .collect())
}

fn build_template(ksize: u8, scaled: usize) -> Sketch {
    let max_hash = max_hash_for_scaled(scaled as u64);
    let template_mh = KmerMinHash::builder()
        .num(0u32)
        .ksize(ksize as u32)
        .max_hash(max_hash)
        .build();
    Sketch::MinHash(template_mh)
}

fn index<P: AsRef<Path>>(
    siglist: P,
    template: Sketch,
    output: P,
) -> Result<(), Box<dyn std::error::Error>> {
    info!("Loading siglist");
    let index_sigs = read_paths(siglist)?;
    info!("Loaded {} sig paths in siglist", index_sigs.len());

    let revindex = RevIndex::new(&index_sigs, &template, 0, None, false);

    info!("Saving index");
    let wtr = niffler::to_path(
        output,
        niffler::compression::Format::Gzip,
        niffler::compression::Level::One,
    )?;
    serde_json::to_writer(wtr, &revindex)?;

    Ok(())
}

fn abundance_matrix<P: AsRef<Path>>(
    siglist: P,
    template: Sketch,
    output: P,
    from_file: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let revindex = if from_file {
        info!("Loading siglist");
        let search_sigs = read_paths(siglist)?;
        info!("Loaded {} sig paths in siglist", search_sigs.len());

        RevIndex::new(&search_sigs, &template, 0, None, true)
    } else {
        RevIndex::load(siglist, None)?
    };

    revindex.abundance_csv(output)?;

    Ok(())
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    env_logger::Builder::from_env(env_logger::Env::default().default_filter_or("info")).init();

    match Cli::from_args() {
        Cli::Index {
            output,
            siglist,
            ksize,
            scaled,
        } => {
            let template = build_template(ksize, scaled);

            index(siglist, template, output)?
        }
        Cli::AbundanceMatrix {
            siglist,
            ksize,
            scaled,
            output,
            from_file,
        } => {
            let template = build_template(ksize, scaled);

            abundance_matrix(siglist, template, output, from_file)?;
        }
    };

    Ok(())
}
