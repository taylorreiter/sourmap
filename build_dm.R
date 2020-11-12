## as implemented here: https://doi.org/10.6084/m9.figshare.12864011.v4
## HELPER FUNCTIONS
## standardize columns func
norm_mat <- function(mat){
  ms <- apply(mat, 2, function(col) (col - mean(col)) / sd(col))
  ms
}
## compute euclid dist and standrdize to similarities
get_euc <- function(mat, n.threads = 1, alt = 'euclidean'){
  if(n.threads > 1){
    eu <- parDist(mat, method = alt, threads = n.threads) %>% as.matrix()
    ieu <- 1 / eu
    diag(ieu) <- 0
  }
  
  if(n.threads == 1){
    eu <- dist(mat, method = alt) %>% as.matrix()
    ieu <- 1 / eu
    diag(ieu) <- 0
  }
  ieu
}
## function to threshold the normalized distance mat
threshold <- function(mat, top_k = 10){
  thr <- mat
  tnr <- nrow(thr)
  ## set similarities outside of the top k to 0
  for(i in 1:tnr){
    ## rank entries in each row in reverse order (so largest value == rank 1), 
    ## and set entries that are outside of 1:k to 0
    thr[i, !rank(-thr[i, ], ties.method = 'random') %in% 1:top_k] <- 0
  }
  
  for(i in 1:tnr){
    for(j in 1:tnr){
      if(thr[i, j] == 0 & thr[j, i] != 0){
        thr[i, j] <- thr[j, i]
      }
    }
  }
  thr
}
## function to calculate norm laplacian
get_laplac <- function(mat){
  L <- -mat
  S <- rowSums(mat)
  nL <- L / S
  diag(nL) <- 1
  nL
}

## build diffusion map
build_dm <- function(mat, k_nbhrs = 10, keep_eig = 10){
  require(magrittr)
  ## as implemented here: https://doi.org/10.6084/m9.figshare.12864011.v4
  nm <- mat %>% norm_mat()        # normalize
  aff <- nm %>%
    get_euc(., n.threads = 1) %>% 
    threshold(., top_k = k_nbhrs) # make affinity matrix
  Lij <- aff %>% get_laplac()     # compute laplacian
  eig <- eigen(Lij)               # smallest keig vectors
  evl <- eig$values %>%           # get eigenvalues
    Re() %>%
    round(., digits = 10)
  evc <- eig$vectors %>%          # get eigenvectors
    Re() %>%
    round(., digits = 10)
  
  # create objects containing diffusion components
  for(d in 1:length(evl)){
    assign(paste('dim', d, sep = '_'),
           Re(evc[, rank(evl) == (d + 1)])
    )
  }
  k_eig <- keep_eig               # specify num of variables you want to keep
  dat <- do.call(mapply, c(FUN = cbind, mget(paste0("dim_", 1:(k_eig))))) %>%
    t() %>%
    as.data.frame()               # merge in array
  
  # make column names pretty
  colnames(dat) <- paste0("DC", 1:ncol(dat))
  return(dat)
}