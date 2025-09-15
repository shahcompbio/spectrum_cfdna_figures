tf_zviran <- function(M = 2, #mutations detected
                      N = 1000, #ntotal mutations
                      mu = 0.001, #mean noise rate, denotes the mean noise rate (number of errors/number  of reads evaluated)
                      R = 1000, #total number of reads
                      cov = 1000 #local coverage in mutation sites
){
  Rp <- N * R
  TF <- 1 - (1 - (M - mu * Rp) / N) ^ (1 / cov)
  return(TF)
}

m_zviran <- function(TF = 0.01, 
                     N = 1000, 
                     mu = 0.001, 
                     R = 1000, 
                     cov = 1000){
  Rp <- N * R
  N * (1 - (1 - TF)^cov) + mu*Rp
}

tlod <- function(cov, N){
  x <- 0.99 / N
  return(1 - (1 - x)^(1/cov))
}
