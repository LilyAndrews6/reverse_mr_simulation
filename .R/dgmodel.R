##Functions for data generating model
gx_to_gp <- function(gx, h2x, prev){
  x_prime <- qnorm(prev, 0, 1)
  p <- pnorm(x_prime, mean=-gx, sd=sqrt(1-h2x), lower.tail = TRUE)
  return(p)
}

dgmodel <- function(
    nid, 
    b_gcc0,
    b_u2c0,
    gr_maf,
    nsnp,
    h2_known,
    h2_unknown,
    b_c0l,
    b_u2l,
    rsq_gr0,
    rsq_lr0,
    b_u1r0,
    b_u1l,
    d_prev,
    b_c0c1,
    b_u2c1,
    b_dc1,
    b_r0r1,
    b_u1r1,
    b_dr1,
    b_dx1, 
    b_x0x1, 
    gprs_rsq, 
    gprs_maf
)  #assuming variance of genetic liability is 1
{
  u1 <- rnorm(nid) #normal distribution of unmeasured confounder
  
  u2 <- rnorm(nid) #normal distribution of unmeasured confounder
  
  gwashits <- tribble(
    ~rsq, ~af,
    0.0007677501, 0.2141#0.351788, 0.0911 #rs75061358
  ) #known causal variants from glioma gwas data
  
  gc <- sapply(gwashits$af, \(x) rbinom(nid, 2, x)) #create genotype matrix for causal variant
  
  c0 <- gc * b_gcc0 + u2*b_u2c0 + rnorm(nid, sd=sqrt(1-b_gcc0^2-b_u2c0^2))
  
  x0 <- rnorm(nid) #consequence of disease protein pre-disease
  
  #gr_maf <- runif(nsnp, 0.01, 0.99)
  gprs_maf <- gprs_maf#rep(gr_maf, nsnp) #create minor allele frequency for each snp
  
  prs_rsq <- gprs_rsq #randomly generated betas for prs
  
  gprs <- make_geno(nid, nsnp, gprs_maf) #create genotype matrix for prs
  
  prs_known <- (gprs %*% prs_rsq) %>% {scale(.) * sqrt(h2_known)} #generate known prs
  
  prs_unknown <- rnorm(nid, 0, sd = sqrt(h2_unknown)) #generate unknown prs
  
  prs <- prs_known + prs_unknown #generate total prs
  
  rsq_prs <- h2_known - (0.1)^2 #generate r-squared prs
  
  l <- c0 * b_c0l + scale(prs) * sqrt(rsq_prs)+ u2 * b_u2l + rnorm(nid, 0, sd = 1 - sqrt(rsq_prs + b_c0l^2 + b_u2l^2)) #total genetic liability to disease
  
  gr <- rbinom(nid, 2, gr_maf) #creating genotype matrix for non causal variant
  
  r0 <- scale(gr) * sqrt(rsq_gr0) + scale(l) * sqrt(rsq_lr0) + u1 * b_u1r0 + rnorm(nid, sd=sqrt(1 - rsq_gr0 - rsq_lr0 - b_u1r0^2)) #create non causal variable
  
  prob_l <- gx_to_gp(gx=scale(l), h2x=rsq_prs + b_c0l^2 + b_u1l^2 + b_u2l^2, prev = d_prev) #generating the probability of liability to disease
  
  d <- rbinom(nid, 1, prob_l) #generate disease in individuals based on liability
  
  c1 <- c0 * b_c0c1 + u2 * b_u2c1 + d * b_dc1 #generate causal variant post disease
  
  r1 <- r0 * b_r0r1 + u1 * b_u1r1 + d * b_dr1 #generate non-causal variant post disease
  
  x1 <- x0 * b_x0x1 + d * b_dx1 + rnorm(nid, sd=sqrt(1 - b_dx1^2*var(d)))  #protein as a consequence of disease - error needs addition of variance of disease
  
  join_mat <- cbind(gprs, gr, gc) #combine genotypes
  
  phen <- tibble(
    u1, u2, r0, r1, c0, c1, prs, l, prob_l, d, x1, x0
  )  
  return(list(geno_gc=gc,geno_prs=gprs, geno_gr=gr, phen=phen, geno_join=join_mat, af=gprs_maf))}
