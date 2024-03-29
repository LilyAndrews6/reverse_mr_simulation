Sims Functions
================
Lily Andrews
2023-11-28

Create data generating model (dgmodel)

``` r
#' Simulation model of DAG for the role of causal and non-causal biomarkers in relation to liability to disease
#'
#' @param nid number of individuals
#' @param nsnp number of SNPs
#' @param gc_maf causal biomarker minor allele frequency input
#' @param b_gl genetic liability from discovered SNPs
#' @param rsq_prs h-squared of y explained
#' @param rsq_z H-squared of y when combined with PRS
#' @param d_prev disease prevalence
#' @param b_lr0 beta of disease liability to non-causal biomarker (pre-diagnosis)
#' @param b_c0l beta of causal biomarker (pre-diagnosis) to disease liability
#' @param gr_maf non-causal biomarker minor allele frequency input
#' @param rsq_gr0 r-squared value for the non-causal pQTL 
#' @param rsq_gc0 r-squared value for the causal pQTL
#' @param b_u1r1 beta of unmeasured confounder 1 to non-causal marker (post-diagnosis)
#' @param b_u1l beta of unmeasured confounder 1 to disease liability
#' @param b_u2c1 beta of unmeasured confounder 2 to causal biomarker (post diagnosis)
#' @param b_u2l beta of unmeasured confounder 2 to disease liability
#' @param b_dr1 beta of disease to non-causal biomarker (post diagnosis)
#' @param b_dc1 beta of disease to causal biomarker (post diagnosis)
#' @param b_c0c1 beta of causal biomarker (pre diagnosis) to causal biomarker (post diagnosis)
#' @param b_r0r1 beta of non-causal biomarker (pre diagnosis) to non-causal biomarker (post diagnosis)
#' @param b_u2c0 beta of unmeasured confounder 2 to causal biomarker (pre diagnosis)
#' @param b_u1r0 beta of unmeasured confounder 1 to non-causal biomarker (pre diagnosis)
#' @param b_gcc0 beta of causal genotype matrix to causal biomarker (pre diagnosis)
#' @param b_lc0 ##CONFUSED BY THIS AS IT WOULD BE GOING IN THE WRONG DIRECTION SHOULD IT BE b_c0l which is already mentioned before
#' @param b_u1d ##would this not be b_u1l*b_ld
#' @param b_u2d ##would this not be b_u2l*b_ld
#'
#' @return
#' @export
#'
#' @examples

dgmodel <- function(nid, nsnp, gc_maf, b_gl, rsq_prs, rsq_z, d_prev, b_lr0, b_c0l, gr_maf, rsq_gr0, rsq_gc0, b_u1r1, b_u1l, b_u2c1, b_u2l, b_dr1, b_dc1, b_c0c1, b_r0r1, b_u2c0, b_u1r0, b_gcc0, b_lc0, b_u1d, b_u2d, b_gcc1)
{
  u1 <- rnorm(nid) #normal distribution of unmeasured confounder
  u2 <- rnorm(nid) #normal distribution of unmeasured confounder
  gc <- make_geno(nid, nsnp, gc_maf) #create genotype matrix for causal variant
  gr <- rbinom(nid, 2, gr_maf) #creating genotype matrix for non causal variant 
  # TODO need to fix this
  c0 <- scale(gc[,1]) * b_gcc0 + u2 * b_u2c0 + rnorm(nid, sd=sqrt(1-rsq_gc0-b_lc0^2-b_u2c0^2)) #we hypothesised one SNP down this route? no need for unmeasured confounder in this case? rnorm introduces the variance into the equation. We are assuming variance is 1
   
  
  # TODO need to fix this
  rsq_prs <- rsq_prs - (0.1)^2 # add on the environmental aspect?
  prs <- scale(gc %*% b_gl) #b_gl does this mean beta genetic liability or route PRS to disease liability, why do we include this line if prs_w is already used
  prs_w <- scale(gc[,-1] %*% b_gl[-1]) #had to remove the SNP as the liability calculated later on would include an extra SNP
  z <- rnorm(nid) #included the rest of the known genetic variants
  l <- prs_w * sqrt(rsq_prs) + z * sqrt(rsq_z) + c0 * b_c0l + u1 * b_u1l + u2 * b_u2l #total genetic liability to disease, decided not to add error into liability to avoid diagram c from happening https://www.nature.com/articles/nrg3377/figures/1  
  r0 <- scale(gr) * sqrt(rsq_gr0) + l * b_lr0 + u1 * b_u1r0 + rnorm(nid, sd=sqrt(1-rsq_gr0-b_lr0^2-b_u1r0^2)) #  We are assuming variance is 1

  # generate probability of disease
  prob_l <- simulateGP::gx_to_gp(gx=scale(l), h2x=rsq_prs + rsq_z + b_c0l^2 + b_u1d^2 + b_u2d^2, prev = d_prev) #translate disease risk from liability to probability scale would h2x be prs_w instead?
  #change b_u1d to b_u1l and same with u2
  
  # generate random disease outcome
  # switch this 0/1 or fix simualteGP function
  d <- rbinom(nid, 1, prob_l) #case as 0 control is 1
  ##can change to 1-prob_l to flip this removing next line
  d <- abs(d-1) #this function flips case and control which makes case 1 and control 0
  #response to disease state
  c1 <- scale(gc[,1]) * b_gcc1 + u2 * b_u2c1 + c0 *b_c0c1 + d * b_dc1 + rnorm(nid, sd=sqrt(1-(b_gcc1)^2 - b_u2c1^2- b_c0c1^2)) #we hypothesised one SNP down this route? We are assuming variance is 1
  
  r1 <- scale(gr) * sqrt(rsq_gr0) + r0 * b_r0r1 + u1 * b_u1r1 + d * b_dr1 + rnorm(nid, sd=sqrt(1-rsq_gc0-b_lr0^2-b_u1r1^2)) #We are assuming variance is 1

  phen <- tibble(
    u1, u2, r0, r1, c0, c1, prs, prs_w, l, prob_l, d
  )
  return(list(geno=gc, phen=phen)) #check genotype would be Gc when it was G
}
```

think about the one snp of PC variable

List of expected variances from the data generating model: v_u1 ~ N(0,1)
v_u2 = 1 v_l = 1 v_pc = 1 v_pr = 1 v_z = 1 v_d = Binomial(n=1, p=d_prev
\* (1-d_prev))

Check output of data generating model (dgmodel_check)

``` r
dgmodel_check <- function(dat)
{
  print(tibble(
    col = names(dat$phen),
    vars = apply(dat$phen, 2, var), ##what does var do here, does it look for the variance between values? 2 means col sums and 1 is row sums
    means = colMeans(dat$phen)
  ))
  print(cor(dat$phen$prs, dat$phen$l)^2) #compare PRS to disease liability
}
```

Assess observational, case control and protein GWAS data
(dgmodel_analysis)

``` r
#looking at different scenerios, all, case control, protein gwas
dgmodel_analysis <- function(dat, ncase, ncontrol, protein_gwas)
{
  d <- dat$phen
  res_all <- bind_rows( 
    fast_assoc(y=d$d, x=d$prs) %>% as_tibble() %>% mutate(x="prs", y="d"), #fast_assoc simple linear regression for prs and disease outcome
    fast_assoc(y=d$c0, x=d$d) %>% as_tibble() %>% mutate(x="c0", y="d"),
    fast_assoc(y=d$r0, x=d$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=d$prs, x=d$c0) %>% as_tibble() %>% mutate(x="prs", y="c0"),
    fast_assoc(y=d$prs, x=d$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=d$c1, x=d$d) %>% as_tibble() %>% mutate(x="c1", y="d"),
    fast_assoc(y=d$r1, x=d$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=d$prs, x=d$c1) %>% as_tibble() %>% mutate(x="prs", y="c1"),
    fast_assoc(y=d$prs, x=d$r1) %>% as_tibble() %>% mutate(x="prs", y="r1")
  ) %>%
    mutate(study="all")
  
  # note - need to switch this once 0/1 is fixed
  obs <- bind_rows(d[d$d == 1,][1:ncontrol,], d[d$d == 0,][1:ncase,]) ##shouldn't this be the other way round as case is 1 and control 0
  
  res_obs <- bind_rows(
    fast_assoc(y=obs$d, x=obs$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=obs$c0, x=obs$d) %>% as_tibble() %>% mutate(x="c0", y="d"),
    fast_assoc(y=obs$r0, x=obs$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=obs$prs, x=obs$c0) %>% as_tibble() %>% mutate(x="prs", y="c0"),
    fast_assoc(y=obs$prs, x=obs$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=obs$c1, x=obs$d) %>% as_tibble() %>% mutate(x="c1", y="d"),
    fast_assoc(y=obs$r1, x=obs$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=obs$prs, x=obs$c1) %>% as_tibble() %>% mutate(x="prs", y="c1"),
    fast_assoc(y=obs$prs, x=obs$r1) %>% as_tibble() %>% mutate(x="prs", y="r1"),
  )  %>%
    mutate(study="observational")

  d <- d[1:protein_gwas, ]
  res_protein <- bind_rows(
    fast_assoc(y=d$d, x=d$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=d$c0, x=d$d) %>% as_tibble() %>% mutate(x="c0", y="d"),
    fast_assoc(y=d$r0, x=d$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=d$prs, x=d$c0) %>% as_tibble() %>% mutate(x="prs", y="c0"),
    fast_assoc(y=d$prs, x=d$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=d$c1, x=d$d) %>% as_tibble() %>% mutate(x="c1", y="d"),
    fast_assoc(y=d$r1, x=d$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=d$prs, x=d$c1) %>% as_tibble() %>% mutate(x="prs", y="c1"),
    fast_assoc(y=d$prs, x=d$r1) %>% as_tibble() %>% mutate(x="prs", y="r1"),
  ) %>%
    mutate(study="protein_gwas")
  
  return(bind_rows(res_all, res_obs, res_protein))
}
```

See reverse MR approach to find non causal markers (ivw_analysis)

``` r
#reverse MR approach using prs to find non causal biomarkers
#IVW and two-stage least squares (should be the same value)
ivw_analysis <- function(dat){
  d <- dat$phen
  gen <- dat$geno_causal
dat <- get_effs(d$prs, d$r0, gen) #might need to change this genetic matrix
print("MR")
mr(dat, metho="mr_ivw") %>% str()
# MR using fixed effects IVW [two-stage least squares method] use system fit package
print("2SLS")
summary(systemfit::systemfit(d$prs ~ d$r0, method="2SLS", inst = ~gen)) #2SLS #might need to change this genetic matrix ##issue with this code as getting NAs
print("Observational")
summary(lm(d$prs ~ d$r0)) #confounded observational estimate
}
```

Parallel

``` r
dgmodel_time <- function(nid=100000, 
  nsnp=99, 
  gc_maf=runif(99, 0.05, 0.95), 
  b_gl=rnorm(99), 
  rsq_prs=0.06, 
  rsq_z=0.19, 
  d_prev=5/1000, 
  b_lr0=sqrt(0.1), 
  b_c0l=sqrt(0.01), 
  gr_maf=0.5, 
  rsq_gr0=0.1,   
  rsq_gc0=0.1,   
  b_u1r1=0.1, 
  b_u1l=0.1, 
  b_u2c1=0.1, 
  b_u2l=0.1, 
  b_dr1=0.1, 
  b_dc1=0.1, 
  b_c0c1=0.1, 
  b_r0r1=0.1, 
  b_u2c0=0.1, 
  b_gcc0=0.1, 
  b_lc0=0.1, 
  b_u1r0=0.1, 
  b_u1d=0.1, 
  b_u2d=0.1)
{
  u1 <- rnorm(nid) 
  u2 <- rnorm(nid) 
  gc <- make_geno(nid, nsnp, gc_maf) 
  gr <- rbinom(nid, 2, gr_maf) 
  c0 <- scale(gc[,1]) * b_gcc0 + u2 * b_u2c0 + rnorm(nid, sd=sqrt(1-rsq_gc0-b_lc0^2-b_u2c0^2)) 
  rsq_prs <- rsq_prs - (0.1)^2 
  prs <- scale(gc %*% b_gl) 
  prs_w <- scale(gc[,-1] %*% b_gl[-1]) 
  z <- rnorm(nid) 
  l <- prs_w * sqrt(rsq_prs) + z * sqrt(rsq_z) + c0 * b_c0l + u1 * b_u1l + u2 * b_u2l 
  r0 <- scale(gr) * sqrt(rsq_gr0) + l * b_lr0 + u1 * b_u1r0 + rnorm(nid, sd=sqrt(1-rsq_gr0-b_lr0^2-b_u1r0^2)) 
  prob_l <- simulateGP::gx_to_gp(gx=scale(l), h2x=rsq_prs + rsq_z + b_c0l^2 + b_u1d^2 + b_u2d^2, prev = d_prev) 
  d <- rbinom(nid, 1, prob_l)
  d <- abs(d-1) 
  c1 <- scale(gc[,1]) * 0.1 + u2 * b_u2c1 + c0 *b_c0c1 + d * b_dc1 + rnorm(nid, sd=sqrt(1-(0.1)^2 - b_u2c1^2- b_c0c1^2)) 
  r1 <- scale(gr) * sqrt(rsq_gr0) + r0 * b_r0r1 + u1 * b_u1r1 + d * b_dr1 + rnorm(nid, sd=sqrt(1-rsq_gc0-b_lr0^2-b_u1r1^2)) 
  phen <- tibble(
    u1, u2, r0, r1, c0, c1, prs, prs_w, l, prob_l, d
  )
  return(list(geno=gc, phen=phen)) 
}
```

Parallel for check

``` r
check_time <- function(dgmodel_time)
{
  print(tibble(
    col = names(dat$phen),
    vars = apply(dat$phen, 2, var), ##what does var do here, does it look for the variance between values? 2 means col sums and 1 is row sums
    means = colMeans(dat$phen)
  ))
  print(cor(dat$phen$prs, dat$phen$l)^2) #compare PRS to disease liability
}
```

For this DAG we haven’t included parameter Z of unknown GWAS variants,
we include all GWAS hits at the same time and look at each reverse MR
result one at a time

``` r
dgmodel <- function(nid, #3301 for protein and 30686 for GWAS
  nsnp,
  b_gl, 
  h2_known, #h-squared value https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667278/
  h2_unknown, #to total 0.25 of total H-squared https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667278/
  d_prev, #~5% http://rarecarenet.istitutotumori.mi.it/analysis.php
  rsq_lr0, 
  b_c0l, 
  gc_maf,
  gr_maf, 
  rsq_gr0,   
  rsq_gc0,   
  b_u1r1, 
  b_u1l, 
  b_u2c1, 
  b_u2l, 
  b_dr1, 
  b_dc1, 
  b_c0c1, 
  b_r0r1, 
  b_u2c0, 
  b_gcc0, 
  b_lc0, 
  b_u1r0, 
  b_u1d, 
  b_u2d, 
  b_gcc1,
  var_l, 
  vp)  #assuming variance of genetic liability is 1
{
  u1 <- rnorm(nid) #normal distribution of unmeasured confounder
  u2 <- rnorm(nid) #normal distribution of unmeasured confounder
gwashits <- tribble(
    ~beta, ~af,
    0.26, 0.54,
    0.37, 0.30,
    0.15, 0.30,
    0.35, 0.09,
    0.21, 0.84,
    0.16, 0.61,
    0.69, 0.08,
    0.25, 0.44,
    0.26, 0.46,
    0.13, 0.33,
    0.16, 0.67,
    0.93, 0.01, 
    0.29, 0.79,
    0.31, 0.80,
    0.16, 0.88,
    0.13, 0.24,
    0.13, 0.73
)
 gc <- sapply(gwashits$af, \(x) rbinom(nid, 2, x)) #creating genotype matrix for causal variant based on GWAS data 

rsq_gcl <- 2*gwashits$af*(1-gwashits$af) * gwashits$beta^2 / var_l
#below we assume equal split of rsq to make gc -> l (route is gc ->c0 ->l)
rsq_gcc0 <- sqrt(rsq_gcl)
rsq_c0l <- sqrt(rsq_gcl)
# rsq_gp = beta^2 var(g) / var(p)
Var_c0 <- rep(1, nrow(gwashits))
 
beta_gcc0 <- sqrt(rsq_gcl * Var_c0 / (2*gwashits$af*(1-gwashits$af)))
beta_c0l <- sqrt(rsq_c0l * var_l / Var_c0)

# gc -> c0
c0 <- gc * b_gcc0 + u2*b_u2c0 + rnorm(nid, sd=sqrt(1-b_gcc0^2-b_u2c0^2))
#c0 -> l
#for now assuming unmeasured confounder is one estimate for all GWAS hits
l <- gc %*% gwashits$beta +u2*b_u2l #remove error element to liability as this is accounted for in the h2x estimate

gr <- rbinom(nid, 2, gr_maf) #creating genotype matrix for non causal variant 
  #to calculate c0 using known GWAS hits
  # TODO need to fix this
  rsq_prs <- h2_known - (0.1)^2 # add on the environmental aspect?
  
  #vp <- 1 #Vp <- rep(1, length(b_gl)) #don't think we need this as it combines SNP values rather than take them individually
  prs <- rnorm(nid, mean=0, sd=sqrt(vp*h2_known)) #check this as vp*h2 is rsq_prs, original code: sd=sqrt(vp*h2) h2 is rsq_prs currently as explains 0.06 of heritability and total heritability is 0.025
  
  #num_c0 <- length(gwashits$af) #number of gwas hits won't need
  
  # won't need to remove gwashits as PRS is total together prs_w <- scale(prs[,-num_c0] %*% b_gl[-num_c0]) #had to remove number of GWAS hits from prs
  l <- prs * sqrt(rsq_prs) + l #combine PRS and liability from known GWAS hits 
   
  r0 <- scale(gr) * sqrt(rsq_gr0) + l * rsq_lr0 + u1 * b_u1r0 + rnorm(nid, sd=sqrt(1-rsq_gr0-rsq_lr0^2-b_u1r0^2)) #  We are assuming variance is 1

  # generate probability of disease
  prob_l <- simulateGP::gx_to_gp(gx=scale(l), h2x=rsq_prs + b_c0l^2 + b_u1d^2 + b_u2d^2, prev = d_prev) #translate disease risk from liability to probability scale would h2x be prs_w instead?
  #change b_u1d to b_u1l and same with u2
  
  # generate random disease outcome
  # switch this 0/1 or fix simualteGP function
  d <- rbinom(nid, 1, prob_l) #case as 0 control is 1
  ##can change to 1-prob_l to flip this removing next line
  d <- abs(d-1) #this function flips case and control which makes case 1 and control 0

  c1 <- scale(gc) * b_gcc1 + u2 * b_u2c1 + c0 *b_c0c1 + d * b_dc1 + rnorm(nid, sd=sqrt(1-(b_gcc1)^2 - b_u2c1^2- b_c0c1^2)) #we hypothesized one SNP down this route? We are assuming variance is 1
    #response to disease state NEED TO ADD ALTERATION OF DISEASE STATE WHETHER THIS IS INCREASED OR DECREASED
  
  r1 <- scale(gr) * sqrt(rsq_gr0) + r0 * b_r0r1 + u1 * b_u1r1 + d * b_dr1 + rnorm(nid, sd=sqrt(1-rsq_gc0-rsq_lr0^2-b_u1r1^2)) #We are assuming variance is 1
#have left this for now but do we need to add element of change for non-causal markers in disease process
   phen <- tibble(
    u1, u2, r0, r1, c0[,1],c0[,2],c0[,3],c0[,4],c0[,5],c0[,6],c0[,7],c0[,8], c0[,9], c0[,10], c0[,11], c0[,12], c0[,13], c0[,14], c0[,15], c0[,16], c0[,17], c1[,1], c1[,2], c1[,3], c1[,4], c1[,5], c1[,6], c1[,7], c1[,8], c1[,9], c1[,10], c1[,11], c1[,12], c1[,13], c1[,14], c1[,15], c1[,16], c1[,17], prs, l, prob_l, d
  )
  return(list(geno=gc, phen=phen))} #check genotype would be Gc when it was G

dgmodel_check <- function(dat)
{
  print(tibble(
    col = names(dat$phen),
    vars = apply(dat$phen, 2, var), ##what does var do here, does it look for the variance between values? 2 means col sums and 1 is row sums
    means = colMeans(dat$phen)
  ), n=25)
  print(cor(dat$phen$prs, dat$phen$l)^2) #compare PRS to disease liability
}
dgmodel_analysis <- function(dat, ncase, ncontrol, protein_gwas, rev_mr_result)
{
  phen <- dat$phen
  res_all <- bind_rows( 
    fast_assoc(y=phen$d, x=phen$prs) %>% as_tibble() %>% mutate(x="prs", y="d"), #fast_assoc simple linear regression for prs and disease outcome
    fast_assoc(y=phen$`c0[, 1]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_1", y="d"),
    fast_assoc(y=phen$`c0[, 2]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_2", y="d"),
    fast_assoc(y=phen$`c0[, 3]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_3", y="d"),
    fast_assoc(y=phen$`c0[, 4]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_4", y="d"),
    fast_assoc(y=phen$`c0[, 5]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_5", y="d"),
    fast_assoc(y=phen$`c0[, 6]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_6", y="d"),
    fast_assoc(y=phen$`c0[, 7]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_7", y="d"),
    fast_assoc(y=phen$`c0[, 8]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_8", y="d"),
    fast_assoc(y=phen$`c0[, 9]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_9", y="d"),
    fast_assoc(y=phen$`c0[, 10]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_10", y="d"),
    fast_assoc(y=phen$`c0[, 11]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_11", y="d"),
    fast_assoc(y=phen$`c0[, 12]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_12", y="d"),
    fast_assoc(y=phen$`c0[, 13]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_13", y="d"),
    fast_assoc(y=phen$`c0[, 14]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_14", y="d"),
    fast_assoc(y=phen$`c0[, 15]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_15", y="d"),
    fast_assoc(y=phen$`c0[, 16]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_16", y="d"),
    fast_assoc(y=phen$`c0[, 17]`, x=phen$d) %>% as_tibble() %>% mutate(x="c0_17", y="d"),
    fast_assoc(y=phen$r0, x=phen$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c0_1"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 2]`) %>% as_tibble() %>% mutate(x="prs", y="c0_2"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 3]`) %>% as_tibble() %>% mutate(x="prs", y="c0_3"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 4]`) %>% as_tibble() %>% mutate(x="prs", y="c0_4"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 5]`) %>% as_tibble() %>% mutate(x="prs", y="c0_5"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 6]`) %>% as_tibble() %>% mutate(x="prs", y="c0_6"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 7]`) %>% as_tibble() %>% mutate(x="prs", y="c0_7"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 8]`) %>% as_tibble() %>% mutate(x="prs", y="c0_8"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 9]`) %>% as_tibble() %>% mutate(x="prs", y="c0_9"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 10]`) %>% as_tibble() %>% mutate(x="prs", y="c0_10"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 11]`) %>% as_tibble() %>% mutate(x="prs", y="c0_11"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 12]`) %>% as_tibble() %>% mutate(x="prs", y="c0_12"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 13]`) %>% as_tibble() %>% mutate(x="prs", y="c0_13"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 14]`) %>% as_tibble() %>% mutate(x="prs", y="c0_14"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 15]`) %>% as_tibble() %>% mutate(x="prs", y="c0_15"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 16]`) %>% as_tibble() %>% mutate(x="prs", y="c0_16"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 17]`) %>% as_tibble() %>% mutate(x="prs", y="c0_17"),
    fast_assoc(y=phen$prs, x=phen$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=phen$`c1[, 1]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_1", y="d"),
    fast_assoc(y=phen$`c1[, 2]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_2", y="d"),
    fast_assoc(y=phen$`c1[, 3]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_3", y="d"),
    fast_assoc(y=phen$`c1[, 4]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_4", y="d"),
    fast_assoc(y=phen$`c1[, 5]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_5", y="d"),
    fast_assoc(y=phen$`c1[, 6]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_6", y="d"),
    fast_assoc(y=phen$`c1[, 7]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_7", y="d"),
    fast_assoc(y=phen$`c1[, 8]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_8", y="d"), 
    fast_assoc(y=phen$`c1[, 9]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_9", y="d"), 
    fast_assoc(y=phen$`c1[, 10]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_10", y="d"),
    fast_assoc(y=phen$`c1[, 11]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_11", y="d"),
    fast_assoc(y=phen$`c1[, 12]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_12", y="d"),
    fast_assoc(y=phen$`c1[, 13]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_13", y="d"),
    fast_assoc(y=phen$`c1[, 14]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_14", y="d"),
    fast_assoc(y=phen$`c1[, 15]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_15", y="d"),
    fast_assoc(y=phen$`c1[, 16]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_16", y="d"),
    fast_assoc(y=phen$`c1[, 17]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_17", y="d"),
    fast_assoc(y=phen$r1, x=phen$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c1_1"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 2]`) %>% as_tibble() %>% mutate(x="prs", y="c1_2"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 3]`) %>% as_tibble() %>% mutate(x="prs", y="c1_3"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 4]`) %>% as_tibble() %>% mutate(x="prs", y="c1_4"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 5]`) %>% as_tibble() %>% mutate(x="prs", y="c1_5"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 6]`) %>% as_tibble() %>% mutate(x="prs", y="c1_6"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 7]`) %>% as_tibble() %>% mutate(x="prs", y="c1_7"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 8]`) %>% as_tibble() %>% mutate(x="prs", y="c1_8"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 9]`) %>% as_tibble() %>% mutate(x="prs", y="c1_9"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 10]`) %>% as_tibble() %>% mutate(x="prs", y="c1_10"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 11]`) %>% as_tibble() %>% mutate(x="prs", y="c1_11"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 12]`) %>% as_tibble() %>% mutate(x="prs", y="c1_12"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 13]`) %>% as_tibble() %>% mutate(x="prs", y="c1_13"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 14]`) %>% as_tibble() %>% mutate(x="prs", y="c1_14"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 15]`) %>% as_tibble() %>% mutate(x="prs", y="c1_15"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 16]`) %>% as_tibble() %>% mutate(x="prs", y="c1_16"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 17]`) %>% as_tibble() %>% mutate(x="prs", y="c1_17"),
    fast_assoc(y=phen$prs, x=phen$r1) %>% as_tibble() %>% mutate(x="prs", y="r1")
  ) %>%
    mutate(study="all")
  
  obs <- bind_rows(phen[phen$d == 0,][1:ncontrol,], phen[phen$d == 1,][1:ncase,]) ##this causes NA when cases doesn't have same number in individuals assessed in DGM
  case <- length(subset(obs$d, obs$d==1))
  paste("Number of cases: ", case) 
  control <- length(subset(obs$d, obs$d==0))
  paste("Number of controls: ", control)
  res_obs <- bind_rows(
    fast_assoc(y=obs$d, x=obs$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=obs$`c0[, 1]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_1", y="d"),
    fast_assoc(y=obs$`c0[, 2]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_2", y="d"),
    fast_assoc(y=obs$`c0[, 3]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_3", y="d"),
    fast_assoc(y=obs$`c0[, 4]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_4", y="d"),
    fast_assoc(y=obs$`c0[, 5]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_5", y="d"),
    fast_assoc(y=obs$`c0[, 6]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_6", y="d"),
    fast_assoc(y=obs$`c0[, 7]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_7", y="d"),
    fast_assoc(y=obs$`c0[, 8]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_8", y="d"),
    fast_assoc(y=obs$`c0[, 9]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_9", y="d"),
    fast_assoc(y=obs$`c0[, 10]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_10", y="d"),
    fast_assoc(y=obs$`c0[, 11]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_11", y="d"),
    fast_assoc(y=obs$`c0[, 12]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_12", y="d"),
    fast_assoc(y=obs$`c0[, 13]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_13", y="d"),
    fast_assoc(y=obs$`c0[, 14]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_14", y="d"),
    fast_assoc(y=obs$`c0[, 15]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_15", y="d"),
    fast_assoc(y=obs$`c0[, 16]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_16", y="d"),
    fast_assoc(y=obs$`c0[, 17]`, x=obs$d) %>% as_tibble() %>% mutate(x="c0_17", y="d"),
    fast_assoc(y=obs$r0, x=obs$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c0_1"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 2]`) %>% as_tibble() %>% mutate(x="prs", y="c0_2"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 3]`) %>% as_tibble() %>% mutate(x="prs", y="c0_3"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 4]`) %>% as_tibble() %>% mutate(x="prs", y="c0_4"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 5]`) %>% as_tibble() %>% mutate(x="prs", y="c0_5"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 6]`) %>% as_tibble() %>% mutate(x="prs", y="c0_6"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 7]`) %>% as_tibble() %>% mutate(x="prs", y="c0_7"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 8]`) %>% as_tibble() %>% mutate(x="prs", y="c0_8"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 9]`) %>% as_tibble() %>% mutate(x="prs", y="c0_9"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 10]`) %>% as_tibble() %>% mutate(x="prs", y="c0_10"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 11]`) %>% as_tibble() %>% mutate(x="prs", y="c0_11"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 12]`) %>% as_tibble() %>% mutate(x="prs", y="c0_12"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 13]`) %>% as_tibble() %>% mutate(x="prs", y="c0_13"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 14]`) %>% as_tibble() %>% mutate(x="prs", y="c0_14"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 15]`) %>% as_tibble() %>% mutate(x="prs", y="c0_15"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 16]`) %>% as_tibble() %>% mutate(x="prs", y="c0_16"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 17]`) %>% as_tibble() %>% mutate(x="prs", y="c0_17"),
    fast_assoc(y=obs$prs, x=obs$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=obs$`c1[, 1]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_1", y="d"),
    fast_assoc(y=obs$`c1[, 2]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_2", y="d"),
    fast_assoc(y=obs$`c1[, 3]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_3", y="d"),
    fast_assoc(y=obs$`c1[, 4]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_4", y="d"),
    fast_assoc(y=obs$`c1[, 5]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_5", y="d"),
    fast_assoc(y=obs$`c1[, 6]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_6", y="d"),
    fast_assoc(y=obs$`c1[, 7]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_7", y="d"),
    fast_assoc(y=obs$`c1[, 8]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_8", y="d"),
    fast_assoc(y=obs$`c1[, 9]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_9", y="d"),
    fast_assoc(y=obs$`c1[, 10]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_10", y="d"),
    fast_assoc(y=obs$`c1[, 11]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_11", y="d"),
    fast_assoc(y=obs$`c1[, 12]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_12", y="d"),
    fast_assoc(y=obs$`c1[, 13]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_13", y="d"),
    fast_assoc(y=obs$`c1[, 14]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_14", y="d"),
    fast_assoc(y=obs$`c1[, 15]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_15", y="d"),
    fast_assoc(y=obs$`c1[, 16]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_16", y="d"),
    fast_assoc(y=obs$`c1[, 17]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_17", y="d"),
    fast_assoc(y=obs$r1, x=obs$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c1_1"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 2]`) %>% as_tibble() %>% mutate(x="prs", y="c1_2"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 3]`) %>% as_tibble() %>% mutate(x="prs", y="c1_3"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 4]`) %>% as_tibble() %>% mutate(x="prs", y="c1_4"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 5]`) %>% as_tibble() %>% mutate(x="prs", y="c1_5"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 6]`) %>% as_tibble() %>% mutate(x="prs", y="c1_6"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 7]`) %>% as_tibble() %>% mutate(x="prs", y="c1_7"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 8]`) %>% as_tibble() %>% mutate(x="prs", y="c1_8"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 9]`) %>% as_tibble() %>% mutate(x="prs", y="c1_9"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 10]`) %>% as_tibble() %>% mutate(x="prs", y="c1_10"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 11]`) %>% as_tibble() %>% mutate(x="prs", y="c1_11"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 12]`) %>% as_tibble() %>% mutate(x="prs", y="c1_12"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 13]`) %>% as_tibble() %>% mutate(x="prs", y="c1_13"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 14]`) %>% as_tibble() %>% mutate(x="prs", y="c1_14"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 15]`) %>% as_tibble() %>% mutate(x="prs", y="c1_15"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 16]`) %>% as_tibble() %>% mutate(x="prs", y="c1_16"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 17]`) %>% as_tibble() %>% mutate(x="prs", y="c1_17"),
    fast_assoc(y=obs$prs, x=obs$r1) %>% as_tibble() %>% mutate(x="prs", y="r1"),
  )  %>%
    mutate(study="observational")

  prot <- phen[phen$d==0,][1:protein_gwas, ] #add element of these individuals not having disease 
  res_protein <- bind_rows(
    fast_assoc(y=prot$d, x=prot$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=prot$`c0[, 1]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_1", y="d"),
    fast_assoc(y=prot$`c0[, 2]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_2", y="d"),
    fast_assoc(y=prot$`c0[, 3]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_3", y="d"),
    fast_assoc(y=prot$`c0[, 4]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_4", y="d"),
    fast_assoc(y=prot$`c0[, 5]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_5", y="d"),
    fast_assoc(y=prot$`c0[, 6]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_6", y="d"),
    fast_assoc(y=prot$`c0[, 7]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_7", y="d"),
    fast_assoc(y=prot$`c0[, 8]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_8", y="d"),
    fast_assoc(y=prot$`c0[, 9]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_9", y="d"),
    fast_assoc(y=prot$`c0[, 10]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_10", y="d"),
    fast_assoc(y=prot$`c0[, 11]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_11", y="d"),
    fast_assoc(y=prot$`c0[, 12]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_12", y="d"),
    fast_assoc(y=prot$`c0[, 13]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_13", y="d"),
    fast_assoc(y=prot$`c0[, 14]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_14", y="d"),
    fast_assoc(y=prot$`c0[, 15]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_15", y="d"),
    fast_assoc(y=prot$`c0[, 16]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_16", y="d"),
    fast_assoc(y=prot$`c0[, 17]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_17", y="d"),
    fast_assoc(y=prot$r0, x=prot$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c0_1"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 2]`) %>% as_tibble() %>% mutate(x="prs", y="c0_2"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 3]`) %>% as_tibble() %>% mutate(x="prs", y="c0_3"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 4]`) %>% as_tibble() %>% mutate(x="prs", y="c0_4"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 5]`) %>% as_tibble() %>% mutate(x="prs", y="c0_5"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 6]`) %>% as_tibble() %>% mutate(x="prs", y="c0_6"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 7]`) %>% as_tibble() %>% mutate(x="prs", y="c0_7"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 8]`) %>% as_tibble() %>% mutate(x="prs", y="c0_8"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 9]`) %>% as_tibble() %>% mutate(x="prs", y="c0_9"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 10]`) %>% as_tibble() %>% mutate(x="prs", y="c0_10"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 11]`) %>% as_tibble() %>% mutate(x="prs", y="c0_11"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 12]`) %>% as_tibble() %>% mutate(x="prs", y="c0_12"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 13]`) %>% as_tibble() %>% mutate(x="prs", y="c0_13"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 14]`) %>% as_tibble() %>% mutate(x="prs", y="c0_14"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 15]`) %>% as_tibble() %>% mutate(x="prs", y="c0_15"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 16]`) %>% as_tibble() %>% mutate(x="prs", y="c0_16"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 17]`) %>% as_tibble() %>% mutate(x="prs", y="c0_17"),
    fast_assoc(y=prot$prs, x=prot$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=prot$`c1[, 1]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_1", y="d"),
    fast_assoc(y=prot$`c1[, 2]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_2", y="d"),
    fast_assoc(y=prot$`c1[, 3]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_3", y="d"),
    fast_assoc(y=prot$`c1[, 4]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_4", y="d"),
    fast_assoc(y=prot$`c1[, 5]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_5", y="d"),
    fast_assoc(y=prot$`c1[, 6]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_6", y="d"),
    fast_assoc(y=prot$`c1[, 7]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_7", y="d"),
    fast_assoc(y=prot$`c1[, 8]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_8", y="d"),
    fast_assoc(y=prot$`c1[, 9]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_9", y="d"),
    fast_assoc(y=prot$`c1[, 10]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_10", y="d"),
    fast_assoc(y=prot$`c1[, 11]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_11", y="d"),
    fast_assoc(y=prot$`c1[, 12]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_12", y="d"),
    fast_assoc(y=prot$`c1[, 13]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_13", y="d"),
    fast_assoc(y=prot$`c1[, 14]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_14", y="d"),
    fast_assoc(y=prot$`c1[, 15]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_15", y="d"),
    fast_assoc(y=prot$`c1[, 16]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_16", y="d"),
    fast_assoc(y=prot$`c1[, 17]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_17", y="d"),
    fast_assoc(y=prot$r1, x=prot$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c1_1"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 2]`) %>% as_tibble() %>% mutate(x="prs", y="c1_2"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 3]`) %>% as_tibble() %>% mutate(x="prs", y="c1_3"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 4]`) %>% as_tibble() %>% mutate(x="prs", y="c1_4"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 5]`) %>% as_tibble() %>% mutate(x="prs", y="c1_5"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 6]`) %>% as_tibble() %>% mutate(x="prs", y="c1_6"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 7]`) %>% as_tibble() %>% mutate(x="prs", y="c1_7"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 8]`) %>% as_tibble() %>% mutate(x="prs", y="c1_8"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 9]`) %>% as_tibble() %>% mutate(x="prs", y="c1_9"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 10]`) %>% as_tibble() %>% mutate(x="prs", y="c1_10"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 11]`) %>% as_tibble() %>% mutate(x="prs", y="c1_11"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 12]`) %>% as_tibble() %>% mutate(x="prs", y="c1_12"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 13]`) %>% as_tibble() %>% mutate(x="prs", y="c1_13"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 14]`) %>% as_tibble() %>% mutate(x="prs", y="c1_14"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 15]`) %>% as_tibble() %>% mutate(x="prs", y="c1_15"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 16]`) %>% as_tibble() %>% mutate(x="prs", y="c1_16"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 17]`) %>% as_tibble() %>% mutate(x="prs", y="c1_17"),
    fast_assoc(y=prot$prs, x=prot$r1) %>% as_tibble() %>% mutate(x="prs", y="r1"),
  ) %>%
    mutate(study="protein_gwas")
  out_all <- paste0(rev_mr_result[row,]$prot, "_res_all.csv")
  write.csv(res_all, file = out_all, row.names = FALSE, quote = FALSE)
  out_obs <- paste0(rev_mr_result[row,]$prot, "_res_obs.csv")
  write.csv(res_obs, file = out_obs, row.names = FALSE, quote = FALSE)
  out_protein <- paste0(rev_mr_result[row,]$prot, "_res_protein.csv")
  write.csv(res_protein, file = out_protein, row.names = FALSE, quote = FALSE)
  return(bind_rows(res_all, res_obs, res_protein))
} 

#reverse MR approach using prs to find non causal biomarkers
#IVW and two-stage least squares (should be the same value)
ivw_analysis <- function(dat){
  phen <- dat$phen
  gen <- dat$geno
mr_dat <- get_effs(phen$prs, phen$r0, gen) #might need to change this genetic matrix
print("Reverse MR")
mr(mr_dat, metho="mr_ivw") %>% str()
mr_dat <- get_effs(phen$c1, phen$d, gen) #might need to change this genetic matrix
print("Forward MR")
mr(mr_dat, metho="mr_ivw") %>% str()
# MR using fixed effects IVW [two-stage least squares method] use system fit package
print("2SLS")
summary(systemfit::systemfit(phen$prs ~ phen$r0, method="2SLS", inst = ~gen)) # 2SLS #might need to change this genetic matrix ##issue with this code as getting NAs
print("Observational")
summary(lm(phen$prs ~ phen$r0)) #confounded observational estimate
}
```

Glioma

Glioblastoma

Non-glioblastoma
