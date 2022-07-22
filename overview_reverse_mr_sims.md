overview\_reverse\_mr\_sims
================
Lily Andrews
7/22/2022

\#Need to add degree to temporality to the model c0 refers to causal
biomarker at t=0 (pre-onset/general population), c1 refers to causal
biomarker at t=1 (post-onset/case-control) \#r0 refers to non-causal
biomarker at t=0 (pre-onset/general population), r1 refers to non-causal
biomarker at t=1 (post-onset/case-control) \#list of expected variances
v\_u1 \~ N(0,1) v\_u2 = 1 v\_l = 1 v\_pc = 1 v\_pr = 1 v\_z = 1 v\_d =
Binomial(n=1, p=d\_prev \* (1-d\_prev))

``` r
#' Title
#'
#' @param nid number of individuals
#' @param nsnp 
#' @param gc_maf number of SNPs
#' @param b_gl 
#' @param rsq_prs 
#' @param rsq_zl 
#' @param d_prev 
#' @param b_lr0 
#' @param b_c0l 
#' @param gr_maf 
#' @param rsq_gr0 
#' @param rsq_gc0 
#' @param b_u1r1 
#' @param b_u1l 
#' @param b_u2c1 
#' @param b_u2l 
#' @param b_dr1 
#' @param b_dc1 
#' @param b_c0c1 
#' @param b_r0r1 
#' @param b_u2c0 
#' @param b_u1r0 
#' @param b_g1c0 
#' @param b_gcc0 
#' @param b_lc0 
#' @param b_u1d 
#' @param b_u2d 
#'
#' @return
#' @export
#'
#' @examples
dgmodel <- function(nid, nsnp, gc_maf, b_gl, rsq_prs, rsq_zl, d_prev, b_lr0, b_c0l, gr_maf, rsq_gr0, rsq_gc0, b_u1r1, b_u1l, b_u2c1, b_u2l, b_dr1, b_dc1, b_c0c1, b_r0r1, b_u2c0, b_u1r0, b_g1c0, b_gcc0, b_lc0, b_u1d, b_u2d)
{
  u1 <- rnorm(nid) #normal distribution of unmeasured confounder
  u2 <- rnorm(nid) #normal distribution of unmeasured confounder
  gc <- make_geno(nid, nsnp, gc_maf) #create genotype matrix for causal variant
  gr <- rbinom(nid, 2, gr_maf) #creating genotype matrix for non causal variant 
  # TODO need to fix this
  c0 <- scale(gc[,1]) * b_gcc0 + u2 * b_u2c0 + rnorm(nid, sd=sqrt(1-rsq_gc0-b_lc0^2-b_u2c0^2)) #we hypothesised one SNP down this route? no need for unmeasured confounder in this case? rnorm introduces the variance into the equation 
  
  # TODO need to fix this
  rsq_prs <- rsq_prs - (0.1)^2 # would 0.1 be b_prsl
  prs <- scale(gc %*% b_gl) #b_gl does this mean beta genetic liability or route PRS to disease liability, why do we include this line if prs_w is already used
  prs_w <- scale(gc[,-1] %*% b_gl[-1]) #had to remove the SNP as the liability calculated later on would include an extra SNP
  z <- rnorm(nid) #included the rest of the known SNPs
  l <- prs_w * sqrt(rsq_prs) + z * sqrt(rsq_zl) + c0 * b_c0l + u1 * b_u1l + u2 * b_u2l #total genetic liability to disease, decided not to add error into liability to avoid diagram c from happening https://www.nature.com/articles/nrg3377/figures/1 
  r0 <- scale(gr) * sqrt(rsq_gr0) + l * b_lr0 + u1 * b_u1r0 + rnorm(nid, sd=sqrt(1-rsq_gr0-b_lr0^2-b_u1r0^2)) # instead of sqrt(rsq_gr0) could this be b_grc0
  # generate probability of disease
  prob_l <- simulateGP::gx_to_gp(gx=scale(l), h2x=rsq_prs + rsq_zl + b_c0l^2 + b_u1d^2 + b_u2d^2, prev = d_prev) #translate disease risk from liability to probability scale would h2x be prs_w instead?
  
  # generate random disease outcome
  # switch this 0/1 or fix simualteGP function
  d <- rbinom(nid, 1, prob_l) #case as 0 control is 1
  d <- abs(d-1) #this function flips case and control which makes case 1 and control 0
  #response to disease state
  c1 <- scale(gc[,1]) * 0.1 + u2 * b_u2c1 + c0 *b_c0c1 + d * b_dc1 + rnorm(nid, sd=sqrt(1-(0.1)^2 - b_u2c1^2- b_c0c1^2)) #we hypothesised one SNP down this route? would 0.1 be sqrt(rsq_gc0)
  
  r1 <- scale(gr) * sqrt(rsq_gr0) + r0 * b_r0r1 + u1 * b_u1r1 + d * b_dr1 + rnorm(nid, sd=sqrt(1-rsq_gc0-b_lr0^2-b_u1r1^2)) 

  phen <- tibble(
    u1, u2, r0, r1, c0, c1, prs, prs_w, l, prob_l, d
  )
  return(list(geno=gc, phen=phen)) #check genotype would be Gc when it was G
}
```

\#think about the one snp of PC variable

``` r
dgmodel_check <- function(dat)
{
  print(tibble(
    col = names(dat$phen),
    vars = apply(dat$phen, 2, var),
    means = colMeans(dat$phen)
  ))
  print(cor(dat$phen$prs, dat$phen$l)^2)
  
}
```

``` r
library(tibble)
library(simulateGP)
dat <- dgmodel(
  nid=100000, #number of individuals
  nsnp=99, #number of snps
  gc_maf=runif(99, 0.05, 0.95), #no rare variants included
  b_gl=rnorm(99), 
  rsq_prs=0.06, #h-squared value
  rsq_zl=0.19, #to total 0.25 of total H-squared
  d_prev=5/1000, 
  b_lr0=sqrt(0.1), 
  b_c0l=sqrt(0.01), 
  gr_maf=0.5, 
  rsq_gr0=0.1,   # this is a pqtl - often large effect size
  rsq_gc0=0.1,   # this is a pqtl - often large effect size 
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

dgmodel_check(dat)
```

    ## # A tibble: 11 × 3
    ##    col       vars     means
    ##    <chr>    <dbl>     <dbl>
    ##  1 u1     1.00    -2.01e- 3
    ##  2 u2     0.998   -5.70e- 3
    ##  3 r0     0.929    5.91e- 3
    ##  4 r1     0.940    1.35e- 3
    ##  5 c0     0.894    2.71e- 4
    ##  6 c1     1.01     2.51e- 3
    ##  7 prs    1        4.78e-18
    ##  8 prs_w  1        1.60e-17
    ##  9 l      0.272   -8.24e- 4
    ## 10 prob_l 0.00508  9.75e- 1
    ## 11 d      0.0241   2.47e- 2
    ##           [,1]
    ## [1,] 0.1795225
