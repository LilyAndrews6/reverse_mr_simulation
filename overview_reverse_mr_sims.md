overview\_reverse\_mr\_sims
================
Lily Andrews
7/22/2022

Data generating model with a degree of temporality where c0 refers to
causal biomarker at t=0 (pre-onset/general population), c1 refers to
causal biomarker at t=1 (post-onset/case-control):

``` r
#' Simulation model of DAG for the role of causal and non-causal biomarkers in relation to liability to disease
#'
#' @param nid number of individuals
#' @param nsnp number of SNPs
#' @param gc_maf causal biomarker minor allele frequency input
#' @param b_gl genetic liability from discovered SNPs
#' @param rsq_prs h-squared of y explained
#' @param rsq_zl H-squared of y when combined with PRS
#' @param d_prev disease prevalence
#' @param b_lr0 beta value of disease liability to non-causal biomarker at pre-onset
#' @param b_c0l beta value of causal biomarker at pre-onset to disease liability
#' @param gr_maf non-causal biomarker minor allele frequency input
#' @param rsq_gr0 r-squared value for the non-causal pQTL 
#' @param rsq_gc0 r-squared value for the causal pQTL
#' @param b_u1r1 beta value of unmeasured confounder 1 to non-causal biomarker at post-onset
#' @param b_u1l beta value of unmeasured confounder 1 to disease liability
#' @param b_u2c1 beta value of unmeasured confounder 2 to causal biomarker at post-onset
#' @param b_u2l beta value of unmeasured confounder 2 to disease liability
#' @param b_dr1 beta value of disease to non-causal biomarker at post-onset
#' @param b_dc1 beta value of disease to causal biomarker at post-onset
#' @param b_c0c1 beta value of causal biomarker at pre-onset to causal biomarker at post-onset
#' @param b_r0r1 beta value of non-causal biomarker at pre-onset to non-causal biomarker at post-onset
#' @param b_u2c0 beta value of unmeasured confounder 2 to causal biomarker at pre-onset 
#' @param b_u1r0 beta value of unmeasured confounder 1 to non-causal biomarker at pre-onset
#' @param b_gcc0 
#' @param b_lc0 beta value of disease liability to causal biomarker at pre-onset 
#' @param b_u1d beta value of unmeasured confounder 1 to disease
#' @param b_u2d beta value of unmeasured confounder 2 to disease
#'
#' @return
#' @export
#'
#' @examples
dgmodel <- function(nid, nsnp, gc_maf, b_gl, rsq_prs, rsq_zl, d_prev, b_lr0, b_c0l, gr_maf, rsq_gr0, rsq_gc0, b_u1r1, b_u1l, b_u2c1, b_u2l, b_dr1, b_dc1, b_c0c1, b_r0r1, b_u2c0, b_u1r0, b_gcc0, b_lc0, b_u1d, b_u2d)
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

List of expected variances from the data generating model: v\_u1 \~
N(0,1) v\_u2 = 1 v\_l = 1 v\_pc = 1 v\_pr = 1 v\_z = 1 v\_d =
Binomial(n=1, p=d\_prev \* (1-d\_prev))

Check output of data generating model

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

Creating simulated data to put into the data generating model (currently
based on protein data from applied study)

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
    ##  1 u1     1.00    -3.84e- 3
    ##  2 u2     0.994   -5.48e- 4
    ##  3 r0     0.931   -3.36e- 4
    ##  4 r1     0.936    3.58e- 3
    ##  5 c0     0.907   -3.23e- 3
    ##  6 c1     1.00     6.22e- 3
    ##  7 prs    1        1.07e-17
    ##  8 prs_w  1        4.65e-17
    ##  9 l      0.271   -2.41e- 3
    ## 10 prob_l 0.00527  9.75e- 1
    ## 11 d      0.0240   2.46e- 2
    ##           [,1]
    ## [1,] 0.1886057

Analysis of data generating model looking at all, case control,
observational and protein gwas

``` r
library(tidyverse)
```

    ## ── Attaching packages ─────────────────────────────────────── tidyverse 1.3.1 ──

    ## ✔ ggplot2 3.3.6     ✔ dplyr   1.0.9
    ## ✔ tidyr   1.2.0     ✔ stringr 1.4.0
    ## ✔ readr   2.1.2     ✔ forcats 0.5.1
    ## ✔ purrr   0.3.4

    ## Warning: package 'tidyr' was built under R version 4.0.5

    ## Warning: package 'readr' was built under R version 4.0.5

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()

``` r
dgmodel_analysis <- function(dat, ncase=1000, ncontrol=1000, protein_gwas=10000)
{
  d <- dat$phen
  res_all <- bind_rows(
    fast_assoc(y=d$d, x=d$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
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
  obs <- bind_rows(d[d$d == 1,][1:ncontrol,], d[d$d == 0,][1:ncase,])
  
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

Test how well the model runs in all, case control, observational and
protein gwas

``` r
library(dplyr)
dgmodel_analysis(dat)
```

    ## # A tibble: 27 × 9
    ##           ahat     bhat       se     fval      pval      n x     y     study    
    ##          <dbl>    <dbl>    <dbl>    <dbl>     <dbl>  <int> <chr> <chr> <chr>    
    ##  1  0.0246     -0.0195  0.000486 1617.    0         100000 prs   d     all      
    ##  2  0.00574    -0.364   0.0194    353.    1.40e- 78 100000 c0    d     all      
    ##  3  0.00830    -0.351   0.0197    319.    3.62e- 71 100000 r0    d     all      
    ##  4 -0.00000790 -0.00244 0.00332     0.541 4.62e-  1 100000 prs   c0    all      
    ##  5  0.0000257   0.0765  0.00327   549.    5.23e-121 100000 prs   r0    all      
    ##  6  0.00558     0.0259  0.0204      1.61  2.05e-  1 100000 c1    d     all      
    ##  7  0.00244     0.0463  0.0197      5.50  1.90e-  2 100000 r1    d     all      
    ##  8  0.0000600  -0.00965 0.00316     9.35  2.23e-  3 100000 prs   c1    all      
    ##  9 -0.0000233   0.00651 0.00327     3.97  4.62e-  2 100000 prs   r1    all      
    ## 10  0.430      -0.180   0.00974   342.    1.10e- 70   2000 prs   d     observat…
    ## # … with 17 more rows

Scenerio 1 - Alter heritability explained

``` r
set.seed(100)
reps = 1000
scenerio <- 1
dat <- dgmodel(
  nid=100000, #number of individuals
  nsnp=99, #number of snps
  gc_maf=runif(99, 0.05, 0.95), #no rare variants included
  b_gl=rnorm(99), 
  rsq_prs=0.06, 
  rsq_zl=0.19, 
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
```

Scenerio 2 - Alteration of causal SNPs

``` r
set.seed(13)
reps = 1000
scenerio <- 2
dat <- dgmodel(
  nid=100000, #number of individuals
  nsnp=99, #number of snps
  gc_maf=runif(99, 0.05, 0.95), #no rare variants included
  b_gl=rnorm(99), 
  rsq_prs=0.06, 
  rsq_zl=0.19, 
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
```
