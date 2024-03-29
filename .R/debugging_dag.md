debugging_dag
================
2024-01-31

Install packages

``` r
library(tibble)
library(simulateGP)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(TwoSampleMR)
```

    ## TwoSampleMR version 0.5.8 
    ## [>] New: Option to use non-European LD reference panels for clumping etc
    ## [>] Some studies temporarily quarantined to verify effect allele
    ## [>] See news(package='TwoSampleMR') and https://gwas.mrcieu.ac.uk for further details

    ## 
    ## Attaching package: 'TwoSampleMR'

    ## The following objects are masked from 'package:simulateGP':
    ## 
    ##     allele_frequency, contingency, get_population_allele_frequency

``` r
dgmodel <- function(nid, #3301 for protein and 30686 for GWAS
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
    0.26, 0.54
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
    u1, u2, r0, r1, c0[,1], c1[,1], prs, l, prob_l, d
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
    fast_assoc(y=phen$r0, x=phen$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=phen$prs, x=phen$`c0[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c0_1"),
    fast_assoc(y=phen$prs, x=phen$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=phen$`c1[, 1]`, x=phen$d) %>% as_tibble() %>% mutate(x="c1_1", y="d"),
    fast_assoc(y=phen$r1, x=phen$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=phen$prs, x=phen$`c1[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c1_1"),
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
    fast_assoc(y=obs$r0, x=obs$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=obs$prs, x=obs$`c0[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c0_1"),
    fast_assoc(y=obs$prs, x=obs$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=obs$`c1[, 1]`, x=obs$d) %>% as_tibble() %>% mutate(x="c1_1", y="d"),
    fast_assoc(y=obs$r1, x=obs$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=obs$prs, x=obs$`c1[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c1_1"),
    fast_assoc(y=obs$prs, x=obs$r1) %>% as_tibble() %>% mutate(x="prs", y="r1"),
  )  %>%
    mutate(study="observational")

  prot <- phen[phen$d==0,][1:protein_gwas, ] #add element of these individuals not having disease 
  res_protein <- bind_rows(
    fast_assoc(y=prot$d, x=prot$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=prot$`c0[, 1]`, x=prot$d) %>% as_tibble() %>% mutate(x="c0_1", y="d"),
    fast_assoc(y=prot$r0, x=prot$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=prot$prs, x=prot$`c0[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c0_1"),
    fast_assoc(y=prot$prs, x=prot$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=prot$`c1[, 1]`, x=prot$d) %>% as_tibble() %>% mutate(x="c1_1", y="d"),
    fast_assoc(y=prot$r1, x=prot$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=prot$prs, x=prot$`c1[, 1]`) %>% as_tibble() %>% mutate(x="prs", y="c1_1"),
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
mr(mr_dat, metho=c("mr_ivw", "mr_wald_ratio")) %>% str()
mr_dat <- get_effs(phen$`c1[, 1]`, phen$d, gen) #might need to change this genetic matrix
print("Forward MR")
mr(mr_dat, metho=c("mr_ivw", "mr_wald_ratio")) %>% str()
# MR using fixed effects IVW [two-stage least squares method] use system fit package
print("2SLS")
summary(systemfit::systemfit(phen$prs ~ phen$r0, method="2SLS", inst = ~gen)) # 2SLS #might need to change this genetic matrix ##issue with this code as getting NAs
print("Observational")
summary(lm(phen$prs ~ phen$r0)) #confounded observational estimate
}
```

``` r
rev_mr_result <- tribble(
  ~beta, ~prot, ~rsq_lr0,
 0.0880, "HMGCS1", 0.2040
)

for(row in 1:nrow(rev_mr_result)){
dat <- dgmodel(nid=10000, #making large sample size so enough individuals for case control analysis
  h2_known=0.06, #h-squared value (form of r-squared) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667278/
  h2_unknown=0.19, #to total 0.25 of total H-squared (form of r-squared) https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667278/
  d_prev=6/100000, #glioma prevalence https://www.ncbi.nlm.nih.gov/books/NBK441874/#:~:text=In%20the%20United%20States%2C%20there,the%20least%20malignant%20brain%20tumors.
  #b_lr0=rev_mr_result[row,]$beta,  #adding rev MR results into data generating model
  b_c0l=0, 
  gc_maf=0, #don't think this is needed anymore
  gr_maf=0.20, #need to estimate this allele frequency
  rsq_gr0=0,   
  rsq_gc0=0, 
  rsq_lr0=0,
  b_u1r1=0, 
  b_u1l=0, 
  b_u2c1=0, 
  b_u2l=0, 
  b_dr1=0, 
  b_dc1=0, 
  b_c0c1=0, 
  b_r0r1=0, 
  b_u2c0=0, 
  b_gcc0=1, 
  b_lc0=0,
  b_u1r0=0, 
  b_u1d=0, 
  b_u2d=0, 
  b_gcc1=0,
  var_l=1, #estimated variance of liability to be 1
  vp=1) #estimated variance of protein to be 1

dgmodel_check(dat)

dgmodel_analysis(dat, ncase=12496, ncontrol=18190, protein_gwas = 3301, rev_mr_result)}
```

    ## # A tibble: 10 × 3
    ##    col          vars     means
    ##    <chr>       <dbl>     <dbl>
    ##  1 u1      1.00       0.000576
    ##  2 u2      0.967     -0.0178  
    ##  3 r0      1.00      -0.0196  
    ##  4 r1      0.986     -0.00665 
    ##  5 c0[, 1] 0.484      1.08    
    ##  6 c1[, 1] 1.00       0.000883
    ##  7 prs     0.0602    -0.00230 
    ##  8 l       0.0357     0.281   
    ##  9 prob_l  0.0000327  0.998   
    ## 10 d       0.00200    0.002   
    ##           [,1]
    ## [1,] 0.0835587

``` r
ivw_analysis(dat)
```

    ## [1] "Reverse MR"

    ## Analysing 'X' on 'Y'

    ## 'data.frame':    1 obs. of  9 variables:
    ##  $ id.exposure: chr "X"
    ##  $ id.outcome : chr "Y"
    ##  $ outcome    : chr "Y"
    ##  $ exposure   : chr "X"
    ##  $ method     : chr "Wald ratio"
    ##  $ nsnp       : num 1
    ##  $ b          : num -44.1
    ##  $ se         : num 32
    ##  $ pval       : num 0.168
    ## [1] "Forward MR"

    ## Analysing 'X' on 'Y'

    ## 'data.frame':    1 obs. of  9 variables:
    ##  $ id.exposure: chr "X"
    ##  $ id.outcome : chr "Y"
    ##  $ outcome    : chr "Y"
    ##  $ exposure   : chr "X"
    ##  $ method     : chr "Wald ratio"
    ##  $ nsnp       : num 1
    ##  $ b          : num 0.247
    ##  $ se         : num 0.0354
    ##  $ pval       : num 2.89e-12
    ## [1] "2SLS"
    ## [1] "Observational"

    ## 
    ## Call:
    ## lm(formula = phen$prs ~ phen$r0)
    ## 
    ## Residuals:
    ##      Min       1Q   Median       3Q      Max 
    ## -0.92358 -0.16614  0.00269  0.16235  0.95189 
    ## 
    ## Coefficients:
    ##               Estimate Std. Error t value Pr(>|t|)
    ## (Intercept) -2.305e-03  2.454e-03  -0.939    0.348
    ## phen$r0     -2.103e-05  2.453e-03  -0.009    0.993
    ## 
    ## Residual standard error: 0.2454 on 9998 degrees of freedom
    ## Multiple R-squared:  7.35e-09,   Adjusted R-squared:  -0.0001 
    ## F-statistic: 7.349e-05 on 1 and 9998 DF,  p-value: 0.9932
