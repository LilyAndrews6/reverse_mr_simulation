Sims Functions
================
Lily Andrews
2023-11-28

Install packages

``` r
library(tibble)
library(simulateGP)
library(dplyr)
library(TwoSampleMR)
library(parallel)
library(data.table)
library(ggplot2)
```

Converting the liability score (gx) to the probability scale (gp)

``` r
gx_to_gp <- function(gx, h2x, prev){
  x_prime <- qnorm(prev, 0, 1)
  p <- pnorm(x_prime, mean=-gx, sd=sqrt(1-h2x), lower.tail = TRUE)
  return(p)
}
```

Logistic association calculations

``` r
logistic_assoc <- function(y, x)
{
    logistic_model <- summary(glm(y ~ x, family="binomial"))$coefficients
    return(list(
        ahat=logistic_model[1,1],
        bhat=logistic_model[2,1],
        se=logistic_model[2,2],
        fval=logistic_model[2,3]^2,
        pval=logistic_model[2,4],
        n=length(y)
    ))
}
```

Create data generating model (dgmodel) - we are assuming a variance of 1
in the variables generated

``` r
#' Simulation model of DAG for the role of causal and non-causal biomarkers in relation to liability to disease
#'
#' @param nid number of individuals
#' @param b_gcc0 beta of causal genotype matrix to causal protein (pre diagnosis)
#' @param b_u2c0 beta of unmeasured confounder 2 to causal protein (pre diagnosis)
#' @param gr_maf non-causal biomarker minor allele frequency input
#' @param nsnp number of non-causal SNPs
#' @param h2_known known h-squared value https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667278/
#' @param h2_unknown total H-squared https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4667278/
#' @param b_c0l beta of causal protein (pre-diagnosis) to disease liability
#' @param b_u2l beta of unmeasured confounder 2 to disease liability
#' @param rsq_gr0 r-squared value for the non-causal protein
#' @param rsq_lr0 r-squared value for disease liability to the non-causal protein
#' @param b_u1r0 beta of unmeasured confounder 1 to non-causal protein (pre diagnosis)
#' @param b_u1d ##would this be b_u1l*b_ld
#' @param b_u2d ##would this be b_u2l*b_ld
#' @param b_c0c1 beta of causal biomarker (pre diagnosis) to causal biomarker (post diagnosis)
#' @param b_u2c1 beta of unmeasured confounder 2 to causal biomarker (post diagnosis)
#' @param b_dc1 beta of disease to causal biomarker (post diagnosis)
#' @param b_r0r1 beta of non-causal protein (pre diagnosis) to non-causal protein (post diagnosis)
#' @param b_u1r1 beta of unmeasured confounder 1 to non-causal protein (post-diagnosis)
#' @param b_dr1 beta of disease to non-causal protein (post diagnosis)
#' @param b_dx1 beta of disease to protein as a consequence of disease
#'
#' @return
#' @export
#'
#' @examples

dgmodel <- function(nid, 
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
b_u1d,
b_u2d,
d_prev,
b_c0c1,
b_u2c1,
b_dc1,
b_r0r1,
b_u1r1,
b_dr1,
b_dx1
)  #assuming variance of genetic liability is 1
{
u1 <- rnorm(nid) #normal distribution of unmeasured confounder
  
u2 <- rnorm(nid) #normal distribution of unmeasured confounder

gwashits <- tribble(
    ~beta, ~af,
    0.26, 0.54
) #known causal variants from glioma gwas data

gc <- sapply(gwashits$af, \(x) rbinom(nid, 2, x)) #create genotype matrix for causal variant

c0 <- gc * b_gcc0 + u2*b_u2c0 + rnorm(nid, sd=sqrt(1-b_gcc0^2-b_u2c0^2))

#gr_maf <- runif(nsnp, 0.01, 0.99)
gprs_maf <- rep(gr_maf, nsnp) #create minor allele frequency for each snp

prs_b <- rnorm(nsnp) #randomly generated betas for prs

gprs <- make_geno(nid, nsnp, gprs_maf) #create genotype matrix for prs

prs_known <- (gprs %*% prs_b) %>% {scale(.) * sqrt(h2_known)} #generate known prs

prs_unknown <- rnorm(nid, 0, sd = sqrt(h2_unknown)) #generate unknown prs

prs <- prs_known + prs_unknown #generate total prs

rsq_prs <- h2_known - (0.1)^2 #generate r-squared prs

l <- c0 * b_c0l + scale(prs) * sqrt(rsq_prs)+ u2 * b_u2l + rnorm(nid, 0, sd = 1 - sqrt(rsq_prs + b_c0l^2 + b_u2l^2)) #total genetic liability to disease

gr <- rbinom(nid, 2, gr_maf) #creating genotype matrix for non causal variant

r0 <- scale(gr) * sqrt(rsq_gr0) + scale(l) * sqrt(rsq_lr0) + u1 * b_u1r0 + rnorm(nid, sd=sqrt(1 - rsq_gr0 - rsq_lr0 - b_u1r0^2)) #create non causal variable

prob_l <- gx_to_gp(gx=scale(l), h2x=rsq_prs + b_c0l^2 + b_u1d^2 + b_u2d^2, prev = d_prev) #generating the probability of liability to disease

d <- rbinom(nid, 1, prob_l) #generate disease in individuals based on liability

c1 <- c0 * b_c0c1 + u2 * b_u2c1 + d * b_dc1 #generate causal variant post disease

r1 <- r0 * b_r0r1 + u1 * b_u1r1 + d * b_dr1 #generate non-causal variant post disease

x1 <- d * b_dx1 + rnorm(nid, sd=sqrt(1 - b_dx1^2*var(d)))  #protein as a consequence of disease - error needs addition of variance of disease

join_mat <- cbind(gprs, gr, gc) #combine genotypes

phen <- tibble(
  u1, u2, r0, r1, c0, c1, prs, l, prob_l, d, x1
)  
return(list(geno_gc=gc,geno_prs=gprs, geno_gr=gr, phen=phen, geno_join=join_mat, af=gprs_maf))} 
```

Check output of data generating model (dgmodel_check)

``` r
dgmodel_check <- function(dat)
{
  print(tibble(
    col = names(dat$phen),
    vars = apply(dat$phen, 2, var), 
    means = colMeans(dat$phen)
  ), n=25)
  print(cor(dat$phen$prs, dat$phen$l)^2)
}
```

Assess observational, case control and protein GWAS data
(dgmodel_analysis)

``` r
#looking at different scenerios, all, case control, protein gwas
dgmodel_analysis <- function(dat, ncase, ncontrol, protein_gwas, rev_mr_result, row, iteration)
{
  phen <- dat$phen
  res_all <- bind_rows( 
    fast_assoc(y=phen$d, x=phen$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=phen$c0, x=phen$d) %>% as_tibble() %>% mutate(x="c0", y="d"),
    fast_assoc(y=phen$r0, x=phen$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=phen$prs, x=phen$c0) %>% as_tibble() %>% mutate(x="prs", y="c0"),
    fast_assoc(y=phen$prs, x=phen$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=phen$c1, x=phen$d) %>% as_tibble() %>% mutate(x="c1", y="d"),
    fast_assoc(y=phen$r1, x=phen$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=phen$prs, x=phen$c1) %>% as_tibble() %>% mutate(x="prs", y="c1"),
    fast_assoc(y=phen$prs, x=phen$c1) %>% as_tibble() %>% mutate(x="prs", y="r1")
  ) %>%
    mutate(study="all")
  
  obs <- bind_rows(phen[phen$d == 0,][1:ncontrol,], phen[phen$d == 1,][1:ncase,])  
  case <- length(subset(obs$d, obs$d==1))
  print("Number of cases: ") 
  print(case)
  control <- length(subset(obs$d, obs$d==0))
  print("Number of controls: ")
  print(control)
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

  prot <- phen[phen$d==0,][1:protein_gwas, ] 
  res_protein <- bind_rows(
    fast_assoc(y=prot$d, x=prot$prs) %>% as_tibble() %>% mutate(x="prs", y="d"),
    fast_assoc(y=prot$c0, x=prot$d) %>% as_tibble() %>% mutate(x="c0", y="d"),
    fast_assoc(y=prot$r0, x=prot$d) %>% as_tibble() %>% mutate(x="r0", y="d"),
    fast_assoc(y=prot$prs, x=prot$c0) %>% as_tibble() %>% mutate(x="prs", y="c0"),
    fast_assoc(y=prot$prs, x=prot$r0) %>% as_tibble() %>% mutate(x="prs", y="r0"),
    fast_assoc(y=prot$c1, x=prot$d) %>% as_tibble() %>% mutate(x="c1", y="d"),
    fast_assoc(y=prot$r1, x=prot$d) %>% as_tibble() %>% mutate(x="r1", y="d"),
    fast_assoc(y=prot$prs, x=prot$c1) %>% as_tibble() %>% mutate(x="prs", y="c1"),
    fast_assoc(y=prot$prs, x=prot$r1) %>% as_tibble() %>% mutate(x="prs", y="r1"),
  ) %>%
    mutate(study="protein_gwas")
  out_all <- paste0(rev_mr_result[row,]$prot, iteration, "_res_all.csv")
  write.csv(res_all, file = out_all, row.names = FALSE, quote = FALSE)
  out_obs <- paste0(rev_mr_result[row,]$prot, iteration, "_res_obs.csv")
  write.csv(res_obs, file = out_obs, row.names = FALSE, quote = FALSE)
  out_protein <- paste0(rev_mr_result[row,]$prot, iteration,"_res_protein.csv")
  write.csv(res_protein, file = out_protein, row.names = FALSE, quote = FALSE)
  return(bind_rows(obs%>%
    mutate(study="observational"), prot%>%
    mutate(study="protein_gwas")))
} 
```

Generaritng two-sample reverse MR estimate

``` r
rev_mr_prot <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_prs[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$r0
  g <- geno_mat
  gwasy <- gwas(y,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_prs[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$d #change this to d and do logistic regression to check this as it's binary, l for fast assoc
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  prs_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)#select SNPs reaching threshold
  print(prs_snp)
  if (dim(prs_snp)[1]>0){
   mr_res <- mr(prs_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
  all_mr_res <- mr(harmonise, metho=c("mr_ivw", "mr_wald_ratio"))
   return(list(pval=mr_res, all=all_mr_res))
  }else{
    mr_res <- tribble(
  ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval, ~mr,
 " ", " ", " ", " "," "," ", " ", " ", " ","rev",
)
all_mr_res <- mr(harmonise, metho=c("mr_ivw", "mr_wald_ratio"))
return(c(pval=mr_res, all=all_mr_res))
  }}
```

Generating nested case control analysis

``` r
nested_case_control <- function(dat, nestedcase, nestedcontrol){
    phen <- dat$phen
 nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$r0) %>% as_tibble() %>% mutate(x="r0", y="d") #case control estimate 
  #fast_nested<- fast_assoc(y=nested$d, x=nested$r0) %>% as_tibble() %>% mutate(x="r0", y="d") #case control estimate 
  return(c(casecontrol=log_nested))} #fast_cc=fast_nested
```

Generating two-sample forward MR analysis

``` r
fwd_mr_prot <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_join[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$c0
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_join[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$d #change this to d and do logistic regression to check this as it's binary
  g <- geno_mat
    gwas_dat <- data.frame()
 gwasy <- gwas(y,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  c_snp <- harmonise[101,] #select simulated causal SNP
  mr_res <- mr(c_snp, metho=c("mr_ivw", "mr_wald_ratio"))
  return(mr_res)}
```

Check two-sample reverse MR in all individuals of population

``` r
all_assoc <- function(dat){
y <- dat$phen$r0
  g <- dat$geno_prs
  gwasy <- gwas(y,g, logistic = FALSE)
  x <- dat$phen$d #change this to d and do logistic regression to check this as it's binary, l for fast assoc
  g <- dat$geno_prs
gwasx <- gwas(x,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  prs_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)#select SNPs reaching threshold
   mr_res <- mr(prs_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
  all_mr_res <- mr(harmonise, metho=c("mr_ivw", "mr_wald_ratio"))
  return(c(real_pval=mr_res, real=all_mr_res))
  }
```

Reverse MR for P1, P2 and P3 (causal, non-causal, consequence of disease
proteins)

``` r
## p1 revmr

rev_mr_prot_p1 <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_gc[ind, ] #change this to gc
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$c0 #change this for c0
  g <- geno_mat
  gwasy <- gwas(y,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_prs[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$d #change this to d and do logistic regression to check this as it's binary, l for fast assoc
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  prs_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)#select SNPs reaching threshold
  print(prs_snp)
  if (dim(prs_snp)[1]>0){
    mr_res <- mr(prs_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
    return(mr_res)
  }else{
    mr_res <- tribble(
      ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval, 
      0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
    return(mr_res)
  }}

## p2 revmr

rev_mr_prot_p2 <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_prs[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  phenosubset<- data.frame(phenosubset)
  y <- phenosubset$r0
  g <- geno_mat
  gwasy <- gwas(y,g, logistic = FALSE)
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_prs[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$d #change this to d and do logistic regression to check this as it's binary, l for fast assoc
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  prs_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)#select SNPs reaching threshold
  print(prs_snp)
  if (dim(prs_snp)[1]>0){
    mr_res <- mr(prs_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
    return(mr_res)
  }else{
    mr_res <- tribble(
      ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval, 
      0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
    return(mr_res)
  }}

## p3 revmr

rev_mr_prot_p3 <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  #d[d=="1"] <- 2 ##check this converting all disease phenotype to 2 as proteins consequence of disease only in diseased individuals?
  genosubset <- dat$geno_prs[ind, ] #CHECK THIS
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$x1 
  g <- geno_mat
  gwasy <- gwas(y,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_prs[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$d #change this to d and do logistic regression to check this as it's binary, l for fast assoc
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  prs_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)#select SNPs reaching threshold
  print(prs_snp)
  if (dim(prs_snp)[1]>0){
    mr_res <- mr(prs_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
    return(mr_res)
  }else{
    mr_res <- tribble(
      ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
    return(mr_res)
  }}
```

Forward MR for P1, P2 and P3 (causal, non-causal, consequence of disease
proteins)

``` r
##p1 fwd

fwd_mr_prot_p1 <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_gc[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$c0
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_join[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$d #change this to d and do logistic regression to check this as it's binary
  g <- geno_mat
  gwas_dat <- data.frame()
  gwasy <- gwas(y,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  p_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)
  if (dim(p_snp)[1]>0){
    mr_res <- mr(p_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
    return(mr_res)
  }else{
    mr_res <- tribble(
      ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval,
     0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
    return(mr_res)
  }}

##p2 fwd

fwd_mr_prot_p2 <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_prs[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$r0
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_join[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$d #change this to d and do logistic regression to check this as it's binary
  g <- geno_mat
  gwas_dat <- data.frame()
  gwasy <- gwas(y,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  p_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)
  if (dim(p_snp)[1]>0){
    mr_res <- mr(p_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
    return(mr_res)
  }else{
    mr_res <- tribble(
      ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
    return(mr_res)
  }}

##p3 fwd

fwd_mr_prot_p3 <- function(dat, n_protein_gwas, ncontrol, ncase){
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_prs[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$x1
  g <- geno_mat
  gwasx <- gwas(x,g, logistic = FALSE)
  
  no_prot_phen <- dat$phen[!ind,] #remove protein individuals so no overlap
  no_prot_gen <- data.matrix(dat$geno_join[!ind,])
  obs <- bind_rows(no_prot_phen[no_prot_phen$d == 0,][1:ncontrol,], no_prot_phen[no_prot_phen$d == 1,][1:ncase,])
  ind_obs <- do.call(paste0, no_prot_phen) %in% do.call(paste0, obs)
  phenosubset <- no_prot_phen[ind_obs, ]
  genosubset <- no_prot_gen[ind_obs, ]
  geno_mat <- data.matrix(genosubset)
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$d #change this to d and do logistic regression to check this as it's binary
  g <- geno_mat
  gwas_dat <- data.frame()
  gwasy <- gwas(y,g, logistic = TRUE)
  harmonise <- merge_exp_out(gwasx, gwasy, "X", "Y")
  p_snp <- subset(harmonise, harmonise$pval.exposure<5e-8)
  if (dim(p_snp)[1]>0){
    mr_res <- mr(p_snp, metho=c("mr_ivw", "mr_wald_ratio")) 
    return(mr_res)
  }else{
    mr_res <- tribble(
      ~id.exposure, ~id.outcome, ~outcome, ~exposure, ~method, ~nsnp,~b, ~se, ~pval,
      0, 0, 0, 0, 0, 0, 0, 0, 0,
    )
    return(mr_res)
  }}
```

CC for P1, P2 and P3 (causal, non-causal, consequence of disease
proteins)

``` r
##cc p1

nested_case_control_p1 <- function(dat, nestedcase, nestedcontrol){
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$c0)  #case control estimate 
  return(log_nested)} #fast_cc=fast_nested

##cc p2

nested_case_control_p2 <- function(dat, nestedcase, nestedcontrol){
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$r0)  #case control estimate 
  return(log_nested)} #fast_cc=fast_nested

##cc p3

nested_case_control_p3 <- function(dat, nestedcase, nestedcontrol){
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$x1) #case control estimate 
  return(log_nested)} #fast_cc=fast_nested
```

Input betas from training data into testing data to output model in
disease prediction

``` r
score_model <- function(betas, testing_dat, nid){
score <- betas$b[1]*testing_dat$phen$c0 + betas$b[2]*testing_dat$phen$r0 + betas$b[3]*testing_dat$phen$x1
lm_model <- lm(testing_dat$phen$d ~ score)
model <- coef(lm_model)[2]
return(model)
}
```

Input betas from testing data into training data to output model in
disease prediction - using data from linear model not MR

``` r
score_model_cc <- function(betas, testing_dat, nid){
betas$bhat <- as.numeric(betas$bhat)
score <- betas$bhat[1]*testing_dat$phen$c0 + betas$bhat[2]*testing_dat$phen$r0 + betas$bhat[3]*testing_dat$phen$x1
lm_model <- lm(testing_dat$phen$d ~ score)
model <- coef(lm_model)[2]
return(model)
}
```
