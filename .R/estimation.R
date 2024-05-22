##to check variables and check variables which should be 1

dgmodel_check <- function(dat)
{
  print(tibble(
    col = names(dat$phen),
    vars = apply(dat$phen, 2, var), 
    means = colMeans(dat$phen)
  ), n=25)
  print(cor(dat$phen$prs, dat$phen$l)^2)
}

## look at associations between different parameters

dgmodel_analysis <- function(dat, ncase, ncontrol, protein_gwas, prot, iteration)
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
  ) %>% mutate(study="all")
  
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
  )  %>% mutate(study="observational")
  
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
  ) %>% mutate(study="protein_gwas")
  return(bind_rows(res_obs%>%
                     mutate(study="observational"), res_protein%>%
                     mutate(study="protein_gwas")))
} 

##logistic association equation

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

## reverse MR of non-causal protein (protein cohort) and disease (case control cohort)

rev_mr_prot <- function(dat, n_protein_gwas, ncontrol, ncase)
{
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

## forward MR of causal protein (protein cohort) and disease (case control cohort)

fwd_mr_prot <- function(dat, n_protein_gwas, ncontrol, ncase)
{
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
  return(mr_res)
}

nested_case_control <- function(dat, nestedcase, nestedcontrol)
{
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$r0) %>% as_tibble() %>% mutate(x="r0", y="d") #case control estimate 
  return(c(casecontrol=log_nested))
}

##reverse MR between disease and r0 in all individuals

all_assoc <- function(dat)
{
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

## p1 revmr

rev_mr_prot_p1 <- function(dat, n_protein_gwas, ncontrol, ncase)
{
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

rev_mr_prot_p2 <- function(dat, n_protein_gwas, ncontrol, ncase)
{
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

rev_mr_prot_p3 <- function(dat, n_protein_gwas, ncontrol, ncase)
{
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  #d[d=="1"] <- 2 ##check this converting all disease phenotype to 2 as proteins consequence of disease only in diseased individuals?
  genosubset <- dat$geno_prs[ind, ] #CHECK THIS
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  y <- phenosubset$x0 
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

##p1 fwd

fwd_mr_prot_p1 <- function(dat, n_protein_gwas, ncontrol, ncase)
{
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

fwd_mr_prot_p2 <- function(dat, n_protein_gwas, ncontrol, ncase)
{
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

fwd_mr_prot_p3 <- function(dat, n_protein_gwas, ncontrol, ncase)
{
  phen <- dat$phen
  prot <- phen[1:n_protein_gwas, ] #first group of protein individuals prot <- phen[phen$d==0,][1:n_protein_gwas, ]
  ind <- do.call(paste0, dat$phen) %in% do.call(paste0, prot)
  genosubset <- dat$geno_prs[ind, ]
  geno_mat <- data.matrix(genosubset)
  phenosubset <- dat$phen[ind, ]
  stopifnot(nrow(genosubset) == nrow(phenosubset))
  x <- phenosubset$x0
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

##cc p1

nested_case_control_p1 <- function(dat, nestedcase, nestedcontrol)
{
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$c0)  #case control estimate 
  return(log_nested)
} 

##cc p2

nested_case_control_p2 <- function(dat, nestedcase, nestedcontrol)
{
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$r0)  #case control estimate 
  return(log_nested)
} 

##cc p3

nested_case_control_p3 <- function(dat, nestedcase, nestedcontrol)
{
  phen <- dat$phen
  nested <- bind_rows(phen[phen$d == 0,][1:nestedcase,], phen[phen$d == 1,][1:nestedcontrol,]) #487 example based on nested case control study
  log_nested<- logistic_assoc(y=nested$d, x=nested$x0) #case control estimate 
  return(log_nested)
} 

##create protein score and create model
score_model <- function(betas, testing_dat, nid)
{
  score <- betas$b[1]*testing_dat$phen$c0 + betas$b[2]*testing_dat$phen$r0 + betas$b[3]*testing_dat$phen$x0
  lm_model <- glm(testing_dat$phen$d ~ score, family = "binomial")
  model <- coef(summary(lm_model))[2,]
  return(model)
}

##create protein score and create model for case control analysis

score_model_cc <- function(betas, testing_dat, nid)
{
  betas$bhat <- as.numeric(betas$bhat)
  score <- betas$bhat[1]*testing_dat$phen$c0 + betas$bhat[2]*testing_dat$phen$r0 + betas$bhat[3]*testing_dat$phen$x0
  lm_model <- glm(testing_dat$phen$d ~ score, family = "binomial")
  model <- coef(summary(lm_model))[2,]
  return(model)
}