library(tibble)
library(simulateGP)
library(dplyr)
library(TwoSampleMR)
library(parallel)
library(data.table)
library(ggplot2)
library(pROC)
library(stringr)
library(readxl)

glioma_revmr_rsq <- read_excel("glioma_revmr_rsq.xlsx")
#glioma_revmr_rsq <- glioma_revmr_rsq[order(glioma_revmr_rsq$b, decreasing = TRUE),] #to get proteins with largest effect first
glioma_revmr_rsq$rsq_gr0 <-0.00574511 #-0.0009103566 HMGSC1#0.00574511 sum of bcan rsq 
glioma_revmr_rsq$af <-0.2614#0.2614859 HMGSC1#0.2614 mean af bcan
glioma_prs <- read.delim("ukbb_prs_val.txt")
fval <- (glioma_prs$beta.exposure/glioma_prs$se.exposure)^2
glioma_prs$rsq <- fval/(30686-2+fval)
#To generate rsq from rev MR results
#fval <- (beta/se)^2
#rsq <- fval/(n-2+fval)
#param  <- expand.grid(
#  af = seq(0.2, 0.4, by=0.1),
#  nestedcasecontrol=c(113, 487,593), # what sample size of cc is needed to match revmr results
#  model = 1:nrow(glioma_revmr_rsq),
#  x1=seq(0, 1, by=0.2), # see effect of consequence of disease - increase and see if prediction/auc decreases
#  sims=1:3,
#  uc= seq(0, 0.4, by=0.1) #see if uc improves prediction of disease
#)
param  <- expand.grid(
  nestedcasecontrol=seq(0, 300, by=50), # what sample size of cc is needed to match revmr results
  model = 1,
  sims=1:15,
  effect = seq(0.05, 0.4, by=0.05),
  uc= seq(0, 0.4, by=0.05) #see if uc improves prediction of disease
)

source("dgmodel.R")
source("estimation.R")

sims <- lapply(1:nrow(param), function(i)
{
  training_dat <- dgmodel(
    nid=150000,
    b_gcc0= param$effect[i],#9.713741e-13,#beta -1.63622e-06#now trying PHLDB1 from fwd MR #0.01004,from ukbb EGFR beta rs75061358 not divided by 99
    b_u2c0=param$uc[i],
    gr_maf=glioma_revmr_rsq$af[param$model[i]],
    nsnp=99,
    h2_known=0.2,
    h2_unknown=0.19,
    b_c0l= param$effect[i],#7.6775e-4,#now trying PHLDB1 from fwd MR #0.341748,gwas SNP beta-ukbb SNP beta not divided by 99
    b_u2l=param$uc[i],
    rsq_gr0=param$effect[i], ##check this
    rsq_lr0=param$effect[i],#glioma_revmr_rsq$rsq[param$model[i]],
    b_u1r0=param$uc[i],
    b_u1l=param$uc[i],
    d_prev=4/10,
    b_c0c1=1,
    b_u2c1=param$uc[i],
    b_dc1=0,
    b_r0r1=1,
    b_u1r1=param$uc[i],
    b_dr1=0, #no change for r1 on disease
    b_dx1=param$effect[i],
    b_x0x1=0,#changed this to 0 check with gib
    gprs_rsq=rep(0.1, 99),#glioma_prs$rsq, 
    gprs_maf=glioma_prs$eaf.exposure
  )
  
  p1_rev <- rev_mr_prot_p1(dat=training_dat, n_protein_gwas=glioma_revmr_rsq$n[param$model[i]], ncontrol=18190, ncase=12496)
  p2_rev <- rev_mr_prot_p2(dat=training_dat, n_protein_gwas=glioma_revmr_rsq$n[param$model[i]], ncontrol=18190, ncase=12496)
  p3_rev <- rev_mr_prot_p3(dat=training_dat, n_protein_gwas=glioma_revmr_rsq$n[param$model[i]], ncontrol=18190, ncase=12496)
  betas_rev <- rbind(p1_rev, p2_rev, p3_rev)
  betas_rev<-as.data.frame(betas_rev)
  rev_prot <- data.frame()
  for (row in 1:nrow(betas_rev))
  {
    if(betas_rev$pval[row]<0.05)
    {
    }
    else
    {
      betas_rev[row,] <- list(0, 0, 0, 0, 0, 0, 0, 0, 0)
    }}
  for (row in 1:nrow(betas_rev))
  {
    if(betas_rev$method[row]!=0)
    {
      rev_prot <- rbind(rev_prot, paste0("prot", row))
      colnames(rev_prot) <- "prot"
    }
    else
    {
    }}
  p1_fwd <- fwd_mr_prot_p1(dat=training_dat, n_protein_gwas=glioma_revmr_rsq$n[param$model[i]], ncontrol=18190, ncase=12496)
  p2_fwd <- fwd_mr_prot_p2(dat=training_dat, n_protein_gwas=glioma_revmr_rsq$n[param$model[i]], ncontrol=18190, ncase=12496)
  p3_fwd <- fwd_mr_prot_p3(dat=training_dat, n_protein_gwas=glioma_revmr_rsq$n[param$model[i]], ncontrol=18190, ncase=12496)
  betas_fwd <- rbind(p1_fwd, p2_fwd, p3_fwd)
  betas_fwd<-as.data.frame(betas_fwd)
  fwd_prot <- data.frame()
  for (row in 1:nrow(betas_fwd))
  {
    if(betas_fwd$pval[row]<0.05)
    {
    }
    else
    {
      betas_fwd[row,] <- list(0, 0, 0, 0, 0, 0, 0, 0, 0)
    }}
  for (row in 1:nrow(betas_fwd))
  {
    if(betas_fwd$method[row]!=0)
    {
      fwd_prot <- rbind(fwd_prot, paste0("prot", row))
      colnames(fwd_prot) <- "prot"
    }
    else
    {
    }}
  p1_ncc <- nested_case_control_p1(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
  p2_ncc <- nested_case_control_p2(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
  p3_ncc <- nested_case_control_p3(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
  ncc <- rbind(p1_ncc, p2_ncc, p3_ncc)
  betas_ncc<- as.data.frame(ncc)
  ncc_prot <- data.frame()
  for (row in 1:nrow(betas_ncc))
  {
    if(betas_ncc$pval[row]<0.05)
    {
    }
    else
    {
      betas_ncc[row,] <- list(0, 0, 0, 0, 0, 0)
    }}
  for (row in 1:nrow(betas_ncc))
  {
    if(betas_ncc$n[row]!=0)  
    {
      ncc_prot <- rbind(ncc_prot, paste0("prot", row))
      colnames(ncc_prot) <- "prot"
    }
    else
    {
    }}
  p1_cc <- case_control_p1(dat=training_dat, case=param$nestedcasecontrol[i], control=param$nestedcasecontrol[i])
  p2_cc <- case_control_p2(dat=training_dat, case=param$nestedcasecontrol[i], control=param$nestedcasecontrol[i])
  p3_cc <- case_control_p3(dat=training_dat, case=param$nestedcasecontrol[i], control=param$nestedcasecontrol[i])
  cc <- rbind(p1_cc, p2_cc, p3_cc)
  betas_cc<- as.data.frame(cc)
  cc_prot <- data.frame()
  for (row in 1:nrow(betas_cc))
  {
    if(betas_cc$pval[row]<0.05)
    {
    }
    else
    {
      betas_cc[row,] <- list(0, 0, 0, 0, 0, 0)
    }}
  for (row in 1:nrow(betas_cc))
  {
    if(betas_cc$n[row]!=0)  
    {
      cc_prot <- rbind(cc_prot, paste0("prot", row))
      colnames(cc_prot) <- "prot"
    }
    else
    {
    }}
  
  testing_dat <- dgmodel(
    nid=150000,
    b_gcc0= param$effect[i],#9.713741e-13,#beta -1.63622e-06#now trying PHLDB1 from fwd MR #0.01004,from ukbb EGFR beta rs75061358 not divided by 99
    b_u2c0=param$uc[i],
    gr_maf=glioma_revmr_rsq$af[param$model[i]],
    nsnp=99,
    h2_known=0.2,
    h2_unknown=0.19,
    b_c0l= param$effect[i],#7.6775e-4,#now trying PHLDB1 from fwd MR #0.341748,gwas SNP beta-ukbb SNP beta not divided by 99
    b_u2l=param$uc[i],
    rsq_gr0=param$effect[i], ##check this
    rsq_lr0=param$effect[i],#glioma_revmr_rsq$rsq[param$model[i]],
    b_u1r0=param$uc[i],
    b_u1l=param$uc[i],
    d_prev=4/10,
    b_c0c1=1,
    b_u2c1=param$uc[i],
    b_dc1=0,
    b_r0r1=1,
    b_u1r1=param$uc[i],
    b_dr1=0, #no change for r1 on disease
    b_dx1=param$effect[i],
    b_x0x1=0,#changed this to 0 check with gib
    gprs_rsq=rep(0.1, 99),#glioma_prs$rsq, 
    gprs_maf=glioma_prs$eaf.exposure
  )
  if(sum(unlist(betas_rev$b))!=0)
  {
    protein_model_rev <- score_model(betas=betas_rev, testing_dat, nid)
  } else
  {
    protein_model_rev <- list(0, 0, 0, 0)
  }
  if(sum(unlist(betas_fwd$b))!=0)
  {
    protein_model_fwd <- score_model(betas=betas_fwd, testing_dat, nid)
  } else
  {
    protein_model_fwd <- list(0, 0, 0, 0)
  }
  if(sum(unlist(betas_ncc$bhat))!=0)
  {
    protein_model_ncc <- score_model_ncc(betas=betas_ncc, testing_dat, nid)
  } else
  {
    protein_model_ncc <- list(0, 0, 0, 0)
  }
  if(sum(unlist(betas_cc$bhat))!=0)
  {
    protein_model_cc <- score_model_cc(betas=betas_cc, testing_dat, nid)
  } else
  {
    protein_model_cc <- list(0, 0, 0, 0)
  }
  
  pscore_df <- rbind(rev=protein_model_rev, fwd=protein_model_fwd, ncc=protein_model_ncc, cc=protein_model_cc)
  pscore_df <- as.data.frame(pscore_df)
  roc <- data.table(rev=0.5, fwd=0.5, ncc=0.5, cc=0.5)
  
  ##REVERSE
  betas <- betas_rev
  if(sum(unlist(betas$b))!=0)
  {
    testing_dat$score <- betas$b[1]*testing_dat$phen$c0 + betas$b[2]*testing_dat$phen$r0 + betas$b[3]*testing_dat$phen$x1
    roc_object <- roc(response=testing_dat$phen$d, predictor=testing_dat$score)
    rev_roc <- pROC::auc(roc_object)
    roc$rev <- rev_roc
  } else
  {
  }
  
  ##FORWARD
  betas <- betas_fwd
  if(sum(unlist(betas$b))!=0)
  {
    testing_dat$score <- betas$b[1]*testing_dat$phen$c0 + betas$b[2]*testing_dat$phen$r0 + betas$b[3]*testing_dat$phen$x1
    roc_object <- roc(response=testing_dat$phen$d, predictor=testing_dat$score)
    fwd_roc <- pROC::auc(roc_object)
    roc$fwd <- fwd_roc
  } else
  {
  }
  
  ## NESTED CASE CONTROLS
  betas <- betas_ncc
  print(betas)
  if(sum(unlist(betas$bhat))!=0)
  {
    betas$bhat <- as.numeric(betas$bhat)
    testing_dat$score <- betas$bhat[1]*testing_dat$phen$c0 + betas$bhat[2]*testing_dat$phen$r0 +   betas$bhat[3]*testing_dat$phen$x1
    roc_object <- roc(response=testing_dat$phen$d, predictor=testing_dat$score)
    ncc_roc <- pROC::auc(roc_object)
    roc$ncc <- ncc_roc} else
    {
    }
  ## CASE CONTROLS
  betas <- betas_cc
  print(betas)
  if(sum(unlist(betas$bhat))!=0)
  {
    betas$bhat <- as.numeric(betas$bhat)
    testing_dat$score <- betas$bhat[1]*testing_dat$phen$c0 + betas$bhat[2]*testing_dat$phen$r0 +   betas$bhat[3]*testing_dat$phen$x1
    roc_object <- roc(response=testing_dat$phen$d, predictor=testing_dat$score)
    cc_roc <- pROC::auc(roc_object)
    roc$cc <- cc_roc} else
    {
    }
  
  print(t(roc))
  out <- tibble(
    gr_af=rep(glioma_revmr_rsq$af[param$model[i]], 4),
    n_cc=rep(param$nestedcasecontrol[i],4),
    nid_prot=rep(glioma_revmr_rsq$n[param$model[i]],4),
    rsq_lr0=rep(glioma_revmr_rsq$rsq[param$model[i]],4),
    x1=rep(glioma_revmr_rsq$rsq[param$model[i]],4),
    sim=rep(param$sims[i],4),
    uc=rep(param$uc[i],4),
    method=row.names(pscore_df),
    score=c(as.numeric(pscore_df[1,1]),as.numeric(pscore_df[2,1]), as.numeric(pscore_df[3,1]), as.numeric(pscore_df[4,1])),
    se=c(as.numeric(pscore_df[1,2]),as.numeric(pscore_df[2,2]), as.numeric(pscore_df[3,2]), as.numeric(pscore_df[4,2])),
    pval=c(as.numeric(pscore_df[1,4]),as.numeric(pscore_df[2,4]), as.numeric(pscore_df[3,4]), as.numeric(pscore_df[4,4])),
    auc=c(as.numeric(t(roc)[1,]), as.numeric(t(roc)[2,]), as.numeric(t(roc)[3,]), as.numeric(t(roc)[4,])),
    n_proteins=c(nrow(rev_prot), nrow(fwd_prot), nrow(ncc_prot), nrow(cc_prot)),
    protein_type=c(str_flatten(rev_prot[1:nrow(rev_prot),], collapse = " "), str_flatten(fwd_prot[1:nrow(fwd_prot),], collapse = " "), str_flatten(ncc_prot[1:nrow(ncc_prot),], collapse = " "), str_flatten(cc_prot[1:nrow(cc_prot),], collapse = " "))
  )
  print(out)
  return(out)
}) %>% bind_rows()

save(sims, file = "sims_run1.RData")
