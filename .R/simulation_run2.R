library(tibble)
library(simulateGP)
library(dplyr)
library(TwoSampleMR)
library(parallel)
library(data.table)
library(ggplot2)
library(pROC)
library(stringr)

n_prot<-c(33436, 33187, 3301, 35362)
b<-c(0.0497226, 0.0447559, 0.0879946, -0.0354690)
rsq_lr0<-c(0.0016068, 0.0013336, 0.0052923, 0.0008731)
real_dat <- data.frame(n_prot, b, rsq_lr0)

param  <- expand.grid(
  af = 0.2,
  nestedcasecontrol=487,
  model = 1:nrow(real_dat),
  x1=seq(0, 0.4, by=0.1),
  sims=1,
  uc=c(0.1, 0.4)
)

source("dgmodel.R")
source("estimation.R")

sims_2 <- lapply(1:nrow(param), function(i)
{
  #print(param[i,])
  training_dat <- dgmodel(
    nid=60000, 
    b_gcc0=0.13,#need to sort this - should be 0.26 or 0.13
    b_u2c0=param$uc[i],
    gr_maf=param$af[i],
    nsnp=99,
    h2_known=0.2,
    h2_unknown=0.19,
    b_c0l=0.13,
    b_u2l=param$uc[i],
    rsq_gr0=0.1, ##check this
    rsq_lr0=real_dat$rsq_lr0[param$model[i]],
    b_u1r0=param$uc[i],
    b_u1l=param$uc[i],
    d_prev=3/10,
    b_c0c1=1,
    b_u2c1=param$uc[i],
    b_dc1=1,
    b_r0r1=1,
    b_u1r1=param$uc[i],
    b_dr1=0, #no change for r1 on disease
    b_dx1=param$x1[i],
    b_x0x1=0
  ) 
  
  p1_rev <- rev_mr_prot_p1(dat=training_dat, n_protein_gwas=real_dat$n_prot[param$model[i]], ncontrol=18190, ncase=12496)
  p2_rev <- rev_mr_prot_p2(dat=training_dat, n_protein_gwas=real_dat$n_prot[param$model[i]], ncontrol=18190, ncase=12496)
  p3_rev <- rev_mr_prot_p3(dat=training_dat, n_protein_gwas=real_dat$n_prot[param$model[i]], ncontrol=18190, ncase=12496)
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
  p1_fwd <- fwd_mr_prot_p1(dat=training_dat, n_protein_gwas=real_dat$n_prot[param$model[i]], ncontrol=18190, ncase=12496)
  p2_fwd <- fwd_mr_prot_p2(dat=training_dat, n_protein_gwas=real_dat$n_prot[param$model[i]], ncontrol=18190, ncase=12496)
  p3_fwd <- fwd_mr_prot_p3(dat=training_dat, n_protein_gwas=real_dat$n_prot[param$model[i]], ncontrol=18190, ncase=12496)
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
  p1_cc <- nested_case_control_p1(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
  p2_cc <- nested_case_control_p2(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
  p3_cc <- nested_case_control_p3(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
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
    nid=60000, 
    b_gcc0=0.13,#need to sort this - should be 0.26 or 0.13
    b_u2c0=param$uc[i],
    gr_maf=param$af[i],
    nsnp=99,
    h2_known=0.2,
    h2_unknown=0.19,
    b_c0l=0.13,
    b_u2l=param$uc[i],
    rsq_gr0=0.1, ##check this
    rsq_lr0=real_dat$rsq_lr0[param$model[i]],
    b_u1r0=param$uc[i],
    b_u1l=param$uc[i],
    d_prev=3/10,
    b_c0c1=1,
    b_u2c1=param$uc[i],
    b_dc1=1,
    b_r0r1=1,
    b_u1r1=param$uc[i],
    b_dr1=0, #no change for r1 on disease
    b_dx1=param$x1[i],
    b_x0x1=0
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
  if(sum(unlist(betas_cc$bhat))!=0)
  {
    protein_model_cc <- score_model_cc(betas=betas_cc, testing_dat, nid)
  } else
  {
    protein_model_cc <- list(0, 0, 0, 0)
  }
  
  pscore_df <- rbind(rev=protein_model_rev, fwd=protein_model_fwd, cc=protein_model_cc)
  pscore_df <- as.data.frame(pscore_df)
  roc <- data.table(rev=0.5, fwd=0.5, cc=0.5)
  
  ##REVERSE
  betas <- betas_rev
  if(sum(unlist(betas$b))!=0)
  {
    testing_dat$score <- betas$b[1]*testing_dat$phen$c0 + betas$b[2]*testing_dat$phen$r0 + betas$b[3]*testing_dat$phen$x0
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
    testing_dat$score <- betas$b[1]*testing_dat$phen$c0 + betas$b[2]*testing_dat$phen$r0 + betas$b[3]*testing_dat$phen$x0
    roc_object <- roc(response=testing_dat$phen$d, predictor=testing_dat$score)
    fwd_roc <- pROC::auc(roc_object)
    roc$fwd <- fwd_roc
  } else
  {
  }
  
  ## CASE CONTROLS
  betas <- betas_cc
  print(betas)
  if(sum(unlist(betas$bhat))!=0)
  {
    betas$bhat <- as.numeric(betas$bhat)
    testing_dat$score <- betas$bhat[1]*testing_dat$phen$c0 + betas$bhat[2]*testing_dat$phen$r0 +   betas$bhat[3]*testing_dat$phen$x0
    roc_object <- roc(response=testing_dat$phen$d, predictor=testing_dat$score)
    cc_roc <- pROC::auc(roc_object)
    roc$cc <- cc_roc} else
    {
    }
  print(t(roc))
  out <- tibble(
    af=rep(param$af[i], 3),
    n_cc=rep(param$nestedcasecontrol[i],3),
    nid_prot=rep(real_dat$n_prot[param$model[i]],3),
    rsq_lr0=rep(real_dat$rsq_lr0[param$model[i]],3),
    x1=rep(param$x1[i],3),
    sim=rep(param$sims[i],3),
    uc=rep(param$uc[i],3),
    method=row.names(pscore_df),
    score=c(as.numeric(pscore_df[1,1]),as.numeric(pscore_df[2,1]), as.numeric(pscore_df[3,1])),
    se=c(as.numeric(pscore_df[1,2]),as.numeric(pscore_df[2,2]), as.numeric(pscore_df[3,2])),
    pval=c(as.numeric(pscore_df[1,4]),as.numeric(pscore_df[2,4]), as.numeric(pscore_df[3,4])),
    auc=c(as.numeric(t(roc)[1,]), as.numeric(t(roc)[2,]), as.numeric(t(roc)[3,])),
    n_proteins=c(nrow(rev_prot), nrow(fwd_prot), nrow(cc_prot)),
    protein_type=c(str_flatten(rev_prot[1:nrow(rev_prot),], collapse = " "), str_flatten(fwd_prot[1:nrow(fwd_prot),], collapse = " "), str_flatten(cc_prot[1:nrow(cc_prot),], collapse = " ")) 
  )
  print(out)
  return(out) 
}) %>% bind_rows()


save(sims_2, file = "sims_run2.RData")