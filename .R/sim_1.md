Simulation 1
================
Lily Andrews
2023-11-28

Install packages Install packages

``` r
library(tibble)
library(simulateGP)
library(dplyr)
library(TwoSampleMR)
library(parallel)
library(data.table)
library(ggplot2)
```

Simulation 1

Data generating model using results from UK Biobank pQTL dataset
(<https://www.nature.com/articles/s41586-023-06592-6>) and glioma GWAS
(<https://www.nature.com/articles/ng.3823>):

``` r
rmarkdown::render("sims_functions.Rmd")
```

    ## 
    ## 
    ## processing file: sims_functions.Rmd

    ##   |                                                           |                                                   |   0%  |                                                           |..                                                 |   3%                     |                                                           |...                                                |   6% [unnamed-chunk-5]   |                                                           |.....                                              |  10%                     |                                                           |.......                                            |  13% [unnamed-chunk-6]   |                                                           |........                                           |  16%                     |                                                           |..........                                         |  19% [unnamed-chunk-7]   |                                                           |............                                       |  23%                     |                                                           |.............                                      |  26% [unnamed-chunk-8]   |                                                           |...............                                    |  29%                     |                                                           |................                                   |  32% [unnamed-chunk-9]   |                                                           |..................                                 |  35%                     |                                                           |....................                               |  39% [unnamed-chunk-10]  |                                                           |.....................                              |  42%                     |                                                           |.......................                            |  45% [unnamed-chunk-11]  |                                                           |.........................                          |  48%                     |                                                           |..........................                         |  52% [unnamed-chunk-12]  |                                                           |............................                       |  55%                     |                                                           |..............................                     |  58% [unnamed-chunk-13]  |                                                           |...............................                    |  61%                     |                                                           |.................................                  |  65% [unnamed-chunk-14]  |                                                           |...................................                |  68%                     |                                                           |....................................               |  71% [unnamed-chunk-15]  |                                                           |......................................             |  74%                     |                                                           |.......................................            |  77% [unnamed-chunk-16]  |                                                           |.........................................          |  81%                     |                                                           |...........................................        |  84% [unnamed-chunk-17]  |                                                           |............................................       |  87%                     |                                                           |..............................................     |  90% [unnamed-chunk-18]  |                                                           |................................................   |  94%                     |                                                           |.................................................  |  97% [unnamed-chunk-19]  |                                                           |...................................................| 100%                   

    ## output file: sims_functions.knit.md

    ## "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc" +RTS -K512m -RTS sims_functions.knit.md --to gfm+tex_math_dollars-yaml_metadata_block --from markdown+autolink_bare_uris+tex_math_single_backslash --output sims_functions.md --template "C:\Users\landr\AppData\Local\R\win-library\4.3\rmarkdown\rmarkdown\templates\github_document\resources\default.md" 
    ## "C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/pandoc" +RTS -K512m -RTS sims_functions.md --to html4 --from gfm+tex_math_dollars --output sims_functions.html --embed-resources --standalone --highlight-style pygments --template "C:\Users\landr\AppData\Local\R\win-library\4.3\rmarkdown\rmarkdown\templates\github_document\resources\preview.html" --variable "github-markdown-css:C:\Users\landr\AppData\Local\R\win-library\4.3\rmarkdown\rmarkdown\templates\github_document\resources\github.css" --metadata pagetitle=PREVIEW --mathjax

    ## 
    ## Preview created: C:\Users\landr\AppData\Local\Temp\RtmpgtSzqw\preview-422c20205f48.html

    ## 
    ## Output created: sims_functions.md

``` r
param  <- expand.grid(
  n_prot=34557,
  af = c(0.2,0.4),
  b=0.0879946, 
  rsq_lr0=0.0052923,
  nestedcasecontrol=2000,
  x1=0.1,
  sims=1:2
)

sims <- lapply(1:nrow(param), function(i){
  
print(param[i,])
print("before training dat")
training_dat <- dgmodel(nid=60000, 
b_gcc0=0.13,
b_u2c0=0,
gr_maf=param$af[i],
nsnp=99,
h2_known=0.2,
h2_unknown=0.19,
b_c0l=0,
b_u2l=0,
rsq_gr0=0,
rsq_lr0=param$rsq_lr0[i],
b_u1r0=0,
b_u1d=0,
b_u2d=0,
d_prev=3/10,
b_c0c1=1,
b_u2c1=0,
b_dc1=1,
b_r0r1=1,
b_u1r1=0,
b_dr1=0,
b_dx1=param$x1[i]) 
print("after training dat")

p1_rev <- rev_mr_prot_p1(dat=training_dat, n_protein_gwas=param$n_prot[i], ncontrol=18190, ncase=12496)
p2_rev <- rev_mr_prot_p2(dat=training_dat, n_protein_gwas=param$n_prot[i], ncontrol=18190, ncase=12496)
p3_rev <- rev_mr_prot_p3(dat=training_dat, n_protein_gwas=param$n_prot[i], ncontrol=18190, ncase=12496)
betas_rev <- rbind(p1_rev, p2_rev, p3_rev)
betas_rev<-as.data.frame(betas_rev)
for (row in 1:nrow(betas_rev)){
  if(betas_rev$pval[row]<0.05){
   
  } else{
    betas_rev[row,] <- list(0, 0, 0, 0, 0, 0, 0, 0, 0)
}}

p1_fwd <- fwd_mr_prot_p1(dat=training_dat, n_protein_gwas=param$n_prot[i], ncontrol=18190, ncase=12496)
p2_fwd <- fwd_mr_prot_p2(dat=training_dat, n_protein_gwas=param$n_prot[i], ncontrol=18190, ncase=12496)
p3_fwd <- fwd_mr_prot_p3(dat=training_dat, n_protein_gwas=param$n_prot[i], ncontrol=18190, ncase=12496)
betas_fwd <- rbind(p1_fwd, p2_fwd, p3_fwd)
betas_fwd<-as.data.frame(betas_fwd)
for (row in 1:nrow(betas_fwd)){
  if(betas_fwd$pval[row]<0.05){
   
  } else{
    betas_fwd[row,] <- list(0, 0, 0, 0, 0, 0, 0, 0, 0)
}}

p1_cc <- nested_case_control_p1(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
p2_cc <- nested_case_control_p2(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
p3_cc <- nested_case_control_p3(dat=training_dat, nestedcase=param$nestedcasecontrol[i], nestedcontrol=param$nestedcasecontrol[i])
cc <- rbind(p1_cc, p2_cc, p3_cc)
betas_cc<- as.data.frame(cc)
for (row in 1:nrow(betas_cc)){
  if(betas_cc$pval[row]<0.05){
  }
  else{
    betas_cc[row,] <- list(0, 0, 0, 0, 0, 0)
}}

print("made betas")
testing_dat <- dgmodel(nid=60000, 
b_gcc0=0.13,
b_u2c0=0,
gr_maf=param$af[i],
nsnp=99,
h2_known=0.2,
h2_unknown=0.19,
b_c0l=0,
b_u2l=0,
rsq_gr0=0,
rsq_lr0=param$rsq_lr0[i],
b_u1r0=0,
b_u1d=0,
b_u2d=0,
d_prev=3/10,
b_c0c1=1,
b_u2c1=0,
b_dc1=1,
b_r0r1=1,
b_u1r1=0,
b_dr1=0,
b_dx1=param$x1[i])

print("made tresting dat")
protein_model_rev <- score_model(betas=betas_rev, testing_dat, nid)
protein_model_fwd <- score_model(betas=betas_fwd, testing_dat, nid)
protein_model_cc <- score_model_cc(betas=betas_cc, testing_dat, nid)
print("made protein models")
pscore_df <- data.frame(rev = protein_model_rev, fwd = protein_model_fwd, cc = protein_model_cc)
print(pscore_df)
return(pscore_df) 
}) 
```

    ##   n_prot  af         b   rsq_lr0 nestedcasecontrol  x1 sims
    ## 1  34557 0.2 0.0879946 0.0052923              2000 0.1    1
    ## [1] "before training dat"
    ## [1] "after training dat"
    ## # A tibble: 0 × 18
    ## # ℹ 18 variables: SNP <int>, exposure <chr>, id.exposure <chr>, outcome <chr>,
    ## #   id.outcome <chr>, beta.exposure <dbl>, beta.outcome <dbl>,
    ## #   se.exposure <dbl>, se.outcome <dbl>, pval.exposure <dbl>,
    ## #   pval.outcome <dbl>, samplesize.exposure <dbl>, samplesize.outcome <dbl>,
    ## #   units.exposure <chr>, units.outcome <chr>, rsq.exposure <dbl>,
    ## #   rsq.outcome <dbl>, mr_keep <lgl>
    ## # A tibble: 14 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1     2 X        X           Y       Y                 -0.134    -0.0208  
    ##  2     6 X        X           Y       Y                 -0.370    -0.0244  
    ##  3    13 X        X           Y       Y                  0.159     0.0267  
    ##  4    15 X        X           Y       Y                  0.138     0.0193  
    ##  5    32 X        X           Y       Y                 -0.171    -0.0178  
    ##  6    40 X        X           Y       Y                 -0.136     0.00182 
    ##  7    42 X        X           Y       Y                 -0.143    -0.0117  
    ##  8    55 X        X           Y       Y                  0.159     0.0195  
    ##  9    57 X        X           Y       Y                 -0.166     0.000481
    ## 10    58 X        X           Y       Y                  0.194     0.0237  
    ## 11    69 X        X           Y       Y                  0.351     0.0167  
    ## 12    73 X        X           Y       Y                  0.142     0.00155 
    ## 13    81 X        X           Y       Y                 -0.158     0.00286 
    ## 14    86 X        X           Y       Y                  0.231     0.0112  
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'

    ## # A tibble: 14 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1     2 X        X           Y       Y                 -0.134    -0.00946 
    ##  2     6 X        X           Y       Y                 -0.370    -0.00259 
    ##  3    13 X        X           Y       Y                  0.159     0.0100  
    ##  4    15 X        X           Y       Y                  0.138     0.0172  
    ##  5    32 X        X           Y       Y                 -0.171     0.00569 
    ##  6    40 X        X           Y       Y                 -0.136     0.0121  
    ##  7    42 X        X           Y       Y                 -0.143    -0.00630 
    ##  8    55 X        X           Y       Y                  0.159     0.0166  
    ##  9    57 X        X           Y       Y                 -0.166    -0.000587
    ## 10    58 X        X           Y       Y                  0.194     0.00920 
    ## 11    69 X        X           Y       Y                  0.351     0.0174  
    ## 12    73 X        X           Y       Y                  0.142    -0.00321 
    ## 13    81 X        X           Y       Y                 -0.158     0.00714 
    ## 14    86 X        X           Y       Y                  0.231     0.00427 
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'

    ## Analysing 'X' on 'Y'

    ## [1] "made betas"
    ## [1] "made tresting dat"
    ## [1] "made protein models"
    ##             rev          fwd        cc
    ## score 0.3293597 0.0002927652 0.1824203
    ##   n_prot  af         b   rsq_lr0 nestedcasecontrol  x1 sims
    ## 2  34557 0.4 0.0879946 0.0052923              2000 0.1    1
    ## [1] "before training dat"
    ## [1] "after training dat"
    ## # A tibble: 0 × 18
    ## # ℹ 18 variables: SNP <int>, exposure <chr>, id.exposure <chr>, outcome <chr>,
    ## #   id.outcome <chr>, beta.exposure <dbl>, beta.outcome <dbl>,
    ## #   se.exposure <dbl>, se.outcome <dbl>, pval.exposure <dbl>,
    ## #   pval.outcome <dbl>, samplesize.exposure <dbl>, samplesize.outcome <dbl>,
    ## #   units.exposure <chr>, units.outcome <chr>, rsq.exposure <dbl>,
    ## #   rsq.outcome <dbl>, mr_keep <lgl>
    ## # A tibble: 17 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1     7 X        X           Y       Y                  0.124      0.0164 
    ##  2    11 X        X           Y       Y                  0.123      0.0134 
    ##  3    12 X        X           Y       Y                  0.168     -0.00462
    ##  4    13 X        X           Y       Y                  0.109     -0.00160
    ##  5    17 X        X           Y       Y                 -0.135     -0.0119 
    ##  6    27 X        X           Y       Y                  0.197      0.0103 
    ##  7    36 X        X           Y       Y                 -0.137     -0.0231 
    ##  8    42 X        X           Y       Y                 -0.140     -0.00953
    ##  9    47 X        X           Y       Y                  0.109      0.00969
    ## 10    51 X        X           Y       Y                  0.170      0.0161 
    ## 11    61 X        X           Y       Y                 -0.107     -0.00574
    ## 12    68 X        X           Y       Y                  0.138      0.0119 
    ## 13    70 X        X           Y       Y                  0.197     -0.00303
    ## 14    75 X        X           Y       Y                 -0.247     -0.0142 
    ## 15    77 X        X           Y       Y                 -0.134     -0.0107 
    ## 16    79 X        X           Y       Y                 -0.139     -0.0152 
    ## 17    97 X        X           Y       Y                 -0.143     -0.0102 
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'

    ## # A tibble: 17 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1     7 X        X           Y       Y                  0.124     0.0101  
    ##  2    11 X        X           Y       Y                  0.123     0.00527 
    ##  3    12 X        X           Y       Y                  0.168    -0.00189 
    ##  4    13 X        X           Y       Y                  0.109    -0.00335 
    ##  5    17 X        X           Y       Y                 -0.135    -0.00539 
    ##  6    27 X        X           Y       Y                  0.197    -0.00522 
    ##  7    36 X        X           Y       Y                 -0.137     0.00240 
    ##  8    42 X        X           Y       Y                 -0.140    -0.00763 
    ##  9    47 X        X           Y       Y                  0.109     0.0120  
    ## 10    51 X        X           Y       Y                  0.170     0.00291 
    ## 11    61 X        X           Y       Y                 -0.107     0.000877
    ## 12    68 X        X           Y       Y                  0.138     0.00641 
    ## 13    70 X        X           Y       Y                  0.197    -0.000978
    ## 14    75 X        X           Y       Y                 -0.247    -0.0273  
    ## 15    77 X        X           Y       Y                 -0.134     0.00362 
    ## 16    79 X        X           Y       Y                 -0.139     0.00332 
    ## 17    97 X        X           Y       Y                 -0.143     0.00168 
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'
    ## Analysing 'X' on 'Y'

    ## [1] "made betas"
    ## [1] "made tresting dat"
    ## [1] "made protein models"
    ##             rev fwd        cc
    ## score 0.3267515  NA 0.1894185
    ##   n_prot  af         b   rsq_lr0 nestedcasecontrol  x1 sims
    ## 3  34557 0.2 0.0879946 0.0052923              2000 0.1    2
    ## [1] "before training dat"
    ## [1] "after training dat"
    ## # A tibble: 0 × 18
    ## # ℹ 18 variables: SNP <int>, exposure <chr>, id.exposure <chr>, outcome <chr>,
    ## #   id.outcome <chr>, beta.exposure <dbl>, beta.outcome <dbl>,
    ## #   se.exposure <dbl>, se.outcome <dbl>, pval.exposure <dbl>,
    ## #   pval.outcome <dbl>, samplesize.exposure <dbl>, samplesize.outcome <dbl>,
    ## #   units.exposure <chr>, units.outcome <chr>, rsq.exposure <dbl>,
    ## #   rsq.outcome <dbl>, mr_keep <lgl>
    ## # A tibble: 17 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1    10 X        X           Y       Y                  0.182     0.00153 
    ##  2    11 X        X           Y       Y                  0.144     0.00754 
    ##  3    16 X        X           Y       Y                  0.196     0.00684 
    ##  4    17 X        X           Y       Y                  0.168     0.00654 
    ##  5    24 X        X           Y       Y                  0.159     0.0184  
    ##  6    27 X        X           Y       Y                 -0.144    -0.00980 
    ##  7    41 X        X           Y       Y                  0.182     0.0150  
    ##  8    42 X        X           Y       Y                  0.170     0.0121  
    ##  9    51 X        X           Y       Y                  0.151     0.0213  
    ## 10    58 X        X           Y       Y                 -0.160     0.000337
    ## 11    60 X        X           Y       Y                  0.169     0.0258  
    ## 12    62 X        X           Y       Y                  0.220     0.00564 
    ## 13    65 X        X           Y       Y                  0.244     0.0103  
    ## 14    78 X        X           Y       Y                 -0.136     0.0274  
    ## 15    80 X        X           Y       Y                  0.136    -0.0112  
    ## 16    81 X        X           Y       Y                 -0.173    -0.00480 
    ## 17    97 X        X           Y       Y                 -0.138    -0.0130  
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'

    ## # A tibble: 17 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1    10 X        X           Y       Y                  0.182      0.00855
    ##  2    11 X        X           Y       Y                  0.144      0.0154 
    ##  3    16 X        X           Y       Y                  0.196      0.0139 
    ##  4    17 X        X           Y       Y                  0.168     -0.00681
    ##  5    24 X        X           Y       Y                  0.159      0.00748
    ##  6    27 X        X           Y       Y                 -0.144     -0.00766
    ##  7    41 X        X           Y       Y                  0.182      0.0121 
    ##  8    42 X        X           Y       Y                  0.170      0.00554
    ##  9    51 X        X           Y       Y                  0.151      0.00389
    ## 10    58 X        X           Y       Y                 -0.160     -0.00111
    ## 11    60 X        X           Y       Y                  0.169      0.00972
    ## 12    62 X        X           Y       Y                  0.220     -0.00690
    ## 13    65 X        X           Y       Y                  0.244      0.0305 
    ## 14    78 X        X           Y       Y                 -0.136      0.00794
    ## 15    80 X        X           Y       Y                  0.136      0.00493
    ## 16    81 X        X           Y       Y                 -0.173     -0.0122 
    ## 17    97 X        X           Y       Y                 -0.138     -0.0163 
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'
    ## Analysing 'X' on 'Y'

    ## [1] "made betas"
    ## [1] "made tresting dat"
    ## [1] "made protein models"
    ##             rev         fwd        cc
    ## score 0.4700381 0.003933513 0.1702689
    ##   n_prot  af         b   rsq_lr0 nestedcasecontrol  x1 sims
    ## 4  34557 0.4 0.0879946 0.0052923              2000 0.1    2
    ## [1] "before training dat"
    ## [1] "after training dat"
    ## # A tibble: 1 × 18
    ##     SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##   <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ## 1     1 X        X           Y       Y                 -0.199        0.136
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'

    ## # A tibble: 22 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1     1 X        X           Y       Y                 -0.199     -0.0153 
    ##  2     7 X        X           Y       Y                  0.117      0.00205
    ##  3     8 X        X           Y       Y                  0.107      0.0171 
    ##  4    14 X        X           Y       Y                  0.152      0.0187 
    ##  5    20 X        X           Y       Y                 -0.149     -0.00913
    ##  6    22 X        X           Y       Y                  0.107      0.0139 
    ##  7    36 X        X           Y       Y                 -0.137      0.00192
    ##  8    38 X        X           Y       Y                 -0.131     -0.0108 
    ##  9    44 X        X           Y       Y                  0.145      0.00956
    ## 10    53 X        X           Y       Y                 -0.130     -0.00390
    ## # ℹ 12 more rows
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'

    ## # A tibble: 22 × 18
    ##      SNP exposure id.exposure outcome id.outcome beta.exposure beta.outcome
    ##    <int> <chr>    <chr>       <chr>   <chr>              <dbl>        <dbl>
    ##  1     1 X        X           Y       Y                 -0.199      0.00358
    ##  2     7 X        X           Y       Y                  0.117      0.00425
    ##  3     8 X        X           Y       Y                  0.107     -0.00207
    ##  4    14 X        X           Y       Y                  0.152      0.00867
    ##  5    20 X        X           Y       Y                 -0.149      0.00311
    ##  6    22 X        X           Y       Y                  0.107      0.00824
    ##  7    36 X        X           Y       Y                 -0.137     -0.00497
    ##  8    38 X        X           Y       Y                 -0.131     -0.00298
    ##  9    44 X        X           Y       Y                  0.145      0.00600
    ## 10    53 X        X           Y       Y                 -0.130      0.0238 
    ## # ℹ 12 more rows
    ## # ℹ 11 more variables: se.exposure <dbl>, se.outcome <dbl>,
    ## #   pval.exposure <dbl>, pval.outcome <dbl>, samplesize.exposure <dbl>,
    ## #   samplesize.outcome <dbl>, units.exposure <chr>, units.outcome <chr>,
    ## #   rsq.exposure <dbl>, rsq.outcome <dbl>, mr_keep <lgl>

    ## Analysing 'X' on 'Y'
    ## Analysing 'X' on 'Y'

    ## [1] "made betas"
    ## [1] "made tresting dat"
    ## [1] "made protein models"
    ##                 rev          fwd        cc
    ## score -0.0006449387 -0.001298712 0.2864647

``` r
pscore <- data.table()
for (i in 1:nrow(param)){
    row <- sims[[i]]
    pscore <- rbind(pscore, row)
}
```
