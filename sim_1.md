Simulation 1
================
Lily Andrews
2023-11-28

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

Simulation 1

Data generating model:

``` r
rmarkdown::render(".R/sims_functions.Rmd")
```

    ## 
    ## 
    ## processing file: sims_functions.Rmd

    ##   |                                                           |                                                   |   0%  |                                                           |....                                               |   8%                     |                                                           |........                                           |  15% [unnamed-chunk-5]   |                                                           |............                                       |  23%                     |                                                           |................                                   |  31% [unnamed-chunk-6]   |                                                           |....................                               |  38%                     |                                                           |........................                           |  46% [unnamed-chunk-7]   |                                                           |...........................                        |  54%                     |                                                           |...............................                    |  62% [unnamed-chunk-8]   |                                                           |...................................                |  69%                     |                                                           |.......................................            |  77% [unnamed-chunk-9]   |                                                           |...........................................        |  85%                     |                                                           |...............................................    |  92% [unnamed-chunk-10]  |                                                           |...................................................| 100%                   

    ## output file: sims_functions.knit.md

    ## /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/pandoc +RTS -K512m -RTS sims_functions.knit.md --to gfm+tex_math_dollars-yaml_metadata_block --from markdown+autolink_bare_uris+tex_math_single_backslash --output sims_functions.md --template /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown/rmarkdown/templates/github_document/resources/default.md 
    ## /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/pandoc +RTS -K512m -RTS sims_functions.md --to html4 --from gfm+tex_math_dollars --output sims_functions.html --embed-resources --standalone --highlight-style pygments --template /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown/rmarkdown/templates/github_document/resources/preview.html --variable 'github-markdown-css:/Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown/rmarkdown/templates/github_document/resources/github.css' --metadata pagetitle=PREVIEW --mathjax

    ## 
    ## Preview created: /var/folders/ml/4j6q26qn7r5b89gy3dssz_d40000gn/T//RtmpWpRxxF/preview-afe277600542.html

    ## 
    ## Output created: sims_functions.md

``` r
dat <- dgmodel(
  nid=100000, #number of individuals
  nsnp=99, #number of snps
  gc_maf=runif(99, 0.05, 0.95), #no rare variants included
  b_gl=rnorm(99), 
  rsq_prs=0.06, #h-squared value
  rsq_z=0.19, #to total 0.25 of total H-squared
  d_prev=5/100000, 
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
  b_u2d=0.1,
  b_gcc1=0.1)

dat$phen
```

    ## # A tibble: 100,000 × 11
    ##        u1     u2  r0[,1]  r1[,1]  c0[,1]  c1[,1] prs[,1] prs_w[,1]   l[,1]
    ##     <dbl>  <dbl>   <dbl>   <dbl>   <dbl>   <dbl>   <dbl>     <dbl>   <dbl>
    ##  1  1.03   2.67  -0.687  -0.410   0.0903 -0.488    0.330    0.238   0.866 
    ##  2 -1.21  -1.94  -1.61   -1.51   -0.417  -0.294   -1.06    -0.977  -0.282 
    ##  3  1.22   0.710 -0.0971  2.38   -0.558  -0.269    1.03     1.03    0.459 
    ##  4 -1.09   0.238 -0.412   0.495  -1.78    0.279   -1.84    -1.93   -0.956 
    ##  5  1.56   1.77   1.06    0.0414  0.210   0.475   -0.549   -0.642   0.242 
    ##  6  0.162 -0.498  0.324  -0.165   0.430  -0.297   -0.138   -0.0517 -0.0831
    ##  7 -0.960 -0.880 -2.15   -0.215  -0.0728 -0.0367   1.47     1.56   -0.163 
    ##  8  0.419  0.556  0.494  -0.861   1.05   -0.244    0.636    0.633   0.226 
    ##  9  0.882 -0.498 -0.492  -0.548  -0.338   1.19    -0.153   -0.0670 -0.434 
    ## 10 -0.259  0.410  0.0449  0.214   0.0713  1.10     0.132    0.129  -0.0334
    ## # ℹ 99,990 more rows
    ## # ℹ 2 more variables: prob_l <dbl[,1]>, d <dbl>

Variance and means of variables from data generating model and
correlation between PRS and disease liability

``` r
dgmodel_check(dat)
```

    ## # A tibble: 11 × 3
    ##    col        vars     means
    ##    <chr>     <dbl>     <dbl>
    ##  1 u1     0.999    -8.97e- 3
    ##  2 u2     0.993    -2.28e- 3
    ##  3 r0     0.939    -1.33e- 4
    ##  4 r1     0.931    -1.11e- 2
    ##  5 c0     0.896     1.14e- 3
    ##  6 c1     1.01     -5.68e- 3
    ##  7 prs    1         1.85e-17
    ##  8 prs_w  1         2.50e-17
    ##  9 l      0.272    -3.79e- 3
    ## 10 prob_l 0.000154  9.98e- 1
    ## 11 d      0.00155   1.55e- 3
    ##           [,1]
    ## [1,] 0.1792494

``` r
dgmodel_analysis(dat, ncase=1000, ncontrol=1000, protein_gwas=1000)
```

    ## # A tibble: 27 × 9
    ##           ahat     bhat       se    fval      pval      n x     y     study     
    ##          <dbl>    <dbl>    <dbl>   <dbl>     <dbl>  <int> <chr> <chr> <chr>     
    ##  1  0.00155    -0.00180 0.000124 210.    1.50e- 47 100000 prs   d     all       
    ##  2  0.00196    -0.534   0.0761    49.4   2.15e- 12 100000 c0    d     all       
    ##  3  0.000558   -0.446   0.0779    32.8   1.01e-  8 100000 r0    d     all       
    ##  4  0.00000962 -0.00847 0.00334    6.43  1.12e-  2 100000 prs   c0    all       
    ##  5  0.0000103   0.0770  0.00325  560.    1.72e-123 100000 prs   r0    all       
    ##  6 -0.00556    -0.0745  0.0807     0.852 3.56e-  1 100000 c1    d     all       
    ##  7 -0.0112      0.0819  0.0775     1.12  2.91e-  1 100000 r1    d     all       
    ##  8 -0.0000241  -0.00425 0.00315    1.82  1.77e-  1 100000 prs   c1    all       
    ##  9  0.0000295   0.00266 0.00328    0.661 4.16e-  1 100000 prs   r1    all       
    ## 10  0.117      -0.117   0.00864  184.    5.27e- 39   1155 prs   d     observati…
    ## # ℹ 17 more rows
