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
rmarkdown::render("sims_functions.Rmd")
```

    ## 
    ## 
    ## processing file: sims_functions.Rmd

    ##   |                                                            |                                                    |   0%  |                                                            |.....                                               |   9%                    |                                                            |.........                                           |  18% [unnamed-chunk-5]  |                                                            |..............                                      |  27%                    |                                                            |...................                                 |  36% [unnamed-chunk-6]  |                                                            |........................                            |  45%                    |                                                            |............................                        |  55% [unnamed-chunk-7]  |                                                            |.................................                   |  64%                    |                                                            |......................................              |  73% [unnamed-chunk-8]  |                                                            |...........................................         |  82%                    |                                                            |...............................................     |  91% [unnamed-chunk-9]  |                                                            |....................................................| 100%                  

    ## output file: sims_functions.knit.md

    ## /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/pandoc +RTS -K512m -RTS sims_functions.knit.md --to gfm+tex_math_dollars-yaml_metadata_block --from markdown+autolink_bare_uris+tex_math_single_backslash --output sims_functions.md --template /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown/rmarkdown/templates/github_document/resources/default.md 
    ## /Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/pandoc +RTS -K512m -RTS sims_functions.md --to html4 --from gfm+tex_math_dollars --output sims_functions.html --embed-resources --standalone --highlight-style pygments --template /Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown/rmarkdown/templates/github_document/resources/preview.html --variable 'github-markdown-css:/Library/Frameworks/R.framework/Versions/4.2/Resources/library/rmarkdown/rmarkdown/templates/github_document/resources/github.css' --metadata pagetitle=PREVIEW --mathjax

    ## 
    ## Preview created: /var/folders/ml/4j6q26qn7r5b89gy3dssz_d40000gn/T//RtmpWpRxxF/preview-a669d83a788.html

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
    ##         u1      u2 r0[,1]  r1[,1]   c0[,1] c1[,1] prs[,1] prs_w[,1]  l[,1]
    ##      <dbl>   <dbl>  <dbl>   <dbl>    <dbl>  <dbl>   <dbl>     <dbl>  <dbl>
    ##  1  0.130   1.76    0.504  0.737   0.448   -0.769 -0.791    -0.778   0.422
    ##  2 -1.24   -1.10   -0.719 -0.0807 -1.08     0.360 -0.252    -0.239  -0.979
    ##  3 -0.668   0.0975  0.579  0.464   0.159   -1.98  -0.0610   -0.0484  0.360
    ##  4  2.40    0.475  -0.330  1.36   -1.47     0.781 -0.483    -0.470  -0.210
    ##  5 -1.37    1.56    0.532 -0.642  -0.664   -1.62   0.824     0.823  -0.457
    ##  6 -0.300  -0.275  -1.08  -0.749   0.274    0.317 -1.73     -1.75   -1.27 
    ##  7  0.344  -1.06    0.233  0.306   0.630    0.799  0.802     0.814  -0.656
    ##  8  1.05   -1.32    0.586 -1.77    0.00892  0.705 -0.844    -0.831  -0.273
    ##  9 -0.0107 -0.413   0.237  0.264   0.358    1.21   1.03      1.04    0.104
    ## 10  1.22   -1.41   -0.360  0.117   1.28     1.71   1.44      1.45    0.837
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
    ##  1 u1     1.00      2.35e- 3
    ##  2 u2     1.00     -2.67e- 3
    ##  3 r0     0.928    -2.35e- 3
    ##  4 r1     0.935    -3.28e- 3
    ##  5 c0     0.894     1.05e- 3
    ##  6 c1     1.01     -2.82e- 3
    ##  7 prs    1         2.41e-17
    ##  8 prs_w  1         6.22e-17
    ##  9 l      0.270    -1.72e- 4
    ## 10 prob_l 0.000150  9.98e- 1
    ## 11 d      0.00145   1.45e- 3
    ##           [,1]
    ## [1,] 0.1849627

``` r
dgmodel_analysis(dat, ncase=1000, ncontrol=1000, protein_gwas=1000)
```

    ## # A tibble: 27 × 9
    ##           ahat     bhat       se     fval      pval      n x     y     study    
    ##          <dbl>    <dbl>    <dbl>    <dbl>     <dbl>  <int> <chr> <chr> <chr>    
    ##  1  0.00145    -0.00173 0.000120 206.     1.11e- 46 100000 prs   d     all      
    ##  2  0.00175    -0.489   0.0786    38.7    5.08e- 10 100000 c0    d     all      
    ##  3 -0.00162    -0.508   0.0800    40.2    2.27e- 10 100000 r0    d     all      
    ##  4 -0.00000723  0.00691 0.00334    4.28   3.86e-  2 100000 prs   c0    all      
    ##  5  0.000177    0.0754  0.00327  530.     4.72e-117 100000 prs   r0    all      
    ##  6 -0.00284     0.00931 0.0835     0.0125 9.11e-  1 100000 c1    d     all      
    ##  7 -0.00314    -0.0959  0.0803     1.43   2.32e-  1 100000 r1    d     all      
    ##  8  0.0000106   0.00376 0.00315    1.43   2.32e-  1 100000 prs   c1    all      
    ##  9  0.0000245   0.00746 0.00327    5.20   2.26e-  2 100000 prs   r1    all      
    ## 10  0.111      -0.119   0.00858  194.     7.52e- 41   1145 prs   d     observat…
    ## # ℹ 17 more rows
