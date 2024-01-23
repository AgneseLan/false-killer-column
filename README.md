# Age and sex differences in the vertebral column of the false killer whale (*Pseudorca crassidens*) 🐬🦴📈
### Testing changes to vertebrae shape and size during growth and between sexes in a population of stranded false killer whales using 3D geometric morphometrics 

Authors: [Mélissa Duflot](mailto:melissa.duflot.21@alumni.ucl.ac.uk?subject=[GitHub]%20Pseudorca%20Vertebrae%20Paper%20Code), Amadine Gillet,
Katrina Jones, Richard Sabin, [Agnese Lanzetti](mailto:agnese.lanzetti@gmail.com?subject=[GitHub]%20Pseudorca%20Vertebrae%20Paper%20Code)

To cite the paper: 

Available at: https://github.com/AgneseLan/false-killer-column

If using any of this code or data please cite the paper above and this repo

To cite this repo: 


## Data :floppy_disk: 

The data are provided in the Data folder. Mesh files (PLY) needed to test postioning when importing landmarks are available at the provided Zenodo link.

- __Landmark data__: *pts folder* <br />
Text files with landmark coordinates for each specimen in PTS format. Unzip folder first.

- __Surface data__: *ply folder* <br />
Empty folder where mesh files from Zenodo need to be saved to reproduce the code. Folder contains txt file with the Zenodo link.

- __Specimens' classifiers, landmark/curves lists__: *curves_caudal.csv, curves_lumbar.csv, curves_thoracic.csv, info_caudal.csv, info_lumbar.csv, info_thoracic.csv, landmark_caudal.csv, landmark_lumbar.csv, landmark_thoracic.csv* <br />
Spreadsheets with additional inforation for analyses: list of curves for each vertebra type, list of landmarks for each vertebra type, classifiers for specimens (vertebra number, vertebra ID, specimen number,	sex,	size (total length, cm),	age) for each vertebra type.

## Analysis :computer:
In this repository you will find raw data (.csv and data files) and code for analyses (code supplied as .R files)

📁 Data

As described above. Meshes used to collect and test landmarks available following the provided Zenodo link. 

⌨ Code for analyses - .R files

*1-import_caudal.R, 1-import_lumbar.R, 1-import_thoracic.R, 1-slider3d_2.r, 2-gpa_pca_caudal.R, 2-gpa_pca_lumbar.R, 2-gpa_pca_thoracic.R, 3-morphoblocks-pca.r, 4-morphoblocks-allometry.r, 5-morphoblocks-phenotypic_trajectory.R*

Code files are numbered providing the order the analyses need to be performed in.
Before running analyses, save Data folder in the same directory as the R project. This will allow to import the data as detailed in the code provided.

## License 📃
This project is licensed under the MIT License - see the LICENSE.md file for details

## Session Info 📋
For reproducibility purposes, here is the output of sessioninfo::session_info() used to perform the analyses in the publication.

```
─ Session info ───────────────────────────────────────────────────────────────────────────────────────────────────────────
 setting  value
 version  R version 4.3.2 (2023-10-31 ucrt)
 os       Windows 10 x64 (build 19045)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United States.utf8
 ctype    English_United States.utf8
 tz       Europe/London
 date     2023-12-18
 rstudio  2023.09.1+494 Desert Sunflower (desktop)
 pandoc   3.1.1 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

─ Packages ───────────────────────────────────────────────────────────────────────────────────────────────────────────────
 package           * version    date (UTC) lib source
 abind             * 1.4-5      2016-07-21 [2] CRAN (R 4.3.1)
 ade4                1.7-22     2023-02-06 [2] CRAN (R 4.3.2)
 adegenet            2.1.10     2023-01-26 [2] CRAN (R 4.3.2)
 adephylo            1.1-16     2023-10-06 [2] CRAN (R 4.3.2)
 ape               * 5.7-1      2023-03-13 [2] CRAN (R 4.3.2)
 backports           1.4.1      2021-12-13 [2] CRAN (R 4.3.1)
 base64enc           0.1-3      2015-07-28 [2] CRAN (R 4.3.1)
 bezier              1.1.2      2018-12-14 [2] CRAN (R 4.3.1)
 borealis          * 2021.02.04 2021-02-10 [2] Github (aphanotus/borealis@cd96174)
 broom               1.0.5      2023-06-09 [2] CRAN (R 4.3.2)
 cachem              1.0.8      2023-05-01 [2] CRAN (R 4.3.2)
 callr               3.7.3      2022-11-02 [2] CRAN (R 4.3.2)
 car               * 3.1-2      2023-03-30 [2] CRAN (R 4.3.2)
 carData           * 3.0-5      2022-01-06 [2] CRAN (R 4.3.2)
 checkmate           2.3.0      2023-10-25 [2] CRAN (R 4.3.2)
 cli               * 3.6.1      2023-03-23 [1] CRAN (R 4.3.2)
 cluster             2.1.4      2022-08-22 [2] CRAN (R 4.3.2)
 clusterGeneration   1.3.8      2023-08-16 [2] CRAN (R 4.3.2)
 coda                0.19-4     2020-09-30 [2] CRAN (R 4.3.2)
 codetools           0.2-19     2023-02-01 [3] CRAN (R 4.3.2)
 colorRamps          2.3.1      2022-05-02 [2] CRAN (R 4.3.1)
 colorspace          2.1-0      2023-01-23 [2] CRAN (R 4.3.2)
 combinat            0.0-8      2012-10-29 [2] CRAN (R 4.3.1)
 corpcor             1.6.10     2021-09-16 [2] CRAN (R 4.3.1)
 crayon              1.5.2      2022-09-29 [2] CRAN (R 4.3.2)
 curl                5.1.0      2023-10-02 [2] CRAN (R 4.3.2)
 data.table          1.14.8     2023-02-17 [2] CRAN (R 4.3.2)
 deSolve             1.38       2023-09-05 [2] CRAN (R 4.3.2)
 devtools          * 2.4.5      2022-10-11 [2] CRAN (R 4.3.2)
 digest              0.6.33     2023-07-07 [2] CRAN (R 4.3.2)
 doParallel          1.0.17     2022-02-07 [2] CRAN (R 4.3.2)
 dplyr             * 1.1.4      2023-11-17 [2] CRAN (R 4.3.2)
 ellipsis            0.3.2      2021-04-29 [2] CRAN (R 4.3.2)
 EMMLi             * 0.0.3      2017-02-17 [2] CRAN (R 4.3.2)
 evaluate            0.23       2023-11-01 [2] CRAN (R 4.3.2)
 evomap            * 0.0.0.9000 2021-04-15 [2] Github (JeroenSmaers/evomap@dfa7dfd)
 expm                0.999-7    2023-01-09 [2] CRAN (R 4.3.2)
 fansi               1.0.5      2023-10-08 [2] CRAN (R 4.3.2)
 fastmap             1.1.1      2023-02-24 [2] CRAN (R 4.3.2)
 fastmatch           1.1-4      2023-08-18 [2] CRAN (R 4.3.1)
 fdrtool             1.2.17     2021-11-13 [2] CRAN (R 4.3.1)
 forcats           * 1.0.0      2023-01-29 [2] CRAN (R 4.3.2)
 foreach             1.5.2      2022-02-02 [2] CRAN (R 4.3.2)
 foreign             0.8-85     2023-09-09 [2] CRAN (R 4.3.1)
 Formula             1.2-5      2023-02-24 [2] CRAN (R 4.3.1)
 fs                  1.6.3      2023-07-20 [2] CRAN (R 4.3.2)
 geiger            * 2.0.11     2023-04-03 [2] CRAN (R 4.3.2)
 generics            0.1.3      2022-07-05 [2] CRAN (R 4.3.2)
 geomorph          * 4.0.6      2023-08-31 [2] CRAN (R 4.3.2)
 ggfortify         * 0.4.16     2023-03-20 [2] CRAN (R 4.3.2)
 gginnards         * 0.1.2      2023-05-24 [2] CRAN (R 4.3.2)
 ggphylomorpho     * 0.1.0      2021-02-03 [2] Github (wabarr/ggphylomorpho@7e1228b)
 ggplot2           * 3.4.4      2023-10-12 [2] CRAN (R 4.3.2)
 ggplotify         * 0.1.2      2023-08-09 [2] CRAN (R 4.3.2)
 ggpubr            * 0.6.0      2023-02-10 [2] CRAN (R 4.3.2)
 ggrepel           * 0.9.4      2023-10-13 [2] CRAN (R 4.3.2)
 ggsignif            0.6.4      2022-10-13 [2] CRAN (R 4.3.2)
 ggthemes          * 4.2.4      2021-01-20 [2] CRAN (R 4.3.2)
 glasso              1.11       2019-10-01 [2] CRAN (R 4.3.1)
 glue                1.6.2      2022-02-24 [1] CRAN (R 4.3.2)
 gridExtra         * 2.3        2017-09-09 [2] CRAN (R 4.3.2)
 gridGraphics        0.5-1      2020-12-13 [2] CRAN (R 4.3.2)
 grImport2           0.3-1      2023-10-27 [2] CRAN (R 4.3.2)
 gtable              0.3.4      2023-08-21 [2] CRAN (R 4.3.2)
 gtools              3.9.4      2022-11-27 [2] CRAN (R 4.3.2)
 Hmisc               5.1-1      2023-09-12 [2] CRAN (R 4.3.2)
 hms                 1.1.3      2023-03-21 [2] CRAN (R 4.3.2)
 htmlTable           2.4.2      2023-10-29 [2] CRAN (R 4.3.2)
 htmltools           0.5.7      2023-11-03 [2] CRAN (R 4.3.2)
 htmlwidgets         1.6.2      2023-03-17 [2] CRAN (R 4.3.2)
 httpuv              1.6.12     2023-10-23 [2] CRAN (R 4.3.2)
 httr                1.4.7      2023-08-15 [2] CRAN (R 4.3.2)
 igraph              1.5.1      2023-08-10 [2] CRAN (R 4.3.2)
 iterators           1.0.14     2022-02-05 [2] CRAN (R 4.3.2)
 jpeg                0.1-10     2022-11-29 [2] CRAN (R 4.3.1)
 jsonlite            1.8.7      2023-06-29 [2] CRAN (R 4.3.2)
 knitr               1.45       2023-10-30 [2] CRAN (R 4.3.2)
 later               1.3.1      2023-05-02 [2] CRAN (R 4.3.2)
 lattice             0.22-5     2023-10-24 [1] CRAN (R 4.3.2)
 lavaan              0.6-16     2023-07-19 [2] CRAN (R 4.3.2)
 lifecycle           1.0.4      2023-11-07 [1] CRAN (R 4.3.2)
 lubridate         * 1.9.3      2023-09-27 [2] CRAN (R 4.3.2)
 magick            * 2.8.1      2023-10-22 [2] CRAN (R 4.3.2)
 magrittr            2.0.3      2022-03-30 [1] CRAN (R 4.3.2)
 maps              * 3.4.1.1    2023-11-03 [2] CRAN (R 4.3.2)
 MASS                7.3-60     2023-05-04 [2] CRAN (R 4.3.2)
 Matrix            * 1.6-3      2023-11-14 [2] CRAN (R 4.3.2)
 memoise             2.0.1      2021-11-26 [2] CRAN (R 4.3.2)
 mgcv                1.9-0      2023-07-11 [2] CRAN (R 4.3.2)
 mime                0.12       2021-09-28 [2] CRAN (R 4.3.1)
 miniUI              0.1.1.1    2018-05-18 [2] CRAN (R 4.3.2)
 mnormt              2.1.1      2022-09-26 [2] CRAN (R 4.3.1)
 Morpho            * 2.11       2023-01-27 [2] CRAN (R 4.3.2)
 morphoBlocks      * 0.1.0      2023-11-20 [1] Github (aharmer/morphoBlocks@00c2051)
 munsell             0.5.0      2018-06-12 [2] CRAN (R 4.3.2)
 mvtnorm             1.2-3      2023-08-25 [2] CRAN (R 4.3.2)
 nlme                3.1-163    2023-08-09 [2] CRAN (R 4.3.2)
 nnet                7.3-19     2023-05-03 [2] CRAN (R 4.3.2)
 numDeriv            2016.8-1.1 2019-06-06 [2] CRAN (R 4.3.1)
 optimParallel       1.0-2      2021-02-11 [2] CRAN (R 4.3.2)
 paleomorph        * 0.1.4      2017-04-19 [2] CRAN (R 4.3.2)
 pbapply             1.7-2      2023-06-27 [2] CRAN (R 4.3.2)
 pbivnorm            0.6.0      2015-01-23 [2] CRAN (R 4.3.1)
 permute             0.9-7      2022-01-27 [2] CRAN (R 4.3.2)
 phangorn            2.11.1     2023-01-23 [2] CRAN (R 4.3.2)
 phylobase           0.8.10     2020-03-01 [2] CRAN (R 4.3.2)
 phytools          * 2.0-3      2023-11-09 [2] CRAN (R 4.3.2)
 pillar              1.9.0      2023-03-22 [2] CRAN (R 4.3.2)
 pkgbuild            1.4.2      2023-06-26 [2] CRAN (R 4.3.2)
 pkgconfig           2.0.3      2019-09-22 [2] CRAN (R 4.3.2)
 pkgload             1.3.3      2023-09-22 [2] CRAN (R 4.3.2)
 plyr                1.8.9      2023-10-02 [2] CRAN (R 4.3.2)
 png               * 0.1-8      2022-11-29 [2] CRAN (R 4.3.1)
 prettyunits         1.2.0      2023-09-24 [2] CRAN (R 4.3.2)
 processx            3.8.2      2023-06-30 [2] CRAN (R 4.3.2)
 profvis             0.3.8      2023-05-02 [2] CRAN (R 4.3.2)
 progress            1.2.2      2019-05-16 [2] CRAN (R 4.3.2)
 promises            1.2.1      2023-08-10 [2] CRAN (R 4.3.2)
 ps                  1.7.5      2023-04-18 [2] CRAN (R 4.3.2)
 psych               2.3.9      2023-09-26 [2] CRAN (R 4.3.2)
 purrr             * 1.0.2      2023-08-10 [2] CRAN (R 4.3.2)
 qgraph            * 1.9.8      2023-11-03 [2] CRAN (R 4.3.2)
 quadprog            1.5-8      2019-11-20 [2] CRAN (R 4.3.1)
 R6                  2.5.1      2021-08-19 [2] CRAN (R 4.3.2)
 RColorBrewer      * 1.1-3      2022-04-03 [2] CRAN (R 4.3.1)
 Rcpp                1.0.11     2023-07-06 [2] CRAN (R 4.3.2)
 readr             * 2.1.4      2023-02-10 [2] CRAN (R 4.3.2)
 remotes             2.4.2.1    2023-07-18 [2] CRAN (R 4.3.2)
 reshape2          * 1.4.4      2020-04-09 [2] CRAN (R 4.3.2)
 RGCCA               2.1.2      2017-05-11 [1] CRAN (R 4.3.2)
 rgl               * 1.2.1      2023-07-06 [2] CRAN (R 4.3.2)
 rlang             * 1.1.2      2023-11-04 [1] CRAN (R 4.3.2)
 rmarkdown           2.25       2023-09-18 [2] CRAN (R 4.3.2)
 rncl                0.8.7      2023-01-08 [2] CRAN (R 4.3.2)
 RNeXML              2.4.11     2023-02-01 [2] CRAN (R 4.3.2)
 rpart               4.1.21     2023-10-09 [2] CRAN (R 4.3.2)
 rphylopic         * 1.2.2      2023-10-28 [2] CRAN (R 4.3.2)
 RRPP              * 1.4.0      2023-08-15 [2] CRAN (R 4.3.2)
 rstatix             0.7.2      2023-02-01 [2] CRAN (R 4.3.2)
 rstudioapi          0.15.0     2023-07-07 [2] CRAN (R 4.3.2)
 rsvg                2.6.0      2023-10-08 [2] CRAN (R 4.3.2)
 Rvcg              * 0.22.1     2023-01-26 [2] CRAN (R 4.3.2)
 scales            * 1.2.1      2022-08-20 [2] CRAN (R 4.3.2)
 scatterplot3d       0.3-44     2023-05-05 [2] CRAN (R 4.3.1)
 seqinr              4.2-30     2023-04-05 [2] CRAN (R 4.3.2)
 sessioninfo         1.2.2      2021-12-06 [2] CRAN (R 4.3.2)
 shiny               1.8.0      2023-11-17 [2] CRAN (R 4.3.2)
 stringi             1.8.1      2023-11-13 [1] CRAN (R 4.3.2)
 stringr           * 1.5.1      2023-11-14 [1] CRAN (R 4.3.2)
 subplex             1.8        2022-04-12 [2] CRAN (R 4.3.1)
 tibble            * 3.2.1      2023-03-20 [2] CRAN (R 4.3.2)
 tidyr             * 1.3.0      2023-01-24 [2] CRAN (R 4.3.2)
 tidyselect          1.2.0      2022-10-10 [2] CRAN (R 4.3.2)
 tidyverse         * 2.0.0      2023-02-22 [2] CRAN (R 4.3.2)
 timechange          0.2.0      2023-01-11 [2] CRAN (R 4.3.2)
 tzdb                0.4.0      2023-05-12 [2] CRAN (R 4.3.2)
 urlchecker          1.0.1      2021-11-30 [2] CRAN (R 4.3.2)
 usethis           * 2.2.2      2023-07-06 [2] CRAN (R 4.3.2)
 utf8                1.2.4      2023-10-22 [2] CRAN (R 4.3.2)
 uuid                1.1-1      2023-08-17 [2] CRAN (R 4.3.1)
 vctrs               0.6.4      2023-10-12 [1] CRAN (R 4.3.2)
 vegan               2.6-4      2022-10-11 [2] CRAN (R 4.3.2)
 withr               2.5.2      2023-10-30 [2] CRAN (R 4.3.2)
 xfun                0.41       2023-11-01 [2] CRAN (R 4.3.2)
 XML                 3.99-0.15  2023-11-02 [2] CRAN (R 4.3.2)
 xml2                1.3.5      2023-07-06 [2] CRAN (R 4.3.2)
 xtable              1.8-4      2019-04-21 [2] CRAN (R 4.3.2)
 yulab.utils         0.1.0      2023-09-20 [2] CRAN (R 4.3.2)

```
