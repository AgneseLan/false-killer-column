# Age and sex differences in the vertebral column of the false killer whale (*Pseudorca crassidens*) 🐬🦴📈
### Testing changes to vertebrae shape and size during growth and between sexes in a population of stranded false killer whales using 3D geometric morphometrics 

Authors: [Mélissa Duflot](mailto:melissa.duflot.21@alumni.ucl.ac.uk?subject=[GitHub]%20Pseudorca%20Vertebrae%20Paper%20Code), Amadine Gillet,
Katrina Jones, Richard Sabin, [Agnese Lanzetti](mailto:agnese.lanzetti@gmail.com?subject=[GitHub]%20Pseudorca%20Vertebrae%20Paper%20Code)

To cite the paper: 

Available at: https://github.com/AgneseLan/false-killer-column

If using any of this code or data please cite the paper above and this repo

To cite this repo: 


## Data :floppy_disk: 

The data are provided in the Data folder. Mesh files (PLY) needed to test postioning when importing landmarks are available at following the provided OneDrive link.

- __Landmark data__: *pts folder* <br />
Text files with landmark coordinates for each specimen in PTS format. Unzip folder first.

- __Surface data__: *ply folder* <br />
Empty folder where mesh files from OneDrive need to be saved to reproduce the code. Unzip folder first.

- __Specimens' classifiers, landmark/curves lists__: *curves_caudal.csv, curves_lumbar.csv, curves_thoracic.csv, info_caudal.csv, info_lumbar.csv, info_thoracic.csv, landmark_caudal.csv, landmark_lumbar.csv, landmark_thoracic.csv* <br />
Spreadsheets with additional inforation for analyses: list of curves for each vertebra type, list of landmarks for each vertebra type, classifiers for specimens (vertebra number, vertebra ID, specimen number,	sex,	size (total length, cm),	age) for each vertebra type.

## Analysis :computer:
In this repository you will find raw data (.csv and data files) and code for analyses (code supplied as .R files)

📁 Data

As described above. Meshes used to collect and test landmarks available following provided OneDrive link. 

⌨ Code for analyses - .R files

*1-import_caudal.R, 1-import_lumbar.R, 1-import_thoracic.R, 1-slider3d_2.r, 2-morphoblocks-gpa_pca_caudal.R, 2-morphoblocks-gpa_pca_lumbar.R, 2-morphoblocks-gpa_pca_thoracic.R, 3-morphoblocks-pca.r, 4-morphoblocks-allometry.r, 5-morphoblocks-phenotypic_trajectory.R*

Code files are numbered providing the order the analyses need to be performed in.
Before running analyses, save Data folder in the same directory as the R project. This will allow to import the data as detailed in the code provided.

## License 📃
This project is licensed under the MIT License - see the LICENSE.md file for details

## Session Info 📋
For reproducibility purposes, here is the output of sessioninfo::session_info() used to perform the analyses in the publication.

```
- Session info -----------------------------------------------------------------------------------------------------------
setting  value
 version  R version 4.1.3 (2022-03-10)
 os       Windows 10 x64 (build 19045)
 system   x86_64, mingw32
 ui       RStudio
 language (EN)
 collate  English_United States.1252
 ctype    English_United States.1252
 tz       Europe/London
 date     2023-07-19
 rstudio  2023.06.0+421 Mountain Hydrangea (desktop)
 pandoc   3.1.1 @ C:/Program Files/RStudio/resources/app/bin/quarto/bin/tools/ (via rmarkdown)

- Packages ---------------------------------------------------------------------------------------------------------------
 package           * version    date (UTC) lib source
 abind             * 1.4-5      2016-07-21 [1] CRAN (R 4.1.1)
 ade4                1.7-22     2023-02-06 [1] CRAN (R 4.1.3)
 adegenet            2.1.10     2023-01-26 [1] CRAN (R 4.1.3)
 adephylo            1.1-13     2022-10-19 [1] CRAN (R 4.1.3)
 ape               * 5.7-1      2023-03-13 [1] CRAN (R 4.1.3)
 backports           1.4.1      2021-12-13 [1] CRAN (R 4.1.2)
 base64enc           0.1-3      2015-07-28 [1] CRAN (R 4.1.0)
 bezier              1.1.2      2018-12-14 [1] CRAN (R 4.1.1)
 borealis          * 2021.02.04 2021-02-10 [1] Github (aphanotus/borealis@cd96174)
 broom               1.0.5      2023-06-09 [1] CRAN (R 4.1.3)
 cachem              1.0.8      2023-05-01 [1] CRAN (R 4.1.3)
 callr               3.7.3      2022-11-02 [1] CRAN (R 4.1.3)
 car               * 3.1-2      2023-03-30 [1] CRAN (R 4.1.3)
 carData           * 3.0-5      2022-01-06 [1] CRAN (R 4.1.2)
 caret               6.0-94     2023-03-21 [1] CRAN (R 4.1.3)
 checkmate           2.2.0      2023-04-27 [1] CRAN (R 4.1.3)
 class               7.3-22     2023-05-03 [1] CRAN (R 4.1.3)
 cli               * 3.6.1      2023-03-23 [1] CRAN (R 4.1.3)
 cluster             2.1.4      2022-08-22 [1] CRAN (R 4.1.3)
 clusterGeneration   1.3.7      2020-12-15 [1] CRAN (R 4.1.1)
 coda                0.19-4     2020-09-30 [1] CRAN (R 4.1.1)
 codetools           0.2-18     2020-11-04 [2] CRAN (R 4.1.3)
 colorRamps          2.3.1      2022-05-02 [1] CRAN (R 4.1.3)
 colorspace          2.1-0      2023-01-23 [1] CRAN (R 4.1.3)
 combinat            0.0-8      2012-10-29 [1] CRAN (R 4.1.1)
 corpcor             1.6.10     2021-09-16 [1] CRAN (R 4.1.1)
 crayon              1.5.2      2022-09-29 [1] CRAN (R 4.1.3)
 curl                5.0.1      2023-06-07 [1] CRAN (R 4.1.3)
 data.table          1.14.8     2023-02-17 [1] CRAN (R 4.1.3)
 Deriv               4.1.3      2021-02-24 [1] CRAN (R 4.1.3)
 deSolve             1.35       2023-03-12 [1] CRAN (R 4.1.3)
 devtools          * 2.4.5      2022-10-11 [1] CRAN (R 4.1.3)
 digest              0.6.32     2023-06-26 [1] CRAN (R 4.1.3)
 doParallel          1.0.17     2022-02-07 [1] CRAN (R 4.1.2)
 dplyr             * 1.1.2      2023-04-20 [1] CRAN (R 4.1.3)
 ellipsis            0.3.2      2021-04-29 [1] CRAN (R 4.1.0)
 EMMLi             * 0.0.3      2017-02-17 [1] CRAN (R 4.1.1)
 evaluate            0.21       2023-05-05 [1] CRAN (R 4.1.3)
 evomap            * 0.0.0.9000 2021-04-15 [1] Github (JeroenSmaers/evomap@dfa7dfd)
 expm                0.999-7    2023-01-09 [1] CRAN (R 4.1.3)
 fansi               1.0.4      2023-01-22 [1] CRAN (R 4.1.3)
 fastmap             1.1.1      2023-02-24 [1] CRAN (R 4.1.3)
 fastmatch           1.1-3      2021-07-23 [1] CRAN (R 4.1.0)
 fdrtool             1.2.17     2021-11-13 [1] CRAN (R 4.1.2)
 forcats           * 1.0.0      2023-01-29 [1] CRAN (R 4.1.3)
 foreach             1.5.2      2022-02-02 [1] CRAN (R 4.1.2)
 foreign             0.8-84     2022-12-06 [1] CRAN (R 4.1.3)
```
