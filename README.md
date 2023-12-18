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
For reproducibility purposes, here is the output of utils::sessionInfo() used to perform the analyses in the publication.

R version 4.3.2 (2023-10-31 ucrt)
Platform: x86_64-w64-mingw32/x64 (64-bit)
Running under: Windows 10 x64 (build 19045)

attached base packages:
[1] stats     graphics  grDevices utils     datasets  methods   base     

other attached packages:
 [1] borealis_2021.02.04 morphoBlocks_0.1.0  scales_1.2.1        reshape2_1.4.4      car_3.1-2           carData_3.0-5      
 [7] evomap_0.0.0.9000   gridExtra_2.3       png_0.1-8           rphylopic_1.2.2     ggplotify_0.1.2     ggpubr_0.6.0       
[13] ggthemes_4.2.4      RColorBrewer_1.1-3  ggfortify_0.4.16    ggphylomorpho_0.1.0 gginnards_0.1.2     ggrepel_0.9.4      
[19] rlang_1.1.2         cli_3.6.1           magick_2.8.1        devtools_2.4.5      usethis_2.2.2       abind_1.4-5        
[25] geiger_2.0.11       phytools_2.0-3      maps_3.4.1.1        ape_5.7-1           qgraph_1.9.8        EMMLi_0.0.3        
[31] paleomorph_0.1.4    Rvcg_0.22.1         geomorph_4.0.6      Matrix_1.6-3        rgl_1.2.1           RRPP_1.4.0         
[37] Morpho_2.11         lubridate_1.9.3     forcats_1.0.0       stringr_1.5.1       dplyr_1.1.4         purrr_1.0.2        
[43] readr_2.1.4         tidyr_1.3.0         tibble_3.2.1        ggplot2_3.4.4       tidyverse_2.0.0    

loaded via a namespace (and not attached):
  [1] splines_4.3.2           later_1.3.1             XML_3.99-0.15           rpart_4.1.21           
  [5] deSolve_1.38            lifecycle_1.0.4         rstatix_0.7.2           doParallel_1.0.17      
  [9] processx_3.8.2          lattice_0.22-5          MASS_7.3-60             backports_1.4.1        
 [13] magrittr_2.0.3          Hmisc_5.1-1             rmarkdown_2.25          remotes_2.4.2.1        
 [17] httpuv_1.6.12           grImport2_0.3-1         sessioninfo_1.2.2       pkgbuild_1.4.2         
 [21] pbapply_1.7-2           ade4_1.7-22             pkgload_1.3.3           expm_0.999-7           
 [25] quadprog_1.5-8          yulab.utils_0.1.0       nnet_7.3-19             vegan_2.6-4            
 [29] adegenet_2.1.10         permute_0.9-7           codetools_0.2-19        xml2_1.3.5             
 [33] adephylo_1.1-16         tidyselect_1.2.0        RNeXML_2.4.11           stats4_4.3.2           
 [37] base64enc_0.1-3         jsonlite_1.8.7          ellipsis_0.3.2          phylobase_0.8.10       
 [41] Formula_1.2-5           iterators_1.0.14        foreach_1.5.2           tools_4.3.2            
 [45] progress_1.2.2          Rcpp_1.0.11             glue_1.6.2              mnormt_2.1.1           
 [49] xfun_0.41               mgcv_1.9-0              withr_2.5.2             numDeriv_2016.8-1.1    
 [53] combinat_0.0-8          fastmap_1.1.1           fansi_1.0.5             callr_3.7.3            
 [57] digest_0.6.33           gridGraphics_0.5-1      timechange_0.2.0        R6_2.5.1               
 [61] mime_0.12               RGCCA_2.1.2             colorspace_2.1-0        rsvg_2.6.0             
 [65] gtools_3.9.4            jpeg_0.1-10             utf8_1.2.4              generics_0.1.3         
 [69] data.table_1.14.8       corpcor_1.6.10          clusterGeneration_1.3.8 prettyunits_1.2.0      
 [73] httr_1.4.7              htmlwidgets_1.6.2       scatterplot3d_0.3-44    pkgconfig_2.0.3        
 [77] gtable_0.3.4            htmltools_0.5.7         lavaan_0.6-16           profvis_0.3.8          
 [81] colorRamps_2.3.1        knitr_1.45              rstudioapi_0.15.0       tzdb_0.4.0             
 [85] rncl_0.8.7              uuid_1.1-1              curl_5.1.0              coda_0.19-4            
 [89] checkmate_2.3.0         nlme_3.1-163            cachem_1.0.8            parallel_4.3.2         
 [93] miniUI_0.1.1.1          foreign_0.8-85          pillar_1.9.0            grid_4.3.2             
 [97] vctrs_0.6.4             urlchecker_1.0.1        promises_1.2.1          xtable_1.8-4           
[101] cluster_2.1.4           htmlTable_2.4.2         evaluate_0.23           pbivnorm_0.6.0         
[105] mvtnorm_1.2-3           compiler_4.3.2          bezier_1.1.2            crayon_1.5.2           
[109] ggsignif_0.6.4          fdrtool_1.2.17          ps_1.7.5                plyr_1.8.9             
[113] fs_1.6.3                stringi_1.8.1           psych_2.3.9             subplex_1.8            
[117] munsell_0.5.0           optimParallel_1.0-2     hms_1.1.3               glasso_1.11            
[121] seqinr_4.2-30           shiny_1.8.0             broom_1.0.5             igraph_1.5.1           
[125] memoise_2.0.1           phangorn_2.11.1         fastmatch_1.1-4   
 
