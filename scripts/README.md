# $\ell_0$ Trend Filtering

This file folder contains all the R scripts that are used to replicate the results in our paper. Before running these scripts, you need to install [R](https://cran.rstudio.com/) and [RStudio](https://posit.co/downloads/).


* SimuL0TF.Rmd: generate Figures 2-7 and include some more illustrative simulated examples. The output HTML file should look like [SimuL0TF.html](result/SimuL0TF.pdf).
* AlgoAnalysis.Rmd: replicate the reuslts and generate Figures 8-11 in Section 4.1.
* utils.R and amiasutils.R: source codes used in AlgoAnalysis.Rmd.
* RealData.R: replicate the results and generate all the graphs in Section 4.3.
* AlgoAnalysis_APP.Rmd: replicate the reuslts and generate Figures B.1-B.7 in Appendix B.1.
* [simu](simu) folder contains R scripts used in Appendix B.2: 
    * nsimu.R and tsimu.R: replicate the results for all methods except for the l0-MIP with large sample size.
    * nsimul0tfc.R and tsimul0tfc.R: replicate the results of the l0-MIP method when sample size is large.
* [simu_plots](simu_plots) folder contains the R Scripts used to generate Figures B.9-B.20 in Appendix B.2.
    * post_plot.R: generate Figures B.9, B.11 and B.13. 
    * pre_plot.R: generate Figures B.10, B.12 and B.14.
    * post_tplot.R: generate Figures B.15, B.17 and B.19. 
    * pre_tplot.R: generate Figures B.16, B.18 and B.20. 
    * combine_RData.R: combine the RData and needed to be run before generating all the figures.
