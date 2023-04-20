# AMIAS: an R package for solving the generalized $\ell_0$ problem

`AMIAS` (Alternating Minimization Induced Active Set) is an R package that aims to solve the generalized $\ell_0$ problem. That is, given an observed vector $y$, the generalized $\ell_0$ problem minimizer the sum of squared residuals with $\ell_0$-regularization on some linear combination of the estimator. The AMIAS method is based on the necessary optimality conditions derived from an augmented Lagrangian framework. The proposed method takes full advantage of the primal and dual variables with complementary supports, and decouples the high-dimensional problem into two sub-systems on the active and inactive sets, respectively. A sequential AMIAS algorithm with warm start initialization is developed for efficient determination of the cardinality parameter, along with the output of solution paths.

The most typical example is the $\ell_0$ trend filtering, which is an effective tool for nonparametric regression with the power of automatic knot detection
in function values or derivatives. For more details, please see the paper **$\ell_0$ Trend Filtering** by C. Wen and X. Wang and A. Zhang. 


## Installation

To install the `AMIAS` R package from Github, just run:
```r
if(!require(devtools)) install.packages('devtools')
devtools::install_github("INFORMSJoC/2021.0313/scripts/AMIAS")
```
or

```r
if(!require(devtools)) install.packages('devtools')
devtools::install_github("C2S2-HF/L0TF/scripts/AMIAS")
```

## File Folder Structure

This folder contains the source files for implementing the AMIAS algorithm:

- [R](R) folder includes the R functions.
    - [amias.R](amias.R) is the main function to implement the AMIAS algorithm. The following three files, i.e., [coef.amias.R](coef.amias.R), [plot.amias.R](plot.amias.R), [print.amias.R](print.amias.R), are used to extract the estimator, plot the estimator and print the results from the AMIAS algorithm.
    - [samias.R](samias.R) is the main function sequential AMIAS algorithm. The following three files, i.e., [coef.samias.R](coef.samias.R), [plot.samias.R](plot.samias.R), [print.samias.R](print.samias.R), are used to extract the estimator, plot the estimator and print the results from the AMIAS algorithm.
    - The files prefixed with `sim`, i.e., [sim.blocks.R](sim.blocks.R), [sim.wave.R](sim.wave.R), [sim.doppler.R](sim.doppler.R), are used to generate simulated data for demonstration, as illustrated in the Examples associated with the `amias` or `samias` function.
    - [genD.R](genD.R) provides functions to generate the difference matrix with any order.
    - [my.rollmean.R](my.rollmean.R) provide function for smooth the final estimator if necessary.
    - [inv2d.R](inv2d.R) and [getM.R](getM.R) contain functions that are called by the main functions and unavailable outside of the AMIAS package.
    
   The detailed documentation can be found in the corresponding `.Rd` file in [R](R) folder. After installing the package, just run `?somefunction` or `help("somefunction")` in R or RStudio to get the detailed documentation.

- [man](man) folder include `.Rd` files for documenting each function listed in [R](R) folder, which is a standard way of documenting a package supported by R. The `.Rd` files use a custom syntax, loosely based on LaTeX.

- [src](src) folder contains the c++ or fortran code for implementing the AMIAS algorithm, which is developed to speed up the computation.

- Other files in this file folder are some necessary file in developing an R package. For example, the [DESCRIPTION](DESCRIPTION) file provides information on this package, including the version, author, maintainer, license et al..


## References

- Canhong Wen, Xueqin Wang, and Aijun Zhang (2023). $\ell_0$ Trend Filtering. 
