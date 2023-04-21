# AMIAS: an R package for solving the generalized $\ell_0$ problem

`AMIAS` (Alternating Minimization Induced Active Set) is an R package that aims to solve the generalized $\ell_0$ problem. That is, given an observed vector $y$, the generalized $\ell_0$ problem minimizer the sum of squared residuals with $\ell_0$-regularization on some linear combination of the estimator. The AMIAS method is based on the necessary optimality conditions derived from an augmented Lagrangian framework. The proposed method takes full advantage of the primal and dual variables with complementary supports, and decouples the high-dimensional problem into two sub-systems on the active and inactive sets, respectively. A sequential AMIAS algorithm with warm start initialization is developed for efficient determination of the cardinality parameter, along with the output of solution paths.

The most typical example is the $\ell_0$ trend filtering, which is an effective tool for nonparametric regression with the power of automatic knot detection
in function values or derivatives. For more details, please see the paper **$\ell_0$ Trend Filtering** by C. Wen and X. Wang and A. Zhang. 


## Installation

If you are going to install the package in Unix, just skip and go to the next paragraph. In Windows, you will need to install some additional software tools to build and install the `AMIAS` package, which implement the main algorithm of our paper. The Rtools is a collection of build tools, a compiler toolchain, headers and pre-compiled static libraries, and is  used for building R packages from source (those that need compilation of C/C++ or Fortran code). You could download Rtools at [Rtools](http://lib.stat.cmu.edu/R/CRAN/bin/windows/Rtools/). After downloading, double click `Rtools.exe` to install it and its accompanying tools: minGW, perl. Rtools automatically recognizes the paths of those relevant softwares and add them to the environment variables of your computer. For more details on the installation of Rtools, please go to the website [http://lib.stat.cmu.edu/R/CRAN/bin/windows/Rtools/rtools43/rtools.html](http://lib.stat.cmu.edu/R/CRAN/bin/windows/Rtools/rtools43/rtools.html).

Once the repo is locally cloned, the user can follow the following steps to install the `AMIAS` R package:

1. Open a terminal window
2. Go to the directory that contains AMIAS/ directory.
   Type
   ```
   R CMD build AMIAS
   ```
   The the user will see something like this
   ```
   * checking for file ‘AMIAS/DESCRIPTION’ ... OK
   * preparing ‘AMIAS’:
   * checking DESCRIPTION meta-information ... OK
   * cleaning src
   * checking for LF line-endings in source and make files and shell scripts
   * checking for empty or unneeded directories
   * building ‘AMIAS_1.0.3.tar.gz’
   ```
   and a *AMIAS_1.0.3.tar.gz* file will be there.

3. Type
   ```
   R CMD INSTALL AMIAS_1.0.3.tar.gz
   ```
   to install the package.

Alternatively, the user could directly install it from Github without downloading it manually. Start R and type
```r
if(!require(devtools)) install.packages('devtools')
devtools::install_github("INFORMSJoC/2021.0313/scripts/AMIAS")
```
or
```r
if(!require(devtools)) install.packages('devtools')
devtools::install_github("C2S2-HF/L0TF/scripts/AMIAS")
```
to install the package.

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
