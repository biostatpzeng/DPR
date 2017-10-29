
DPR: Non-parametric genetic prediction of complex traits with latent Dirichlet process regression models
========================================================================================================
# Introduction

**DPR** is the software written in [C++](http://www.cplusplus.com/) for implementing the latent Dirichlet Process Regression (DPR) model for robust genetic prediction of complex traits with typed genotypes. While most existing GWAS prediction models use parametric priors, DPR makes use of a flexible non-parametric prior on the SNPs effect sizes. By using a non-parametric model, DPR is adaptive to a broad spectrum of genetic architectures and can achieve robust predictive performance for a variety of complex traits. DPR can be fitted using two complementary algorithms: the Monte Carlo Markov Chain (MCMC) algorithm, and the mean filed variational Bayesian (VB) approximation algorithm.

Currently, DPR is only supported for Linux platform.

# Installing and Compiling DPR

If you have downloaded a binary executable [procedure](https://github.com/biostatpzeng/DPR/blob/master/DPR), no installation is necessary. In some cases, you may need to use “chmod a+x DPR” before using the binary executable. In addition, notice that the end-of-line coding in Windows (DOS) is different from that in Linux, and so you may have to convert input files using the utility "dos2unix" or "unix2dos" in order to use them in a different platform.

If you want to compile DPR by yourself, you will need to download the source code, and you will need a standard C/C++ compiler such as GNU gcc, as well as the GSL and LAPACK libraries. You will need to change the library paths in the Makefile accordingly.
A sample Makefile is provided along with the source code. For details on installing GSL library, please refer to [GSL](http://www.gnu.org/s/gsl/). For details on installing LAPACK library, please refer to [LAPACK](http://www.netlib.org/lapack/). The DPR software also relies on the [Eigen library](http://eigen.tuxfamily.org/index.php?title=Main_Page), which is included in the source code folder.

## Cite
#### [Ping Zeng](https://github.com/biostatpzeng) and Xiang Zhou (2017). Non-parametric genetic prediction of complex traits with latent Dirichlet process regression models. Nature Communications, 8(1): 456, Doi: 10.1038/s41467-017-00470-2.

# Contact
We are very grateful to any questions, comments, or bugs reports; and please contact [Ping Zeng](https://github.com/biostatpzeng) via zpstat@xzhmu.edu.cn or ~~pingzeng@umich.edu~~.



