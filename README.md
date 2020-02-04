# BESIDE-MR
BayEsian Set IDEntification Mendelian randomization

## Installation
To install `BESIDEMR` directly from the GitHub repository, first make sure you have the `devtools` package installed:

    install.packages("devtools")

Then the `BESIDEMR` package can be installed using:

    library(devtools)
    install_github("CYShapland/BESIDEMR")
    
To update the package just run the `install_github("CYShapland/BESIDEMR")` command again.

## Description

The main function is `BESIDE_MR` to perform BayEsian Set IDEntification Mendelian randomization (BESIDE-MR). We develop a bespoke Metropolis-Hasting algorithm to perform the search using the recently developed Robust Adjusted Profile Likelihood (MR-RAPS) of Zhao et al as the basis for defining a posterior distribution that efficiently accounts for pleiotropic and weak instrument bias. BESIDE-MR can be extended from a standard one-parameter causal model to a two-parameter model, to allow a large proportion of SNPs to violate the Instrument Strength Independent of Direct Effect (InSIDE) assumption.

`BESIDE_MR` returns an object of class `beside`, consists of the posterior of effect estimate, pleiotropy variance and instrument inclusion indicator variable from each iteration. The estimation of variance is challenging, hence we have included `tau_estimate` with the options of `DL estimate` and `Full_Bayes`, where the former is a plug-in estimate for pleiotropy variance.

## Citation

The corresponding paper can be accessed at:

Profile-likelihood Bayesian model averaging for two-sample summary data Mendelian randomization in the presence of horizontal pleiotropy.

## License

This project is licensed under GNU GPL v3.
