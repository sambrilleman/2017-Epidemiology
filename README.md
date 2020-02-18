# Supporting materials for published article

This repository contains the supplementary materials, including computing R code, for fitting the models described in the following paper:

Brilleman SL, Howe L, Wolfe R, Tilling K. [Bayesian piecewise linear mixed models with a random change point: an application to BMI rebound in childhood](https://journals.lww.com/epidem/fulltext/2017/11000/Bayesian_Piecewise_Linear_Mixed_Models_With_a.10.aspx). *Epidemiology* 2017;**(6)**:827-833. doi: [10.1097/EDE.0000000000000723](https://doi.org/10.1097/EDE.0000000000000723)

These materials however are subject to change, to make sure they stay up-to-date and working with their dependencies, and hopefully also with some improvements to the code over time. The ultimate aim will be to develop an R package that will provide a general purpose user-friendly interface for fitting these random change point models. 

## Getting Started

Estimation of the random change point models is done using the Bayesian software [Stan](http://mc-stan.org). In order to use Stan from R, you will need to install RStan, which is the R iterface to Stan. The details for installing RStan can be found [here](https://github.com/stan-dev/rstan/wiki/RStan-Getting-Started).

## Bug Reports

If you find any bugs, please report them via email to [Sam Brilleman](mailto:sam.brilleman@monash.edu).

## References

1. Brilleman SL, Howe L, Wolfe R, Tilling K. [Bayesian piecewise linear mixed models with a random change point: an application to BMI rebound in childhood](https://journals.lww.com/epidem/fulltext/2017/11000/Bayesian_Piecewise_Linear_Mixed_Models_With_a.10.aspx). *Epidemiology* 2017;**(6)**:827-833. doi: [10.1097/EDE.0000000000000723](https://doi.org/10.1097/EDE.0000000000000723)

2. Stan Development Team (2015) Stan Modeling Language Users Guide and Reference Manual. [http://mc-stan.org/documentation/](http://mc-stan.org/documentation/)
