all: 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18

%.html: %.Rmd
	Rscript -e "rmarkdown::render(\"$<\")"

0: basic_mcmc/nimble_basic_mcmc.html
1: build_a_model/nimble_build_a_model.html
2: CAR/CAR.html
3: customized_mcmc/nimble_customizing_mcmc.html
4: Ecology_Examples/Ecology_Examples.html
5: gaussian_process/gaussian_process.html
6: irt_models_example/IRT_example.html
7: linear_predictors/linpred.html
8: logistic_regression/nimble_logistic_regression.html
9: MLE/MLE.html
10: Partial_Pooling/PartialPooling.html
11: pumpMCEM/pumpMCEM.html
12: RJMCMC_example/RJMCMC_example.html
13: simulation_from_model/simulation_from_model.html
14: stochastic_volatility/stochastic_volatility.html
15: zero_inflated_poisson/zero_inflated_poisson.html
16: parallelization/parallelizing_NIMBLE.html
17: bnp/intro_bnp_density.html
18: bnp/intro_bnp_raneff.html
