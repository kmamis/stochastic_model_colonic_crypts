# stochastic_model_colonic_crypts
Codes for paper "A simple stochastic model for cell population dynamics in colonic crypts" (2023)
by Konstantinos Mamis, Ruibo Zhang, and Ivana Bozic.

MATLAB application distribution_calculator.mlapp takes as imputs the stem cell number,
stem and TA cell division rates, FD cell apoptosis rate, and average total number of cells per crypt,
calulculates TA differentiation rate, and plots the probability distributions of TA and FD and total cell populations.
It also calculates the average TA and FD numbers, the coefficient of variation (CV) for TA, FD and total cell populations,
and the Pearson's correlation coefficient of TA and FD cell populations.

MATLAB .m files contain the codes for the figures in paper denoted by each file's name.
Parallel computing toolbox has been employed to speed up the execution of for-loops.

Data .mat files contain the results from 100,000 direct Gillespie simulations of the branching process. 
Gillespie simulations have been performed by using the code https://github.com/nvictus/Gillespie
TA_Gillespie.mat and FD_Gillespie.mat contain the steady state TA and FD cell populations respectively of each simulation, 
for parameter values typical of healthy human colonic crypts (see paper for details).
TA_Gillespie_varying_gamma.mat and FD_Gillespie_varying_gamma.mat contain the steady state TA and FD cell populations of each simulation, 
for apoptosis rate gamma randomly varying from crypt to crypt: 1/gamma~Uniform[1d, 6d]. 
