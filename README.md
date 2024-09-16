# Generalized Hyperbolic Models for Exchange Rate Log-Returns

## Description

This project explores the generalized hyperbolic fitting functions provided by the "ghyp" package in R. Specifically, the log-returns of three different exchange rates are fitted to the following distributions:
- Generalized Hyperbolic
- Normal Inverse Gaussian
- skewed Student-t
- Variance-Gamma distributions.
  
Having obtained parameters from the fitting functions for each distribution, the EM algorithm is applied by having the Generalized Inverse Gaussian act as a mixing distribution to generate and analyze the volatility of the log-returns.

A table is also generated, presenting the log-likelihood estimates of each model.

## Installations

Prior to running this code, the "ghyp" is required to be installed in order to access the fitting functions of the different generalized hyperbolic distributions.
The "quantmod" function is also necessary to use function "getSymbols" and retrieve historical data of the exchange currencies.

## Usage

The libraries from both installed packages must be loaded to run this code.
While the code presented performs the experiment using exchange rate from American dollar to Euro, Ruble, and Chinese Yuan, other exchange currencies may be used when retrieving data using "getSymbols".
Note that the code provided seeks to observe the changes in volatility between a five year interval from 2017 to 2022. In order to alter this interval, three parts must be edited:
- the start and end dates of the new interval for the "getSymbols" functions
- the start and end dates of the Training and Test for the log-likelihood table
- the loop function used to generate the results from the EM algorithm (make the function loop a number of times equal to the number of log-return entries)
