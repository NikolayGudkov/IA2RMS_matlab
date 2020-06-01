# Overview

These MATLAB codes implement the Independent Doubly Adaptive Rejection Metropolis Sampling (IA2RMS) withing Gibbs Sampling algorithm.

The original algorithm is described in [L. Martino, J. Read, D. Luengo, "Independent Doubly Adaptive Rejection Metropolis Sampling within Gibbs Sampling",
IEEE Transactions on Signal Processing, Volume 63, Issue 12, Pages: 3123-3138, 2015.](http://www.lucamartino.altervista.org/TSP_IA2RMS.pdf)

This version of the codes (v0.5) aims to improve the computational performance of the version (v0.4) from (https://sourceforge.net/projects/a2rms/files/IA2RMS/)

The main difference between (v0.5) and (v0.4) is threefold:
* a new function "update_proposal.m' for updating the proposal distribution is added. This function "inserts" a point into an interval from which it has been simulated and updated the proposal accordingly.
* the approach to select an interval to simulate from is made more efficient;
* the new type of proposal construction (exponential pieces/linear in log-scale) is added;

We implement three types of proposal construction with:
*  uniform pieces (type '0');
*  linear pieces (type '1').
*  exponential pieces (type '2');

The code was tested on two examples from [L. Martino, J. Read, D. Luengo (2015)](http://www.lucamartino.altervista.org/TSP_IA2RMS.pdf):
1. Multimodal Target: Mixture of Gaussians,
1. Heavy-Tailed Distribution: Levy distribution.

# Content

* __build_proposal.m__ - builds a proposal from a set of support points
* __update_proposal.m__ - updates the proposal distribution by inserting a point in an interval
* __checkInitPoints.m__ - checks the initial points to ensure that there are no potential
* __Eval_Proposal.m__ - evaluates proposal distribution at a given point
* __Sample_Eval_Proposal.m__ - samples from a proposal distribution and evaluates proposal at the simulated point
* __simulate_index.m__ - simulates index of the interval from which we want to sample the proposal
* __IA2RMS.m__ - the realisation of the IA2RMS algorithm

# Usage 

## Description
Function **IA2RMS** simulates *N* points from the target continuous distribution provided by *target_function*.
### Input:
* __*f*__ - continuous probability distribution function known up to a constant of proportionality.
* __*S*__ - initial set of support points. 
* __*M*__ - number of random samples required.
* __*type*__ - the type of construction of proposal.

### Output: 
* vector of simulated points.
