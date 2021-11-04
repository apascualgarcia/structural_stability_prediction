# Structural stability in model ecosystems 

This repository contains code to compute the critical perturbation on the growth rates that a Lotka-Volterra mutualistic system can afford, as described in:

> _Pascual-García, A. & Bastolla, U._ (2017). Mutualism supports biodiversity when the direct competition is weak. Nature Communications, 8(1), 1-13. 

Please note that this is not the code used in the article but an adaptation to illustrate the method, and it has not been exhaustively tested. 


## Preliminaries

The idea behind our approach is that we observe a specific system at a steady state and we aim to quantify which is the maximum amplitude of an environmental perturbation affecting the species growth rates that the system can afford until a species get extinct. This implies the following considerations:

* This approach measures the structural stability of _a system_, as opposed to the structural stability of _an interaction matrix_, which is the approximation followed by other groups when they compute the feasibility volume. See [this commentary](https://www.sjscience.org/article?id=585) for further clarifications. This means that the specific abundances at equilibrium are required as an input.

* This code considers a system hosting two pools of species (e.g. plants and their pollinator animals) with intraspecific and interspecific competition within each pool and mutualistic interactions between pools. A similar computation can be derived for prey-predator networks or purely competitive systems.

* The critical perturbation is a semi-empirical prediction meaning that, for each metaparameter specified in the model, we generate a large number of random perturbations and we estimate the minimal amplitude of the perturbation having a probability equal to 0.5 of observing at least one extinction (see MS for details).

* The model considers a Holling type II saturating response for the mutualistic effects and, then, we linearize the system to perform our predictions.

* Given the model parameters and species abundances the code specifies the growth rates compatible with the existence of a fixed point. 

## Code description

The code was developed in Matlab and hence it does not require installation, you just need to be sure that all the folders under `src` are in Matlab's path. _I ask for forgiveness to Stallman's spirit for developing scientifc code in a language requiring a non-free license, this is just a relique of my branches as a student in Physics and I promise that all my new code is developed in open and free languages :-)._

* The main script is `src/MainScripts/DeltaCritical_SingleNetwork_Mutualism_Median.m`. You will have to specify a number of parameters in that file and, in addition, you should provide a  mutualistic matrix (binary or not). You can additionally provide some input abundances, which will be internally generated otherwise following one of the different prossibilities available.


#### Input parameters

You should fix the following parameters in the main code. Default values are indicated. In most cases, these are metaparameters (let's say X0) indicating the mean of a uniform distribution:  Unif(X0*(1-Xwidth),X0*(1+Xwidth)). 

** Parameters with noise **

* Metaparameter mean
   * rho=0.15 ; Interspecific competition parameter between plant species or animal species
   * gamma0=0.15;  Mutualistic strength
   * beta0=1.0;   Intraspecific competition parameter
   * N0p=1;  Average biomass parameter for plants
   * N0a=1;  Average biomass parameter for animals

* Width of the distribution, as indicated above
   * Nwidth=0.05;
   * betaWidth=0.05; % width of the interval around the mean value of beta0  (uniform distribution)
   * gammaWidth=0.05; % width around gamma0
   * rhoWidth=0.05; % width around rho

** Parameters with no noise** 

* Model-related
   * NmidP=1;  Carrying capacity normalization
   * NmidA=1;
   * Sa=0;  Number of animal species. It must be >0 unless distro = exp
   * Sp=0;  Number of plant species. It must be > 0 unless distro = exp
   * h=0.1;  handling times. In the obligatory regime h is function of other parameters so it will estimated internally.   

* Computation-related
   * DeltaRegion=[]; A specific region to look for the DeltaCritical perturbation. Leave blank for a full exploration
   * Nrnd=100;  Number of perturbations on the growth rates

* Specific modelling assumptions and handy options

   * distro='exp';  This parameter specifies the abundances, either you can read the abundances from a file  (`exp`, i.e. "experimental") or you can generate them randomly following a `uniform` or `lognorm` distribution.
   * assign='rand'; %  If a lognormal distribution is generated since it is very skewed it may have a deep impact in the structural stability. To test this extent you can assign the abundances generated to the different species either randomly (`rand`), assigning the species with the highest mutualistic degree the highest abundances (`direct`) or the species with the lowest degrees the highest abundances (`inverse`).
   * MutType='Facultative'; You can test `Facultative` or `Obligatory` mutualism with this choice. Since feasibility conditions for obligatory mutualism are challenging, you will not be allowed to read from file the abundances if `MutType=Obligatory`. Instead, the abundances and h will be computed automatically asuming N0A=1 following a uniform distribution and estimating N0P accordingly.
   * binary=0; If your input matrix is binary (or if you want to override the parameters and create new ones) set `binary=1`.  `binary=0` will instead read your matrix with pre-computed parameters.

#### Input files

The code requires a mutualistic matrix (obligatory) and optionally abundances vectors. You have to fix the name of your files in the main script. If you generate random abundances internally, the name of the file included will be ignored.

* Mutualistic matrix: should be located in `data/Gamma` and has the format (see the file `Mutualistic_matrix_long.txt` for an example):

```
This should be a matrix in long format in which each line contains the matrix entry i,j for animal i and plant j:
    matrix(1,1)
    matrix(1,2)
    ...
    matrix(N,M-1)
    matrix(N,M)
```

* Abundances vectors: A single column file for plants and another one for animal species with a single header line (see files `Abundances.txt` for examples).


#### Output file

The code will return an output file with a single line containing several quantities. This format was chosen because the user may want to transform the script in a function and loop it through a large set  of networks, so she will retrieve one line for each network. Output metrics:

* 1Network: Name of the input file of the network
* 2Nest: Nestedness as defined in Bastolla et al. 2009 and subsequent work.
* 3Conn: Connectance of the network.
* 4NDNSp and 5DNCos: Two measures of nestedness that consider the parameters.
* 6Assort: Assortativity of the network.
* 7StdvDeg: Standard deviation of degree.
* 8NestOrder2P and 9NestOrder2A: A corrected version of the nestedness metric.
* 10RhoEffP and 11RhoEffA: Effective competition
* 12EtaPrimeP and 13EtaPrimeA: Propagation of perturbations
* 14InterceptP and 15InterceptA: Vulnerability
* 16DeltaCp, 17DeltaCa and 18DeltaC: Critical perturbations of plants, animals and total.
* 19NsingletonsA and 20NsingletonsP: Number of species with a single mutualistic connection (when is high sometimes explains bad predictions).






