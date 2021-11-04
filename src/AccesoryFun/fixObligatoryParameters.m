

function [NmidP,NmidA]=fixObligatoryParameters(n,m,beta0,betaWidth,rho)
% This function estimates the average values of the abundances for
% obligatory mutualistic systems to ensure that the growth rates lies in
% this regime. Some parameters can be modified, because the parameters
% chosen are rather conservative. See Pascual-Garcia & Bastolla, Nature
% Comm (2017) for further details.
% A. Pascual-Garcia
% Bioinformatics Unit (CBMSO)
gammaMin=0.01; % Conservative minimum mutualistic strength expected in the simulations
bup=beta0+beta0*(betaWidth/2); % Conservative competition values to ensure feasibility
blow=beta0-beta0*(betaWidth/2);
hmaxA=1/(bup*((m-1)*rho+1));
hmaxP=1/(bup*((n-1)*rho+1));
hA=0.75*hmaxA; % Based in exploratory simulations and the computations above
hP=0.25; %hmaxP;
Nrate=(bup*(1+(m-1)*rho))^2/(gammaMin*blow*(1-hA*(bup*(1+(m-1)*rho))))^2;
N0a=1;
N0p=Nrate;
NmidA=N0a;
NmidP=N0p;
end