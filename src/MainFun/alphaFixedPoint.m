


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% alphaFixedPoint.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function computes the growth rates, its mean and variance
% given the pure competition matrix, the adjacency term already with the
% metaparameters included (not the binary matrix), handling times and abundances.
%
% Madrid, March 25th, 2014
% Bioinformatics Unit (CBMSO)
% A.P-G.
%


function [AlphaP,AlphaA,MeanAlphaP,MeanAlphaA,VarAlphaP,VarAlphaA]=alphaFixedPoint(BetaP,BetaA,gamma0P,gamma0A,hP,hA,Np,Na)

HollingP=1+hP.*gamma0P*Na';
HollingA=1+hA.*gamma0A*Np';

GammaP=gamma0P*Na'./HollingP;
GammaA=gamma0A*Np'./HollingA;

AlphaP=BetaP*Np'-GammaP;
AlphaA=BetaA*Na'-GammaA;

MeanAlphaP=mean(AlphaP);
VarAlphaP=var(AlphaP);
MeanAlphaA=mean(AlphaA);
VarAlphaA=var(AlphaA);

end
