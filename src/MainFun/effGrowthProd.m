function [AlphaEffP,AlphaEffA,ProdP,ProdA]=effGrowthProd(AlphaP,AlphaA,GammaP,GammaA,InvBetaP,InvBetaA,hP,hA,Np,Na)
% This function computes the effective growth rates, which are the
% growth rates obtained after linearizing the bare growth rates evaluated
% at the stationary abundances, and the effective productivity vectors.
% INPUT
% - Growth rates vectors Alpha for both pools of plants and animals
% - Mutualistic matrices
% - Vector abundances
% - 
h=hP; % Arbitrary choice, just to check if there are saturating terms
if(h ~= 0) % we will work with effective alpha
    zP=hP.*GammaP*Na'; %Saturating terms
    zA=hA.*GammaA*Np';
    AlphaEffP=AlphaP+(1/hP).*((zP./(1+zP)).^2)';
    AlphaEffA=AlphaA+(1/hA).*((zA./(1+zA)).^2)';
    [GammaEffP,GammaEffA] = EffGammaMatrix(GammaP,hP,hA,Np',Na');
else
    AlphaEffP=AlphaP;
    AlphaEffA=AlphaA;
    GammaEffP=GammaP;
    GammaEffA=GammaA;
end
ProdP=AlphaEffP'+GammaEffP*InvBetaA*AlphaEffA';
ProdA=AlphaEffA'+GammaEffA*InvBetaP*AlphaEffP';

end