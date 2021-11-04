%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% etaA.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function computes the value of eta for a given perturbation Delta0.
% This is the minimum value of the projection of the productivity vector
% over the eigenvector of the effective competition matrix, and it
% determines the feasibility condition.
% 
% Madrid, February 16th, 2016
% Bioinformatics Unit CBMSO (A.P-G-)
%

function eta=etaA(Alpha0P,Alpha0A,GammaEffA,InvBetaP,randvecP,randvecA,Delta0,hP,hA,zP,zA,e1A)

h=hA;
% ... Perturb the alpha values
AlphaP=Alpha0P.*(1+(randvecP-0.5)*Delta0*2);
AlphaA=Alpha0A.*(1+(randvecA-0.5)*Delta0*2);
% ... And repeat the procedure, compute the productivities
if(h ~= 0) % we will work with effective alpha
    AlphaEffP=AlphaP+(1/hP).*((zP./(1+zP)).^2)';
    AlphaEffA=AlphaA+(1/hA).*((zA./(1+zA)).^2)';
    ProdA=AlphaEffA'+GammaEffA*InvBetaP*AlphaEffP';
else
    ProdA=AlphaA'+GammaA*InvBetaP*AlphaP';
end
% ... Normalize productivities
%Prod0A=Prod0A./CdiagA; % the effective competition, comment
%ProdA=ProdA./CdiagA;

% --- Project the final productivities over the eigenvector

Prod1eA=ProdA'*e1A; %ProdeA(idxA);

% --- This is the normalized projection

eta_i=1-ProdA./(e1A.*Prod1eA); %
eta=max(eta_i);

end