
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CompMatrices.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% This function generates the pure competition matrix and its inverse and the
% effective and normalized effective competition matrices for plants and
% animals for a given set of parameters, including its randomization.
%
% Madrid



function [BetaP,BetaA,InvBetaP,InvBetaA,Cp,Ca,Bp,Ba]=CompMatrices(n,m,beta0,betaWidth,rho,rhoWidth,Np,Na,GammaP,GammaA)


% -- Pure and inverse competition matrix for plants
BetaP=rho.*(1+(rand(n,n)-0.5)*rhoWidth); 
aux=ones(n,n)-eye(n,n);
BetaP=BetaP.*aux;
aux=beta0.*(1+(rand(n,n)-0.5)*betaWidth);
aux=eye(n,n).*aux;
BetaP=BetaP+aux; % Pure competition matrix
BetaP=BetaP./Np;
InvBetaP=BetaP^(-1); % Inverse pure competition matrix
% -- Pure and inverse competition matrix (Animals):
BetaA=rho.*(1+(rand(m,m)-0.5)*rhoWidth);
aux=ones(m,m)-eye(m,m);
BetaA=BetaA.*aux;
aux=beta0.*(1+(rand(m,m)-0.5)*betaWidth);
aux=eye(m,m).*aux;
BetaA=BetaA+aux; % Pure competition matrix
BetaA=BetaA./Na;
InvBetaA=BetaA^(-1); % Inverse pure competition matrix

% -- Effective competition matrices (and normalized effective competition)
Cp=BetaP-GammaP*InvBetaA*GammaA; % Effective competition matrix
for vv=1:n % Animals
    for ww=1:n
        Bp(vv,ww)=Cp(vv,ww)/(sqrt(Cp(vv,vv)*Cp(ww,ww))); % Normalized effective competition matrix
    end
end
Ca=BetaA-GammaA*InvBetaP*GammaP; % Effective competition matrix
for vv=1:m % Animals
    for ww=1:m
        Ba(vv,ww)=Ca(vv,ww)/(sqrt(Ca(vv,vv)*Ca(ww,ww))); % Normalized effective competition matrix
    end
end

end
