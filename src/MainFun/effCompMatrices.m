
function [Cp,Ca,Bp,Ba]=effCompMatrices(BetaP,BetaA,InvBetaP,InvBetaA,GammaP,GammaA,n,m)


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