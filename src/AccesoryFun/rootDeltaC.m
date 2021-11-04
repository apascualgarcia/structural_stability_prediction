function [DeltaC,m,b]=rootDeltaC(k,RhoEff,etaVec,DeltaVec,etaCritical)
% Find the root taking the two points surrounding DeltaC and fitting
% them to a line
%%%% you may want to USE a solver here? (interp+froot+predict?)
if((k~=1)&&(RhoEff>0)) % The minimum value chosen is lower than etaCritical, this can be modified in the DeltaVec input
    m=(etaVec(k-1)-etaVec(k))/(DeltaVec(k-1)-DeltaVec(k));
    b = etaVec(k-1)-(DeltaVec(k-1)*m);
    DeltaC=(etaCritical-b)/m;
elseif(RhoEff<=0.00001) % there is a singularity here, the system is structurally stable for any perturbation (up to Delta0~1) but this
    %is likely an artifact due to the existence of modules in the competition. This cases are the subject of further research.
    m=0;
    b=0;
    DeltaC=0.99;
else % The system is dynamically unstable, typically this happens if etaCritical<eta(Delta=0) (the feasibility condition does not hold for infinitesimaly small perturbations).
    m=10000; % Should be infinite (a very small perturbation completely leads to extinctions ).
    b=0;
    DeltaC=0;
end
end