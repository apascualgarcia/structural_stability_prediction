function [DeltaC,DeltaCp,DeltaCa,mP,mA,bP,bA,RhoEffP,RhoEffA]=deltaCritical_Mutualism_Median(DeltaRegion,BetaP,BetaA,InvBetaP,InvBetaA,Cp,Ca,GammaP,GammaA,hP,hA,Np,Na,Nrnd)

% Check user supplied optional parameters
if(isempty(DeltaRegion))
    DeltaVec=[0:0.01:0.99]; % Sample the whole interval
else
    DeltaVec=DeltaRegion; % Or a specific one
end
Ndelta=size(DeltaVec,2);
if(isempty(Nrnd)) Nrnd=100; end

n=size(Cp,1);
m=size(Ca,1);

[Ep,Dp]=eig(Cp);  [Ea,Da]=eig(Ca);
[rp,e1P]=eigenMax(Dp,Ep);   [ra,e1A]=eigenMax(Da,Ea);

traceCp=trace(Cp);
traceCa=trace(Ca);
RhoEffP=(n*rp-traceCp)/(traceCp*(n-1)); % recover the effective competition, formula if C is used
RhoEffA=(m*ra-traceCa)/(traceCa*(m-1));
S0p=(1-RhoEffP)/RhoEffP; % Compute S0, needed for DeltaCritical
S0a=(1-RhoEffA)/RhoEffA;
etaCriticalP=S0p/(S0p+n);
etaCriticalA=S0a/(S0a+m);

% --Find the feasible alpha for the model given
[Alpha0P,Alpha0A,MeanAlphaP,MeanAlphaA,VarAlphaP,VarAlphaA]=alphaFixedPoint(BetaP,BetaA,GammaP,GammaA,hP,hA,Np,Na);
Alpha0P=Alpha0P';
Alpha0A=Alpha0A';

fprintf('%s \n','   ');
fprintf('%s %f \n','    ~~~~ Looking for Delta Critical... ');
fprintf('%s \n','   ');
etaVecP=zeros(1,Ndelta);
etaVecA=zeros(1,Ndelta);
keyP=0;
keyA=0;
for k=1:Ndelta % Delta = 1 should be the maximum value unless weird things happen, see below
    Delta0=DeltaVec(k);
    etaRndP=zeros(1,Nrnd);
    etaRndA=zeros(1,Nrnd);
    for j=1:Nrnd % For each value of Delta we perform a large number of perturbations         
        %  --- Compute the vulnerability
        randvecP=rand(1,n);
        randvecA=rand(1,m);
        AlphaInP=Alpha0P.*(1+(randvecP-0.5)*Delta0*2);
        AlphaInA=Alpha0A.*(1+(randvecA-0.5)*Delta0*2);
        [AlphaEffP,AlphaEffA,ProdP,ProdA]=effGrowthProd(AlphaInP,AlphaInA,GammaP,GammaA,InvBetaP,InvBetaA,hP,hA,Np,Na);        
        eta1P=eta(ProdP,e1P);
        etaRndP(j)=eta1P;                
        eta1A=eta(ProdA,e1A);
        etaRndA(j)=eta1A;
    end 
    etaVecP(k)=1-median(etaRndP);
    etaVecA(k)=1-median(etaRndA);
    if((keyP==0)&&(etaVecP(k) >= etaCriticalP))
        keyP=1;
        idxP=k;        
    end
    if((keyA==0)&&(etaVecA(k) >= etaCriticalA))
        keyA=1;
        idxA=k;        
    end
    if((keyP==1)&&(keyA==1))
        break;
    end
end

if((keyP==0)&&(keyA==0)) % This should not happen, warn about it
    DeltaCp=1; DeltaCa=1; DeltaC=1; mP=1000; bP=0; mA=1000; bA=0;
    fprintf('%s \n','*** WARNING');
    fprintf('%s %f \n','    ~~~ Delta Critical not found, fixing it equal to one!');
elseif(keyP==0)
    DeltaCp=1; mP=1000; bP=0;
    [DeltaCa,mA,bA]=rootDeltaC(idxA,RhoEffA,etaVecA,DeltaVec,etaCriticalA);
    fprintf('%s \n','*** WARNING');
    fprintf('%s %f \n','    ~~~ Delta Critical of plants not found, fixing it equal to one!');
elseif(keyA==0)
    DeltaCa=1; mA=1000; bA=0;
    [DeltaCp,mP,bP]=rootDeltaC(idxP,RhoEffP,etaVecP,DeltaVec,etaCriticalP);
    fprintf('%s \n','*** WARNING');
    fprintf('%s %f \n','    ~~~ Delta Critical of animals not found, fixing it equal to one!');
else
    [DeltaCa,mA,bA]=rootDeltaC(idxA,RhoEffA,etaVecA,DeltaVec,etaCriticalA); 
    [DeltaCp,mP,bP]=rootDeltaC(idxP,RhoEffP,etaVecP,DeltaVec,etaCriticalP);
end

if(DeltaCp<DeltaCa)
    DeltaC=DeltaCp;
    mT=mP;
    bT=bP;
else
    DeltaC=DeltaCa;
    mT=mA;
    bT=bA;
end

fprintf('%s \n','   ');
fprintf('%s %f \n','    ~~~ Delta Critical found: ',DeltaC);
fprintf('%s \n','   ');
end