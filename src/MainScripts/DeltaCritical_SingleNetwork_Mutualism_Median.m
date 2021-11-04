%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           DeltaCritical_SingleNetwork_Mutualism_Median.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This script computes a prediction of the critical structural perturbation that
% a given system can afford. 
% We start with the alpha obtained when we estimate a fixed point for a given distribution of
% abundances. Then, we compute the productivity vector and the main eigenvector of the
% effective competition matrix, which determines the
% optimal productivity and abundances vectors. Both measures are computed
% considering effective parameters, obtained after linearizing the
% non-linear terms evaluated at the fixed point.
% The next step computes the effective interspecific parameter and the
% propagation of perturbations. 
% Once we get both measures we estimate the critical perturbation.
%%%%%%%
%
% ## Input files: 
% * Mutualistic matrix: we need a file containing the mutualistic matrix to explore. It
%  is considered that the matrix for plants is the transpose of the matrix
%  for animals. It should be located at the directory data/Gamma, see the file
%  Mutualistic_matrix_long.txt for an example. 
%  This should be a matrix in long format in which each line
%  contains the matrix entry i,j for animal i and plant j:
%    matrix(1,1)
%    matrix(1,2)
%    ...
%    matrix(N,M-1)
%    matrix(N,M)
%
% * Abundances vectors (optional): A single column file for plants and another one
%   for animal species with a single header line describing the species
%   abundances at fixed point.
% ## Abundances parameters:
%  * distro: A parameter to determine if the abundances are read from files
%  or generated internally. If:
%       ** 'exp' = read from file
%       ** 'uniform' = are generated from a uniform distribution.
%       ** 'lognorm' = are generated from a lognormal distribution.
%  * assign: If distro = 'lognorm' was used the abundances can be generated
%  in three ways:
%       ** 'random' = there is no relation between the abundances and the
%       species degrees
%      ** 'direct' = more connected species will have higher abundances
%      ** 'inverse' = more connected species will have smaller abundances.
%
%  ## Model parameters
%
%  The remaining parameters corresspond to the dynamical model and are described at the Parameters section
%
%
% ## Output: File with the top props and the minimum critical perturbation predicted as
% a function of the architecture. It is expected an output directory called data/DeltaCritical. Please note that
% the values are returned in a single line. This is because the script
% could be converted in a function to run it several times on a list of
% networks, so each line will contain the properties of each network.
% 

clear all
close all

%%%%%% START EDITING  --------------------------------------
%% Fix paths:

% Loading data
fileGamma='Mutualistic_matrix_long.txt'; % Include here the name of the mutualistic matrix. It should be located at the folder data/Gamma
fileAbundP='AbundancesP.txt'; % Name of the file containing the abundances of plants. It  should be located at the folder data/Abundances
fileAbundA='AbundancesA.txt';  % Name of the file containing the abundances of animals (pollinators, seed dispersals). It  should be located at the folder data/Abundances
SimulDescr='Prediction of DeltaC after computing alpha at fixed point and estimating eta(Delta) with the median of a set of perturbations'; % You can write here a description of your simulation, it will be written in the header of the output file.

%% Parameters
% Unless specified otherwise, parameters that are generated internally follow a
% uniform distibution with mean X0:  Unif(X0*(1-Xwidth),X0*(1+Xwidth)).

% .... Parameters typically explored
rho=0.15 ; % Intraspecific competition parameter and mutualistic strength, typical parameters to be tuned
gamma0=0.15; % mutualistic strength
beta0=1.0;  % Direct competition parameter
N0p=1; % Average biomass parameter for plants
N0a=1; % Average biomass parameter for animals
DeltaRegion=[]; % A specific region to look for the DeltaCritical perturbation. Leave blank for full exploration
Nrnd=100; % Number of perturbations on growth rates
distro='exp'; % 'uniform', 'lognorm' or 'exp'. Distributions of abundances. If "exp" it will be read from file, "uniform" and "lognorm" will generate a random distro with the specified distribution
assign='rand'; % 'rand', 'direct' or 'inverse'. Ways in which you can assign the lognorm abundances to species, see generateLogNormDirectedBipartite.m
MutType='Facultative'; % "Facultative" or "Obligatory". If you fix it as obligatory, N and h will be computed automatically so the input values below will be overwritten
binary=0; % 1 = you input a binary matrix and you want to generate parameters, 0 to read a matrix with pre-computed parameters
Sa=0; % Number of animal species. It must be >0 unless distro = exp
Sp=0; % Number of plant species. It must be > 0 unless distro = exp

%.... Typically fixed parameters
NmidP=1; % Carrying capacity normalization
NmidA=1;
Nwidth=0;
betaWidth=0; % width of the interval around the mean value of beta0  (uniform distribution)
gammaWidth=0; % width around gamma0
rhoWidth=0; % width around rho
h=0.1; % handling times. In the obligatory regime h is function of other parameters, it is fixed later. No noise is considered.
hA=h; hP=h; %

%%%%% STOP EDITING  --------------------------------------

%%  Start computation

fprintf('%s \n','   ');
fprintf('%s \n','***************************************************');
fprintf('%s \n','*  DeltaCritical_SingleNetwork_Mutualism_Median.m');
fprintf('%s \n','***************************************************');
fprintf('%s \n','   ');

%% Open input and output files

% ... fix permanent directories structure
dirScript=mfilename('fullpath'); tmp=strsplit(dirScript,'src');
dirRoot=tmp{1};
dirGamma='data/Gamma/';
dirAbund='data/Abundances/';
dirOut='data/DeltaCritical/';
under='_';

%... Prepare output files
FileOut=['DeltaCmedian-' MutType '_'];
if(strcmp(distro,'lognorm'))
    distroLabel=[distro under assign under 'stdv-' num2str(stdGauss)];
else
    distroLabel=distro;
end
if(strcmp(distro,'exp'))
    Labelout=['_rho' num2str(rho) '-gammaExp' under 'N' under distroLabel  '.txt'];
else
    string=['N' num2str(N0p) under distroLabel];% under 'Corr'];%['h0.1_N' num2str(N0p)]; % Any additional label for the output?
    if(string)
        Labelout=['_rho' num2str(rho) '-gamma' num2str(gamma0) under string '.txt'];
    else
        Labelout=['_rho' num2str(rho) '-gamma' num2str(gamma0) '.txt'];
    end
end
FileOut=[FileOut Labelout]; % Prepare a global output file
cd(dirRoot);
cd(dirOut);
fidOut=fopen(FileOut,'w');
%headout=['SpectrProps-' MutType];

% ... Print Header

fprintf(fidOut,'# Computed with DeltaCritical_allNetworks_Mutualism_Median_Exp.m \n');
fprintf(fidOut,'%s %s \n','# At date and time: ',datestr(now));
fprintf(fidOut,'%s \n','# Author:, apascualgarcia.github.io');
fprintf(fidOut,'%s \n','# Computing the effect of the starting relative angles between the productivities, abundances and their optimal values');
fprintf(fidOut,'%s %s \n','# Simulation specific details: ',SimulDescr);
%fprintf(fidOut,'%s %s \n','# Input File: ',fileGamma);
fprintf(fidOut,'%s %s \n','# Distribution for abundances: ',distroLabel);
fprintf(fidOut,'%s %f %s %f %s %f %s %f %s %f %s %f %s %f %s %f  %s %f %s %f \n','#NmidP= ',NmidP,' NmidA= ',NmidA,' N0p= ',N0p,'N0a= ',N0a,' rho= ',rho,' beta0= ',beta0,...
    ' hA= ',hA,' hP= ',hP,' gamma0=',gamma0,' Delta=[EtaCritical-eps,EtaCritical+eps] with eps=',eps);
fprintf(fidOut,'%s %s \n','#1Network, 2Nest, 3Conn, 4NDNSp, 5DNCos, 6Assort, 7StdvDeg, 8NestOrder2P, 9NestOrder2A, 10RhoEffP, 11RhoEffA, ',...
    '12EtaPrimeP, 13EtaPrimeA, 14InterceptP, 15InterceptA, 16DeltaCp, 17DeltaCa, 18DeltaC, 19NsingletonsA, 20NsingletonsP');

%% START computation
% ... start reading abundances if are provided
cd(dirRoot);
if(strcmp(distro,'exp')) % "experimental" abundances, i.e. from simulations
    % --- different  options to read in the following
    cd(dirAbund); % 
    Na=readAbund(fileAbundA); % Na=Na'; % 
    Sa=length(Na);
    Np=readAbund(fileAbundP); % Np=Np';
    Sp=length(Np);
end
% ... Only if we work with obligatory mutualism, we need to estimate Nmid
if(strcmp(MutType,'Obligatory'))
    [NmidP,NmidA]=fixObligatoryParameters(n,m,beta0,betaWidth,rho);
end
cd(dirRoot);
cd(dirGamma);
data=readFortranGammaOut(fileGamma,Sa,Sp);
adjObs=data';
[n,m]=size(adjObs);
cd(dirRoot);
adj=zeros(n,m);
degree=zeros(n,m);
if(binary == 1) % Generate parameters
    adj(adjObs>0)=1; % Convert in binary if it is not.
    adj=adj.*(gamma0*(1+(rand(n,m)-0.5)*gammaWidth))./sqrt(NmidA*NmidP);
    FileNetOut='Mutualism_matrix_Out.txt'; 
    cd(dirGamma); % If we generate parameters we will print the new matrix
    dlmwrite(FileNetOut,adj);
else % keep the matrix obtained.
    adj=adjObs;
end

%-- Compute degrees, identify singletons
degree(adj>0)=1;
degreeA=sum(degree,1);
degreeP=sum(degree,2);
IDsingleA=degreeA;
IDsingleA(degreeA>1)=0; % Cancel the index for those species with more than one connection
singletonsA=sum(IDsingleA); % Count how many species have one connection
IDsingleP=degreeP;
IDsingleP(degreeP>1)=0;
singletonsP=sum(IDsingleP);

%--- Biomasses (just if distro != exp )
if(strcmp(distro,'uniform'))
    [Np,Na]=generateAbundancesBipartite(distro,n,m,NmidP,NmidA,Nwidth);
elseif(strcmp(distro,'lognorm'))
    [Np,Na]=generateLogNormDirectedBipartite(n,m,Nwidth,adjObs,assign);
end

%--- Gamma (effective)
GammaP=adj;
GammaA=adj';
if(h ~=0) % We work with effective gamma matrices
    [GammaEffP,GammaEffA] = EffGammaMatrix(GammaP,hP,hA,Np',Na');
else
    GammaEffP=GammaP;
    GammaEffA=GammaA;
end

%--- Effective/Competition matrices
[BetaP,BetaA,InvBetaP,InvBetaA,Cp,Ca,Bp,Ba]=CompMatrices(n,m,beta0,betaWidth,rho,rhoWidth,NmidP,NmidA,GammaEffP,GammaEffA);

% --- DeltaCritical. The growth rates and effective productivities are
% computed within this function
[DeltaC,DeltaCp,DeltaCa,mP,mA,bP,bA,RhoEffP,RhoEffA]=deltaCritical_Mutualism_Median(DeltaRegion,BetaP,BetaA,InvBetaP,InvBetaA,Cp,Ca,GammaP,GammaA,hP,hA,Np,Na,Nrnd);

%-- Compute topological properties
Nest = Nestedness_Overlap(adjObs);
NDNSp=DoubleNestedness_Spectral(adjObs);
DNCos=DoubleNestedness_Cosine(adjObs);
DNOv=DoubleNestedness_Overlap(adjObs);
Assort=Assortativity_Und(adjObs);
Conn=Connectance(adjObs);
StdvDeg=Stdv_Degree(adjObs);
[NO2P,NO2A]=Nestedness_Order2(adjObs,rho);

fprintf(fidOut,'%s %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f \n',...
    fileGamma,Nest,Conn,NDNSp,DNCos,Assort,StdvDeg,NO2P,NO2A,...
    RhoEffP,RhoEffA,mP,mA,bP,bA,DeltaCp,DeltaCa,DeltaC,singletonsA,singletonsP);

fprintf('%s \n','*                                            * ');
fprintf('%s \n','* ... Done! Check the results...             *');
fprintf('%s \n','**********************************************');
fprintf('%s \n','   ');
