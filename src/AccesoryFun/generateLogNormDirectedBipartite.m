function [Np,Na]=generateLogNormDirectedBipartite(n,m,Nwidth,adjObs,assign)
% This function generates abundances from a lognormal distribution, and
% then assign the values to the species following three different
% possibilities. Either higher abundances are assigned to species with high
% connectances (direct), to species with low connectances (inverse), or
% randomly (rand).
% INPUT
% - number of species of plants (n) and animals (m)
% - Nwidth, the std of the gaussian that will generate the lognormal
% distribution (see details in lognormalGenerator.m
% - adjObs, the adjacency matrix linking the two pools of a species in a
% bipartite network (mutualistic or prey-predator)
% - assing, the way in which abundances are assigned (rand, direct, inverse).
% OUTPUT
% - Two vectors of abundances for each pool of species.
%%%%%%%%%%%%
% Bioinformatics Unit (CBMSO)
% A. Pascual-Garcia
%%%%%%%%%%%%

NtmpP=lognormalGenerator(n,Nwidth);
NtmpA=lognormalGenerator(m,Nwidth);
adjDegree=zeros(n,m);
adjDegree(adjObs>0)=1; % Control that you work with a binary matrix
degP=sum(adjDegree,2);
degA=sum(adjDegree,1);
[degSortP,idDegP]=sort(degP);
[degSortA,idDegA]=sort(degA);
[NsortP,idNp]=sort(NtmpP);
[NsortA,idNa]=sort(NtmpA);
Np=zeros(1,n);
Na=zeros(1,m);
%fprintf('%s \n','uu, Np(uu), degP(uu)');
for uu=1:n
    if(strcmp(assign,'rand'))
        Np(uu)=NtmpP(uu);
    elseif(strcmp(assign,'direct'))
        Np(idDegP(uu))=NtmpP(idNp(uu));
    else
        Np(idDegP(n-uu+1))=NtmpP(idNp(uu));
    end
    %fprintf('%i %f %i \n',uu, Np(uu), degP(uu));
end
%fprintf('%s \n','uu, Na(uu), degA(uu)');
%plot(degP,Np,'o');
for uu=1:m
    if(strcmp(assign,'rand'))
        Na(uu)=NtmpA(uu);
    elseif(strcmp(assign,'direct'))
        Na(idDegA(uu))=NtmpA(idNa(uu));
    else
        Na(idDegA(m-uu+1))=NtmpA(idNa(uu));
    end
    %fprintf('%i %f %i \n',uu, Na(uu), degA(uu));
end


end