%********************************************%
%            Nestedness_Overlap.m            %
%********************************************%

% function Nestedness = Nestedness_Overlap(input);

 %function [NestednessP,NestednessA] = Nestedness_Overlap(input);
function [NestednessP,NestednessA,Nestedness] = Nestedness_Overlap(input);
% Computes a nestedness measure as defined in Bastolla et al. (2009)
% Nature


% Transform an arbitrary input into a binary matrix:
[m,n]=size(input);
adj=zeros(m,n);
input=input./(max(max(input))); % Normalize values between zero and one
adj(input>0)=1; % Turn to one those cells of adj where you find in input positive values

Limit=0;
Nestedness=0;
NshareTotA=0; % Connectivity between elements from the same community.

for i=1:m-1;
    for j=i+1:m;
        Nshare=sum(adj(i,:).*adj(j,:)); 
        NshareTotA=NshareTotA+Nshare;
        % Nshare<=min(Ni,Nj) for each i,j by nestedness definition.
        Ni=sum(adj(i,:));
        Nj=sum(adj(j,:));
        Limit=Limit+min(Ni,Nj);
        Nestedness=Nestedness+Nshare;                                 
    end
end
LimitT=Limit;
NestednessT=Nestedness;
NestednessA=Nestedness/Limit;

Limit=0;
Nestedness=0;
NshareTotP=0;

for i=1:n-1;
    for j=i+1:n;
        Nshare=sum(adj(:,i).*adj(:,j)); 
        NshareTotP=NshareTotP+Nshare;
        % Nshare<=min(Ni,Nj) for each i,j by nestedness definition.
        Ni=sum(adj(:,i));
        Nj=sum(adj(:,j));
        Limit=Limit+min(Ni,Nj);
        Nestedness=Nestedness+Nshare;                                 
    end
end
LimitT=LimitT+Limit;
NestednessT=NestednessT+Nestedness;
NestednessP=Nestedness/Limit;
%Nestedness=NestednessT/LimitT;
Nestedness=(NestednessP+NestednessA)/2;

end
