%********************************%
%   Nestedness_Order2.m          %
%********************************%
% Computes the Order2 Nestednes of a bipartite adjacency undirected matrix
% i.e. the nestedness corrected by higher order terms found in the
% expansion of the effective competition for weak mutualism.
% Input: adjacency matrix, rho Output: Feedback.
% Madrid, February 14th, 2014
% Bioinformatics Unit CBMSO, A.P-G.
%

function [NestO2P,NestO2A]=Nestedness_Order2(input,rho)
Shat=(1-rho)/rho; % It is the same value for plants and animals just if rho is the same for both
[Sp,Sa]=size(input);

SSa=1/(Shat+Sa);
SSp=1/(Shat+Sp);

adj=zeros(Sp,Sa);

input=input./(max(max(input))); % Normalize values between zero and one
adj(input>0)=1; % Turn to one those cells of adj where you find in input positive values

adjP=adj;
adjA=adj';

degP=sum(adj,2);
degA=sum(adj,1);

for i=1:Sp
    for j=1:Sp
        Mp(i,j)=0;
        for k=1:Sa
            Mp(i,j)=Mp(i,j)+((adjP(i,k)-degA(k)*SSp)*(adjA(k,j)-degP(j)*SSp));                               
        end
    end
end
for i=1:Sa
    for j=1:Sa
        Ma(i,j)=0;
        for k=1:Sp
            Ma(i,j)=Ma(i,j)+((adjA(i,k)-degP(k)*SSa)*(adjP(k,j)-degA(j)*SSp));                               
        end
    end
end
Mp=(Mp+Mp')./2; Ma=(Ma+Ma')./2; % Symmetrization

NestO2P=max(eig(Mp));
NestO2A=max(eig(Ma));

end
