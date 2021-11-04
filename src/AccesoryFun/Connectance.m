

%**************************************%
%            Connectance.m             %
%**************************************%
% Computes the connectance of a given adjacency undirected matrix, defined
% as the fraction of non-empty cells.
% Input: adjacency matrix. Output: Connectance.
% Madrid, January 27th, 2014
% Bioinformatics Unit CBMSO, A.P-G.
%

function Conn=Connectance(input)

[m,n]=size(input);
adj=zeros(m,n);
%input=input./(max(max(input))); % Normalize values between zero and one
adj(input>0)=1; % Turn to one those cells of adj where you find in input positive values

Conn=sum(sum(adj))/(m*n);
end
