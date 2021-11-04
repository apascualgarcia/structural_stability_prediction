%********************************************%
%            Assortativity_Und.m             %
%********************************************%
% Computes the assortativity of a given bipartite adjacency undirected matrix with
% the Pearson correlation coefficient, following the measure proposed by
% Newman (Assortativity mixing in networks, PRL 2002, equation 4 arxiv
% version).
%
% Input: adjacency matrix. Output: assortativity coeff.
% Madrid, 30 dec, 2013
% Bioinformatics Unit CBMSO, A.P-G.
%

function r=Assortativity_Und(input)

% Transform an arbitrary input into a binary matrix:
[m,n]=size(input);
adj=zeros(m,n);
input=input./(max(max(input))); % Normalize values between zero and one
adj(input>0)=1; % Turn to one those cells of adj where you find in input positive values
% Compute degrees
degcol=sum(adj,1); % degree for columns (sum row values)
degrow=sum(adj,2); % degree for rows (sum column values)
M=sum(sum(adj)); % total number of links
% Compute assortativity
term1=0;
term2=0;
term3=0;
for i=1:m
    for j=1:n
        if(adj(i,j) == 0) continue; end
        term1=term1+degrow(i)*degcol(j);
        term2=term2+degrow(i)+degcol(j);
        term3=term3+degrow(i)^2+degcol(j)^2;
    end
end
term1=term1/M;
term2=(term2/(2*M))^2;
term3=term3/(2*M);

r=(term1-term2)/(term3-term2);

end
        


