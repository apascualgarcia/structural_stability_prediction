%********************************************%
%            DoubleNestedness_Cosine.m      %
%********************************************%

 function DNcos = DoubleNestedness_Cosine(adj)

% Computes a generalized nestedness measure that includes the specific
% mutualistic values of the matrix.
%
% Bioinformatics Unit CBMSO
% A.P-G. Madrid, Dec 23rd, 2013 
[m,n]=size(adj);

Limit=0;
Nest=0;
NestA=0;
countA=0;
for i=1:m-1;
    for j=i+1:m;
        Nshare=sum(adj(i,:).*adj(j,:)); 
        Limit=sqrt(sum(adj(i,:).^2)*sum(adj(j,:).^2));
        NestA=NestA+Nshare/Limit;
        Nest=Nest+Nshare/Limit;
        countA=countA+1;
    end
end
NestA=NestA/countA;
NestP=0;
countP=0;
for i=1:n-1;
    for j=i+1:n;
        countP=countP+1;
        Nshare=sum(adj(:,i).*adj(:,j)); 
        Limit=sqrt(sum(adj(:,i).^2)*sum(adj(:,j).^2));
        NestP=NestP+Nshare/Limit;
        Nest=Nest+Nshare/Limit;
    end
end
NestP=NestP/countP;
Nestedness=Nest/(countA+countP);
DNcos=(NestA+NestP)/2;

 end
