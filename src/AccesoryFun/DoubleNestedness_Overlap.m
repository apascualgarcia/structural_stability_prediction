%********************************************%
%            DoubleNestedness_Overlap.m      %
%********************************************%

 function DNOv = DoubleNestedness_Overlap(adj)

% Computes a generalized nestedness measure that includes the specific
% mutualistic values of the matrix.
%
% Bioinformatics Unit CBMSO
% A.P-G. Madrid, Dec 23rd, 2013 
[m,n]=size(adj);

Nest=0;

NestA=0;
countA=0;
for i=1:m-1
    for j=i+1:m
        Nshare=sum(min(adj(i,:),adj(j,:)));                       
        Limit=sum(min(sort(adj(i,:)),sort(adj(j,:))));  
        NestA=NestA+Nshare/Limit;
        Nest=Nest+Nshare/Limit;
        countA=countA+1;
    end
end
NestA=NestA/countA;

NestP=0;
countP=0;
for i=1:n-1
    for j=i+1:n
        Nshare=sum(min(adj(:,i),adj(:,j)));
        Limit=sum(min(sort(adj(:,i)),sort(adj(:,j))));
        NestP=NestP+Nshare/Limit;
        Nest=Nest+Nshare/Limit;
        countP=countP+1;
    end
end
NestP=NestP/countP;

Nestedness=Nest/(countA+countP);
DNOv=(NestA+NestP)/2;

 end
