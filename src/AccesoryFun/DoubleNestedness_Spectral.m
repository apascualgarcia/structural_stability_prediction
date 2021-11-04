%********************************************%
%     DoubleNestedness_Spectral.m            %
%********************************************%

 function DNSp = DoubleNestedness_Spectral(A);

% Computes the spectral radius of a mutualistic matrix following:
%  ``The ghost of nestedness in ecological networks''
% by Phillip P.A. Staniczenko, Jason C. Kopp and Stefano Allesina
% Nature Communications 4:1391 doi: 10.1038/ncomms2422 (2013)
%
% Bioinformatics Unit CBMSO
% A.P-G. Madrid, Dec 23rd, 2013 

[Sp,Sa]=size(A);
St=Sp+Sa;
adj=zeros(St,St);
adj(1:Sp,Sp+1:St)=A;
adj=adj+adj';
[V,D]=eig(adj); % Extract eigevectors and eigenvalues
[Lambda,index]=max(max(D)); % Identify the index
Lambda=Lambda/norm(adj,'fro'); % Normalize by the Frobenius norm
DNSp=Lambda;
 end
