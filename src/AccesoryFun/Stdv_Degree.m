

function Return = Stdv_Degree(input)
% Given a matrix, it turns it into a binary matrix (positive values will be
% equal to one) and it computes the standard deviation of the degrees,
% considering both rows and columns.

% Transform an arbitrary input into a binary matrix:
[m,n]=size(input);
adj=zeros(m,n);
input=input./(max(max(input))); % Normalize values between zero and one
adj(input>0)=1; % Turn to one those cells of adj where you find in input positive values

Drow=sum(adj,1); Dcol=sum(adj,2)'; % Degrees
Deg=cat(2,Drow,Dcol);
Return=std(Deg);

end