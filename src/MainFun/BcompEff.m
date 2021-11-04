%***********************************************************%
%    Normalized Effective Competition Matrix.m (BcompEff.m) %
%***********************************************************%



function [Ba,Bp] = BcompEff(adj,Gamma,Rho,Beta0)


[m,n]=size(adj);

Smax=(1-Rho)/Rho;
% Terms for the inverse of the pure competition matrix
Term1=(1/Beta0)*(1-Rho);
Term2m=1/(m+Smax); % For species m
Term2n=1/(n+Smax); % For species n
% Apply a mean field value to the binary matrix
adjA=Gamma*adj; adjP=Gamma*adj';
Ca=zeros(m); Cp=zeros(n);

% Compute pure competition matrix and its inverse

for i=1:m % Animals
    for j=1:m
        if i==j
            Comp_PureA(i,j)=Beta0;
            Inv_Comp_PureA(i,j)=Term1*(1+Term2m);
        else
            Comp_PureA(i,j)=Beta0*Rho;
            Inv_Comp_PureA(i,j)=Term1*Term2m;
        end
    end
end
        
             
for i=1:n % Plants
    for j=1:n
        if i==j
            Comp_PureP(i,j)=Beta0;
            Inv_Comp_PureP(i,j)=Term1*(1+Term2n);
        else
            Comp_PureP(i,j)=Beta0*Rho;
            Inv_Comp_PureP(i,j)=Term1*Term2n;
        end
    end
end

% Compute the effective competition matrix and its normalized matrix.


for i=1:m % Animals
    for j=1:m
        Ca(i,j)=Comp_PureA(i,j);
        for k=1:n % You must use the animals inverse competition matrix
            for l=1:n
                Ca(i,j)=Ca(i,j)-adjA(i,k)*Inv_Comp_PureP(k,l)*adjP(l,j);
            end
        end
    end
end
for i=1:m % Animals
    for j=1:m
        Ba(i,j)=Ca(i,j)/(sqrt(Ca(i,i)*Ca(j,j)));
    end
end
for i=1:n % Plants
    for j=1:n
        Cp(i,j)=Comp_PureP(i,j);
        for k=1:m % You must use the plants inverse competition matrix
            for l=1:m
                Cp(i,j)=Cp(i,j)-adjP(i,k)*Inv_Comp_PureA(k,l)*adjA(l,j);
            end
        end
    end
end
for i=1:n % Plants
    for j=1:n
        Bp(i,j)=Cp(i,j)/(sqrt(Cp(i,i)*Cp(j,j)));
    end
end
% imagesc(Ca);
% figure
% imagesc(Ba);
end
                    
                
                
                
                