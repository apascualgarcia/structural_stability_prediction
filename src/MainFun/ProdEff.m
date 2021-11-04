function [Pa,Pp] = ProdEff(alphaA,alphaP,adj,Gamma,Rho,Beta0)




[m,n]=size(adj);

Smax=(1-Rho)/Rho;
% Terms for the inverse of the pure competition matrix
Term1=(1/Beta0)*1/(1-Rho);
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
            Inv_Comp_PureA(i,j)=Term1*(1-Term2m);
        else
            Comp_PureA(i,j)=Beta0*Rho;
            Inv_Comp_PureA(i,j)=Term1*(-Term2m);
        end
    end
end
        
             
for i=1:n % Plants
    for j=1:n
        if i==j
            Comp_PureP(i,j)=Beta0;
            Inv_Comp_PureP(i,j)=Term1*(1-Term2n);
        else
            Comp_PureP(i,j)=Beta0*Rho;
            Inv_Comp_PureP(i,j)=Term1*(-Term2n);
        end
    end
end

% Compute effective productivities

for i=1:m % Animals
    cum=0;
    for k=1:n
        for l=1:n
            cum=cum+adjA(i,k)*Inv_Comp_PureP(k,l)*alphaP(l);
        end
    end
    Pa(i)=alphaA(i)+cum;
end

for i=1:n % Animals
    cum=0;
    for k=1:m
        for l=1:m
            cum=cum+adjP(i,k)*Inv_Comp_PureA(k,l)*alphaA(l);
        end
    end
    Pp(i)=alphaP(i)+cum;
end


            