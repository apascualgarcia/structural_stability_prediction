function [GammaEffP,GammaEffA] = EffGammaMatrix(matrixP,hP,hA,Np,Na)

% This function is an extended version of the function
% EffNestedness.m, it includes the biomasses as input and
% it returns the effective gamma matrices, instead of the effective
% nestedness.
% The following function takes a list of gamma matrices and computes the
% effective gamma matrix. The effective gamma matrix is the Lotka-Volterra
% mutualistic matrix, i.e. the matrix we would consider in a linearized
% system (therefore close enough to the equilibrium fixed point). In this
% way, if the perturbation introduced is small, we expected a behaviour
% for the non linear system similar to
% the one of the linear system with the effective matrix.
% This effective matrix is estimated through the derivative of the mutualistic
% term with respect to the variation in the abundances, evaluated at the stationary
% state.
% Input:
% The resultant term for the effective matrix requires for its computation
% the adjacency matrix (including the gamma parameters) and the handling times. The
% steady state abundance values are estimated internally, but you will need
% to fix other parameters considered in the model.
% Output:
% Effective gamma matrices for plants and animals and/or Effective
% nestedness.
% Madrid, February 06th, 2014
% Bioinformatics Unit CBMSO, A.P-G.
%


[Sp,Sa]=size(matrixP);
matrixA=matrixP';

gammasumP=matrixP*Na; %sum(matrixP,2);
gammasumA=matrixA*Np; %sum(matrixA,2);

for i=1:Sp
    for j=1:Sa
        GammaEffP(i,j)=matrixP(i,j)/(1+hP*gammasumP(i))^2;
    end
end
for i=1:Sa
    for j=1:Sp
        GammaEffA(i,j)=matrixA(i,j)/(1+hA*gammasumA(i))^2;
    end
end

end
