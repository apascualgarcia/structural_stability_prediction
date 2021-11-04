function [r,e1]=eigenMax(D,E)
% This function simply gets the diagonal form of a matrix and of its
% eigenvectors and returns the maximum eigenvalue and the absolute value of
% the main eigenvector, making some controls to ensure the applicability of
% Perron Frobenius
[r,idx]=max(max(real(D))); % Identify the maximum eigenvalue
e1=E(:,idx);    % And eigenvector, it is already normalized, i.e. \sum_i (e^1_i)^2=1
minE1=min(e1);
maxE1=max(e1);
if(abs(minE1)<1e-20) minE1=0; end % Be aware of this condition, influences the next one
if(abs(maxE1)<1e-20) maxE1=0; end % Be aware of this condition, influences the next one
if(abs(sign(minE1)-sign(maxE1))==2) % This condition should hold
    fprintf('%s %f \n','   !! Error: Perron-Frobenius does not hold -- leaving the function');
    return;
end
e1=abs(e1); % To apply Perron-Frobenious should have the same sign, we choose a positive one
%%% hasta aquí

end