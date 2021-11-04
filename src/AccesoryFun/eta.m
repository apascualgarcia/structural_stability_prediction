function eta1=eta(Alpha,e1)
% Project a vector on the main eigenvector and retrieve the 
% minimum of the quotient of the vector times the eigenvector
        Alpha1e=Alpha'*e1; % control this transposition and the one made for rand
        % Is it needed the absolute value here?
        eta_i=Alpha./(e1.*Alpha1e); %  Compute the vulnerability of species i
        eta1=min(eta_i);
end