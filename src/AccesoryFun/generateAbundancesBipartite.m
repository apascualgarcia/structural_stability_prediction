function [Np,Na]=generateAbundancesBipartite(distro,n,m,N0p,N0a,Nwidth)
    % This function generates abundances for a bipartite system either
    % following a uniform distribution or a lognormal one. Nwidth controls
    % either the width of the uniform distribution:
    % N=N0*(1+eps*Nwidth) with eps \in [-0.5,0.5]
    % or the standard deviation of a N(0,Nwidth) from which a lognormal
    % distribution is generated, with variance  [\exp(Nwidth ^{2})-1]*\exp(2 +Nwidth ^{2})
    if(strcmp(distro,'uniform'))
        Np=N0p.*(1+(rand(1,n)-0.5)*Nwidth);
        Na=N0a.*(1+(rand(1,m)-0.5)*Nwidth);
    else
        Np=lognormalGenerator(n,Nwidth);
        Na=lognormalGenerator(m,Nwidth);
    end
end
