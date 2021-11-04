
    function etaTmp=delta2Eta(Delta0)
    % This function is included in
    % deltaCritical_Competition_MedianWithSolver.m, this version never
    % worked properly (see that function for details)
        global DeltaVec
        global etaVec
        global etaCritical
        % Handle function to use in fzero. Returns the position of Delta0
        % in the vec DeltaVec, and look for its response variable in etaVec
        k=DeltaVec==Delta0;
        etaTmp=etaVec(k)-etaCritical;
    end