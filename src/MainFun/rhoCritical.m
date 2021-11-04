%% Function to compute rho critical

function [rhoP,rhoA]=rhoCritical(Conn,XiP,XiA,Sp,Sa)

        rho0A=rand; % starting conditions
        rho0P=rand;
        rhoA=rho0A;
        rhoP=rho0P;    
        key=true; keyA=true; keyP=true; % keys to control convergence
        k=0; 
        while(key) % start procedure to look for convergence
            k=k+1; % counter for plots 
            S0a=(1-rhoA)/rhoA;  
            S0p=(1-rhoP)/rhoP;
            SSA=Sa/(Sa+S0a);
            SSP=Sp/(Sp+S0p);
            rhoCp=(XiA-Conn^2*SSA)/(Conn-XiP*SSA);
            rhoCa=(XiP-Conn^2*SSP)/(Conn-XiA*SSP); % rho critical
            if(abs(rhoA-rhoCa)>eps) % test if there is convergence for animals
                rhoA=rhoCa; % rewrite and look for a new value
            else
                keyA=false; % convergence achieved
            end
            if(abs(rhoP-rhoCp)>eps) % test if there is convergence for plants
                rhoP=rhoCp;
            else
                keyP=false;
            end
            if((keyA == false)&&(keyP == false)) % if both achieve convergence
                key=false; % and look for the newt realization
            end
        end
end
