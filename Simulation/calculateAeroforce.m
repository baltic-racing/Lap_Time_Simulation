%% calculateAeroforce.m
% Calculates the aerodynamic forces produced by the wings of the race car. 
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function Faero = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vV, ConstantDownforce, c_l_DRS, DRS_status)
%% Aerodynamic force
        if (DRS_status == 1)
            Faero = downforce_multiplier * ((c_l_DRS * A_S * rho_L * (vV)^2 / 2) + ConstantDownforce);  
        else
            Faero = downforce_multiplier * ((c_l * A_S * rho_L * (vV)^2 / 2) + ConstantDownforce);  
        end
end