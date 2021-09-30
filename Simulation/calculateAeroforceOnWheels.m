%% calculateAeroforceOnWheels.m
% Function to calculate the aero force on each individual tire.
% rl = rear left, rr = rear right, fl = front left, fr = front right
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [dFWZrl_aero, dFWZrr_aero, dFWZfl_aero, dFWZfr_aero] = calculateAeroforceOnWheels(Faero, aero_ph, aero_pv)
%% Individual aerodynamic forces acting on each wheel
%   For calculation of wheel load transfer

        dFWZrl_aero = Faero/2*aero_ph;   % [N] Aerodynamic force on rear left wheel (Aerokraft auf linkes Hinterrad)
        dFWZrr_aero = Faero/2*aero_ph;   % [N] Aerodynamic force on rear right wheel (Aerokraft auf rechtes Hinterrad)
        dFWZfl_aero = Faero/2*aero_pv;   % [N] Aerodynamic force on front left wheel (Aerokraft auf linkes Vorderrad)
        dFWZfr_aero = Faero/2*aero_pv;   % [N] Aerodynamic force on front right wheel (Aerokraft auf rechtes Vorderrad)
            
end