%% calculateVtirestiff.m
% Calculates the interpolated tire stiffness.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [cZ_fl, cZ_fr, cZ_rl, cZ_rr] = calculateVtirestiff(Fz, cZ_tire, FWZ_fl, FWZ_fr, FWZ_rl, FWZ_rr)
%% Vertical tire stiffness

        cZ_fl = interp1(Fz,cZ_tire,FWZ_fl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_fr = interp1(Fz,cZ_tire,FWZ_fr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_rl = interp1(Fz,cZ_tire,FWZ_rl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_rr = interp1(Fz,cZ_tire,FWZ_rr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)       
end