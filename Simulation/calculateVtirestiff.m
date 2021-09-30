%% calculateVtirestiff.m
% Calculates the interpolated tire stiffness.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [cZ_fl, cZ_fr, cZ_rl, cZ_rr] = calculateVtirestiff(Fz, cZ_tire, FWZ_vl, FWZ_vr, FWZ_hl, FWZ_hr)
%% Vertical tire stiffness

        cZ_fl = interp1(Fz,cZ_tire,FWZ_vl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_fr = interp1(Fz,cZ_tire,FWZ_vr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_rl = interp1(Fz,cZ_tire,FWZ_hl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_rr = interp1(Fz,cZ_tire,FWZ_hr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)       
end