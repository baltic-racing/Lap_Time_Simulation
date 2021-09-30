%% calculateDynRadii.m
% Calculates the dynamic tire radius.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [Rdyn_fl, Rdyn_fr, Rdyn_rl, Rdyn_rr] = calculateDynRadii(R0, FWZ_vl, FWZ_vr, FWZ_hl, FWZ_hr, cZ_vl, cZ_vr, cZ_hr, cZ_hl)
%% Dynamic tire radii calculation

        Rdyn_fl = R0 - FWZ_vl/cZ_vl;    % [m] Dynamic front left tire radius (Dynamischer Reifenradius)
        Rdyn_fr = R0 - FWZ_vr/cZ_vr;    % [m] Dynamic front right tire radius (Dynamischer Reifenradius)
        Rdyn_rl = R0 - FWZ_hl/cZ_hl;    % [m] Dynamic rear left tire radius (Dynamischer Reifenradius)
        Rdyn_rr = R0 - FWZ_hr/cZ_hr;    % [m] Dynamic rear right tire radius (Dynamischer Reifenradius)        
end