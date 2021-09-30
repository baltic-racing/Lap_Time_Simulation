%% calculateLatTireforces.m
% Calculates the latitudinal tire forces.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [FWYmax_fl, FWYmax_fr, FWYmax_rl, FWYmax_rr, FWYmax_f, FWYmax_r] = calculateLatTireforces(FWZ_vl, FWZ_vr,FWZ_hl, FWZ_hr, GAMMA, TIRparam)
%% Maximum transmissible tire forces in lateral direction

    FWYmax_fl = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_vl,GAMMA,0,TIRparam))); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
    FWYmax_fr = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_vr,GAMMA,0,TIRparam))); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
    FWYmax_rl = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_hl,GAMMA,0,TIRparam))); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
    FWYmax_rr = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_hr,GAMMA,0,TIRparam))); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)

    FWYmax_f = FWYmax_fl + FWYmax_fr;    % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
    FWYmax_r = FWYmax_rl + FWYmax_rr;    % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
end