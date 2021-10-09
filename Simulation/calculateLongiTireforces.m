%% calculateLongiTireforces.m
% Calculates the longitudinal tire forces.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [FWXmax_fl, FWXmax_fr, FWXmax_rl, FWXmax_rr, FWXmax_f, FWXmax_r] = calculateLongiTireforces(FWZ_vl, FWZ_vr,FWZ_hl, FWZ_hr, GAMMA, TIRparam, alpha_f, alpha_r)
%% Maximum transmissible tire forces in longitudinal direction
    if nargin == 8
        FWXmax_fl = max(abs(MF52_Fx_cs(alpha_f,FWZ_vl,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
        FWXmax_fr = max(abs(MF52_Fx_cs(alpha_f,FWZ_vr,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
        FWXmax_rl = max(abs(MF52_Fx_cs(alpha_r,FWZ_hl,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
        FWXmax_rr = max(abs(MF52_Fx_cs(alpha_r,FWZ_hr,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)        
    else
        FWXmax_fl = max(abs(MF52_Fx_cs(0,FWZ_vl,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
        FWXmax_fr = max(abs(MF52_Fx_cs(0,FWZ_vr,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
        FWXmax_rl = max(abs(MF52_Fx_cs(0,FWZ_hl,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
        FWXmax_rr = max(abs(MF52_Fx_cs(0,FWZ_hr,GAMMA,0:0.01:0.2,TIRparam)));   % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)        
    end
    
    FWXmax_f = FWXmax_fl + FWXmax_fr;                                       % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
    FWXmax_r = FWXmax_rl + FWXmax_rr;                                       % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
end