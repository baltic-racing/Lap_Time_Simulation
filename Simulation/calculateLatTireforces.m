%% calculateLatTireforces.m
% Calculates the latitudinal tire forces.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [FWYmax_fl, FWYmax_fr, FWYmax_rl, FWYmax_rr, FWYmax_f, FWYmax_r, alpha_f, alpha_r] = calculateLatTireforces(FWZ_vl, FWZ_vr, FWZ_hl, FWZ_hr, GAMMA, TIRparam, alpha_f, alpha_r)
%% Maximum transmissible tire forces in lateral direction
    if nargin == 8
        FWYmax_fl = max(abs(MF52_Fy_cs(alpha_f,FWZ_vl,GAMMA,0,TIRparam))); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
        FWYmax_fr = max(abs(MF52_Fy_cs(alpha_f,FWZ_vr,GAMMA,0,TIRparam))); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
        FWYmax_rl = max(abs(MF52_Fy_cs(alpha_r,FWZ_hl,GAMMA,0,TIRparam))); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
        FWYmax_rr = max(abs(MF52_Fy_cs(alpha_r,FWZ_hr,GAMMA,0,TIRparam))); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)
    else  
        x = 0:0.025:12; % Slip angle range [deg]

        % Front Left
        FWYmax_fl = abs(MF52_Fy_cs(x,FWZ_vl,GAMMA,0,TIRparam)); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)

        % Find the maximum output and its corresponding index
        [FWYmax_fl, index] = max(FWYmax_fl); 

        % Save slip angle
        alpha_fl = x(index);
        
        % Front Right
        FWYmax_fr = abs(MF52_Fy_cs(x,FWZ_vr,GAMMA,0,TIRparam)); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)

        % Find the maximum output and its corresponding index
        [FWYmax_fr, index] = max(FWYmax_fr); 

        % Save slip angle
        alpha_fr = x(index);

        % Rear Left
        FWYmax_rl = abs(MF52_Fy_cs(x,FWZ_hl,GAMMA,0,TIRparam)); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)

        % Find the maximum output and its corresponding index
        [FWYmax_rl, index] = max(FWYmax_rl); 

        % Save slip angle
        alpha_rl = x(index);

        % Rear Right
        FWYmax_rr = abs(MF52_Fy_cs(x,FWZ_hr,GAMMA,0,TIRparam)); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)

        % Find the maximum output and its corresponding index
        [FWYmax_rr, index] = max(FWYmax_rr); 

        % Save slip angle
        alpha_rr = x(index);

        alpha_f = max(alpha_fr, alpha_fl);
        alpha_r = max(alpha_rr, alpha_rl);
    end
    
    FWYmax_f = FWYmax_fl + FWYmax_fr;    % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
    FWYmax_r = FWYmax_rl + FWYmax_rr;    % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
end