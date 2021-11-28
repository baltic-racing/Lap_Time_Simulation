%% calculateSkidPad.m
% Calculates the skidpad time and skidpad velocity.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [t_skidpad, vV_skidpad] = calculateSkidPad(downforce_multiplier, c_l, A_S, rho_L, ConstantDownforce, c_l_DRS, DRS_status, m_tot, lr, lf, wheelbase, track_f, track_r, aero_ph, aero_pv, h_COG, GAMMA, TIRparam, FWZ_fl_stat, FWZ_fr_stat, FWZ_rl_stat, FWZ_rr_stat)
    FWYf = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
    FWYr = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
    FWYmax_f = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal übertragbare Querkraft Vorderachse)
    FWYmax_r = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal übertragbare Querkraft Hinterachse)
    vV = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)
   
    R = (15.25 + max(track_f,track_r) / 1000 / 2); % Skidpad
    aVX = 0; % Skidpad

    while  FWYf < FWYmax_f && FWYr < FWYmax_r && vV < 30

        vV = vV + 0.01;   % [m/s] Increaing vehicle speed (Erhöhen der Fahrzeuggeschwindigkeit)

        Faero = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vV, ConstantDownforce, c_l_DRS, DRS_status); % [N] Aerodynamic force

        FVY = m_tot*vV^2/R;    % [N] Centrifugal force (Zentrifugalkraft)

        aVY = vV^2/R;  % [m/s²] Lateral acceleration (Querbeschleunigung)

        % Lateral forces to be applied on front and rear axle (Aufzubringende Querkräfte an Vorder- und Hinterachse)
        FWYf = lr/wheelbase*abs(FVY);   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
        FWYr = lf/wheelbase*abs(FVY);   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

        % Wheel load transfer due to drag forces (Radlastverlagerung in Folge von Aerokräften) 
        [dFWZrl_aero, dFWZrr_aero, dFWZfl_aero, dFWZfr_aero] = calculateAeroforceOnWheels(Faero, aero_ph, aero_pv);

        % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in Längsrichtung = 0 angenommen)
        [dFWZfl_x, dFWZfr_x, dFWZrl_x, dFWZrr_x] = calculateWheelloadLongDisp(h_COG, m_tot, aVX, wheelbase); % Loads = 0 assumed

        % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
        [dFWZfl_y, dFWZfr_y, dFWZrl_y, dFWZrr_y] = calculateWheelloadLatDisp(h_COG, track_f, track_r, lr, lf, wheelbase, FVY);

        % Wheel loads (Radlasten)
        FWZ_fl = FWZ_fl_stat + dFWZfl_aero + dFWZfl_x + dFWZfl_y; % [N] Front left wheel load (Radlast vorne links)
        FWZ_fr = FWZ_fr_stat + dFWZfr_aero + dFWZfr_x + dFWZfr_y; % [N] Front right wheel load (Radlast vorne rechts)
        FWZ_rl = FWZ_rl_stat + dFWZrl_aero + dFWZrl_x + dFWZrl_y; % [N] Rear left wheel load (Radlast hinten links)
        FWZ_rr = FWZ_rr_stat + dFWZrr_aero + dFWZrr_x + dFWZrr_y; % [N] Rear right wheel load (Radlast hinten rechts)   

        % Maximum transmissible tire forces in longitudinal direction = 0 assumed (because longitudinal wheel loads = 0 assumed) 

        % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)    
        [FWYmax_f, FWYmax_r] = calculateLatTireforces(FWZ_fl, FWZ_fr,FWZ_rl, FWZ_rr, GAMMA, TIRparam);

    end

    vV_skidpad = vV;            % [m/s] Maximum speed for Skidpad
    
    t_skidpad = pi * R / vV;    % [s] tinme for Skidpad
end