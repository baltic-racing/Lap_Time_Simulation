%% calculateVehicleResistancesForces.m
% Calculates all resistance forces of the race car.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [FR, FL, Fdr, FVY, aVX, aVY] = calculateVehicleResistancesForces(k_R, FWZges, rho_L, vV, c_d, A_S, m_tot, R, FVX, FVX_f, c_d_DRS, DRS_status, rpm, n_Mmax, FB_fl, FB_fr, FB_rl, FB_rr)
%% Driving resistances and vehicle

        FR = k_R*FWZges;                % [N] Rolling resistance (Rollwiderstand)
        
        if DRS_status == 1
            FL = rho_L*(vV^2/2)*c_d_DRS*A_S;
        else
            FL = rho_L*(vV^2/2)*c_d*A_S;      % [N] Air Resistance/Drag (Luftwiderstand)
        end

        Fdr = FR + FL;                  % [N] Total resistance (Gesamtwiderstand)
        FVY = m_tot*vV^2/R;             % [N] Centrifugal force (Zentrifugalkraft)
        
        % Check if rpm limiter is reached if so don't accelerate further
        if rpm >= n_Mmax
            aVX = 0;
        else
            aVX = ((FVX+FVX_f)-Fdr-(FB_fl+FB_fr+FB_rl+FB_rr))/m_tot;       % [m/s²] Longitudinal acceleration (Längsbeschleunigung)
        end      
            
        aVY = vV^2/R;                   % [m/s²] Lateral acceleration (Querbeschleunigung)
end