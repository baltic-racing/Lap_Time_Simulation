%% calculateGearbox.m
% Calculates the data of the gearbox and calculates parameters such as the
% motor and wheel rpm.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [rpm, gear, t_x] = calculateGearbox(gearbox, idleRPM, n_shift, n_downshift, vV, gr, gear, Rdyn_rl, Rdyn_rr, i_G, n_Mmax, t_x, rpm_in, t, t_last)
%calculateGearbox Calculates the motor rpm and the current gear and handles
%the shift dead time.
%   
%   rpm_in = last motor rpm (when i = 0 idle rpm or 0 (electric motor))
    Rdyn_r = (Rdyn_rl+Rdyn_rr)/2;                                   % Median of rear axle tire radius
    
    %% Gearbox
    if nargin == 15                                                 % Checks if number of input arguments is 14 or 11 (no time and rpm_in)   
        if gearbox                                                  % Checks if the car uses a gearbox
            num_gears = length(gr);                                 % highest gear -> number of gears in gearbox

            if rpm_in < n_downshift && gear > 1                     % Shift down
                gear = gear - 1;                    
                t_x = t-t_last;                                     % Dead time for shift time
            elseif rpm_in >= n_shift && gear < num_gears            % Shift one gear up when shifting rpm is reached
                gear = gear + 1;                                    
                t_x = t-t_last;                                     % step size of the time variable       
            end

            if t_x > 0                                              % Checks if gearshift dead time is active (>0)
                t_x = t_x+t-t_last;                                 % Adds past time to shift dead time
            end

            rpm = vV * 30 / pi * gr(gear) * i_G / Rdyn_r;           % [1/min] calculate rpm with gearbox    
        else                                                        % Else stament if there is no gearbox      
            rpm = vV * 30 / pi * i_G / Rdyn_r;                      % [1/min] Determine current motor speed 
        end
    else                                                            % if nargin == 11 (no time and no rpm_in (first iteration)
        if gearbox
            rpm = vV * 30 / pi * gr(gear) * i_G / Rdyn_r;           % [1/min] calculate rpm with gearbox 
        else                                                        % Determination of motor speed (no gearbox)
            rpm = vV * 30 / pi * i_G / Rdyn_r;                      % [1/min] Determine current motor speed 
        end
        
        t_x = 0;                                                    % Set shift dead time to zero for first iteration
    end
    
    
    %% RPM Limiter
    if rpm >= n_Mmax                                                % RPM-Limiter
        rpm = n_Mmax;                                               % if RPM is above limiter rpm set to RPM to limiter RPM
    elseif rpm < idleRPM                                            % Checks if rpm droppes below idle rpm (anti stall)
        rpm = idleRPM;                                              % Sets rpm to idle rpm
    end  
end