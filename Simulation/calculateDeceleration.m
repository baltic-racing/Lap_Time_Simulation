%% calculateDeceleration.m
% Function to calculate the braking deceleration and brakebias.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

%% Calculates the braking deceleration with given tire parameters and brake setup. The BrakeBias can be set static and is adjusted accordingly or can be set dynamic (ABS)
function [FB_fl, FB_fr, FB_rl, FB_rr, aRev, BrakeBias, ABS] = calculateDeceleration(FB, m_tot, Fdr, FWXmax_fl, FWXmax_fr, FWXmax_rl, FWXmax_rr, BrakeBias_setup)
    
    ABS = 0; % Value for ABS status at beginning = 0 (no ABS activated)

    % Front Brake
    if (FB / 4 > FWXmax_fl)
        FB_fl = FWXmax_fl;
        ABS = 1; % Set ABS to on (force on tire needs to be reduced)
    else
        FB_fl = FB / 4;
    end

    if (FB / 4 > FWXmax_fr)
        FB_fr = FWXmax_fr;
        ABS = 1;
    else
        FB_fr = FB / 4;
    end

    % Rear Brake
    if (FB / 4 > FWXmax_rl)
        FB_rl = FWXmax_rl;
        ABS = 1;
    else
        FB_rl = FB / 4;
    end

    if (FB / 4 > FWXmax_rr)
        FB_rr = FWXmax_rr;
        ABS = 1;
    else
        FB_rr = FB / 4;
    end

    % Add Brakeforce on front and rear axle
    FB_f = FB_fl + FB_fr;
    FB_r = FB_rl + FB_rr;
    
    % Calculate Brake Bias in Percent Front
    BrakeBias = (FB_f) / ((FB_f)+(FB_r)) * 100;
    
    % Dynamic ideal Brake Bias (ABS) no BrakeBias as parameter
    if BrakeBias_setup ~= -1                   

        % Adjust Brakebias to setup (front brake threshold)
        if BrakeBias_setup > BrakeBias
            
            % Calculate new rear brake force
            FB_r = (-BrakeBias_setup / 100 * FB_f + FB_f) / BrakeBias_setup;
            
            % Adjust brake forces on left and right
            FB_right = FB_rl /(FB_rr + FB_rl);
            
            FB_rl = (1-FB_right) * FB_r;
            FB_rr = FB_right * FB_r;
            
        % Adjust Brakebias to setup (rear brake threshold)
        elseif BrakeBias_setup < BrakeBias
            
            % Calculate new front brake force
            FB_f = (BrakeBias_setup / 100 * FB_r)/(-BrakeBias_setup + 1);
            
            % Adjust brake forces on left and right
            FB_right = FB_fl /(FB_fr + FB_fl);
            
            FB_fl = (1-FB_right) * FB_f;
            FB_fr = FB_right * FB_f;
            
        end
            
        % Final BrakeBias
        BrakeBias = (FB_f) / ((FB_f)+(FB_r)) * 100;            
    end
    
    % Calculate Deacceleration
    aRev = (-Fdr -(FB_fl + FB_fr + FB_rl + FB_rr)) / m_tot;
end
