%% calculateSteeringData.m
% Calculates steering data such as delta and psi1 as well as slip angles.
%
% By Eric Dornieden / Patrick Siuta, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [delta, beta, psi1, alpha_f, alpha_r, alpha_fr, alpha_fl, alpha_rr, alpha_rl] = calculateSteeringData(wheelbase, R, lr, lf, vV, FWZ_fl, FWZ_rl)
    
    %% Check corner direction (left / right corner)
    if R > 0
        f = -1;
    else
        f = 1;
    end
    
    %% old
%     % Lenkwinkel, Schwimmwinkel, Gierrate, Schr�glaufwinkel  
%     delta = f*atan((wheelbase/1000)/sqrt(R^2-(lr/1000)^2));   % [rad] Lenkwinkel
%     beta = f*atan((lr/1000)/sqrt(R^2-(lr/1000)^2));  % [rad] Schwimmwinkel
%     psi1 = vV/R;                               % [rad/s] Gierrate
%     alpha_f = 180/pi*(delta-(lf/1000)/vV*psi1-beta);  % [�] Schr�glaufwinkel vorne
%     alpha_r = 180/pi*((lr/1000)/vV*psi1-beta);           % [�] Schr�glaufwinkel hinten
%     
%     alpha_fr = 0;
%     alpha_fl = 0;
%     alpha_rr = 0;
%     alpha_rl = 0;
    
    delta = f*atan((wheelbase/1000)/sqrt(R^2-(lr/1000)^2));         % [rad] steerung angel (Lenkwinkel)
    beta = f*atan(((-lr/1000)/wheelbase)*delta);                    % [rad] sideslip (Schwimmwinkel)
    psi1 = vV/R;                                                    % [rad/s] yaw rate(Gierrate)
    
    %% Slipangle front
    alpha_fl = 180/pi*(beta+(lf/1000)/vV-delta);                    % [degree] Slipangle rear left
    alpha_fr = 180/pi*(beta+(lf/1000)/vV-delta);                    % [degree] Slipangle rear left
    alpha_f = (alpha_fr+ alpha_fl)/2;                               % [degree] Slipangle rear left
    
    %% Slipangle rear
    alpha_rl = 180/pi*((1-(FWZ_fl/FWZ_rl)*alpha_f));                % [degree] Slipangle rear left
    alpha_rr = 180/pi*((1-(FWZ_fl/FWZ_rl)*alpha_f));                % [degree] Slipangle rear right
    alpha_r = (alpha_rr+ alpha_rl)/2;                               % [degree] Slipangle rear
end