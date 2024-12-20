%% calculateSteeringData.m
% Calculates steering data such as delta and psi1 as well as slip angles.
%
% By Eric Dornieden / Patrick Siuta, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [delta, beta, psi1, alpha_f, alpha_r, alpha_fr, alpha_fl, alpha_rr, alpha_rl, delta_fl, delta_fr, delta_sw, ackermann, ackermannPercent] = calculateSteeringData(wheelbase, R, lr, lf, vV, FWZ_fl, FWZ_rl, FWZ_fr, FWZ_rr, track_f, cZ_fl, cZ_fr, cZ_rl, cZ_rr, alpha_f, alpha_r)
    
    alpha_f = alpha_f;
    alpha_r = alpha_r;

    alpha_fl = alpha_f;
    alpha_fr = alpha_f;
    alpha_rl = alpha_r;
    alpha_rr = alpha_r;


    %% Check corner direction (left / right corner)
    if R > 0
        f = -1;
    else
        f = 1;
    end

    % Convert rad to degree by multypling by pi/180 convert to rad by
    % multypling by 180/pi

    % Steering Angle,
    % Angle, 
    %delta = f*atan((wheelbase/1000)/sqrt(R^2-(lr/1000)^2));     % [rad] steering angle (Lenkwinkel)
    delta = wheelbase/1000/R; %+(alpha_f-alpha_r);              % [rad] steering angle (Formula 5.5.2 p. 110 / Performance Vehicle Dynamics, James Balkwill)
    beta = f*atan((lr/1000)/sqrt(R^2-(lr/1000)^2));             % [rad] Schwimmwinkel
    psi1 = vV/R;                                                % [rad/s] Gierrate
    % alpha_f = 180/pi*((beta+((lf/1000) *psi1 )/vV)-delta);      % [degree] Slipangle Front (Formula 5.4 S148 / Race Car Vehicle Dynamics) 
    % alpha_r = 180/pi*((lr/1000)/vV*psi1-beta);                  % [degree] Slipangle rear
    delta_sw = delta * 180/pi * 5.093;                          % [degree] steering angle at the steering wheel
    
    % alpha_fl = FWZ_fl * (vV^2 / R) / cZ_fl;

    %% Calculate ideal steering angles for inner and outer wheel
    if R > 0
        % R is postive == right direction corner == fr inner wheel and fl
        % outer wheel for the corner.
        delta_fl = atan((wheelbase/1000)/(R+(track_f/1000/2)));            % [rad] steering front left wheel (outer wheel) 
        delta_fr = atan((wheelbase/1000)/(R-(track_f/1000/2)));            % [rad] steering front left wheel (inner wheel)

        ackermann = atan(wheelbase/((wheelbase/tan(delta_fl))-track_f));   % 
        ackermannPercent = delta_fr/ackermann*100;                         % [%] calculated ackermmann in percent (for ideal steering always 100 (checked))
    else
        delta_fl = atan((wheelbase/1000)/(R-(track_f/1000/2)));            % [rad] steering front left wheel (inner wheel)
        delta_fr = atan((wheelbase/1000)/(R+(track_f/1000/2)));            % [rad] steering front left wheel (outer wheel)

        ackermann = atan(wheelbase/((wheelbase/tan(delta_fr))-track_f));
        ackermannPercent = delta_fl/ackermann*100;                         % [%] calculated ackermmann in percent (for ideal steering always 100 (checked))
    end   
    
    % %% Check fr and fl @PatrickSiuta
    % if FWZ_fr/FWZ_rr >1                                     % [degree] Slipangle rear right
    %     alpha_fr = ((1/(FWZ_fr/FWZ_rr))*alpha_f);
    % else
    %     alpha_fr = ((1-(FWZ_fr/FWZ_rr))*alpha_f);        
    % end
    % 
    % if FWZ_fl/FWZ_rl >1                                     % [degree] Slipangle rear right
    %     alpha_fl = ((1/(FWZ_fl/FWZ_rl))*alpha_f);
    % else
    %     alpha_fl = ((1-(FWZ_fl/FWZ_rl))*alpha_f);        
    % end
    % 
    % %% rr and rl
    % if FWZ_fr/FWZ_rr >1                                     % [degree] Slipangle rear right
    %     alpha_rr = ((1/(FWZ_fr/FWZ_rr))*alpha_r);
    % else
    %     alpha_rr = ((1-(FWZ_fr/FWZ_rr))*alpha_r);        
    % end
    % 
    % if FWZ_fl/FWZ_rl >1                                     % [degree] Slipangle rear left
    %     alpha_rl = ((1/(FWZ_fl/FWZ_rl))*alpha_r);
    % else
    %     alpha_rl = ((1-(FWZ_fl/FWZ_rl))*alpha_r);        
    % end

    
    
%     delta = f*atan((wheelbase/1000)/sqrt(R^2-(lr/1000)^2));         % [rad] steering angle (Lenkwinkel)
%     beta = f*atan(((-lr/1000)/wheelbase)*delta);                    % [rad] sideslip (Schwimmwinkel)
%     psi1 = vV/R;                                                    % [rad/s] yaw rate (Gierrate)
%     
%     %% Slipangle front
%     alpha_fl = 180/pi*(beta+(lf/1000)/vV-delta);                    % [degree] Slipangle rear left
%     alpha_fr = 180/pi*(beta+(lf/1000)/vV-delta);                    % [degree] Slipangle rear left
%     alpha_f = (alpha_fr+ alpha_fl)/2;                               % [degree] Slipangle rear left
%     
%     if alpha_f > 100 
%         alpha_f = 0;
%     end
%     
%     %% Slipangle rear
%     alpha_rl = 180/pi*((1-(FWZ_fl/FWZ_rl)*alpha_f));                % [degree] Slipangle rear left
%     alpha_rr = 180/pi*((1-(FWZ_fl/FWZ_rl)*alpha_f));                % [degree] Slipangle rear right
%     alpha_r = (alpha_rr+ alpha_rl)/2;                               % [degree] Slipangle rear
end