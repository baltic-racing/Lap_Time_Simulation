%% calculateRollMoment.m
% Calculates the Roll Moment of the racecar for the front and the rear
% axle.
%
% Documentation for roll centre and roll moment calculation:
% https://suspensionsecrets.co.uk/roll-centre-and-roll-moment/
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [rollMoment_f, rollMoment_r, rollAngleChassis] = calculateRollMoment(aVY, m_tot, h_cog, m_pr, h_rc_f, h_rc_r, g, x_va, x_cog, wheelbase)
    
    % h_rc_f = height roll center front [mm]
    % h_rc_r = height roll center rear [mm]

    roll_stiffness_f = 871; % Roll stifness front axle [Nm/deg]
    roll_stiffness_r = 864; % Roll stifness rear axle [Nm/deg]
    roll_stiffness = roll_stiffness_f + roll_stiffness_r; % overall Roll stifness [Nm/deg]

    % Calculate distance from roll center to cog
    alpha = tan(abs(h_rc_f-h_rc_r)/wheelbase);
    h_rc_cog = alpha*(x_cog - x_va) + min(h_rc_f, h_rc_r);
    h_rc = h_cog-h_rc_cog;

    % Sprung mass should be used here! GET SPRUNG UNSPRUNG MASS OF CAR
    m_r = m_tot * m_pr/100;                                 % Mass on rear axle [kg]
    m_f = m_tot - m_r;                                      % Mass on front axle [kg]
    
    rollMoment_f = m_f*aVY*g*(h_cog-h_rc_f)/1000;             % Rollmoment on front axle [Nm]
    rollMoment_r = m_r*aVY*g*(h_cog-h_rc_r)/1000;             % Rollmoment on rear axle [Nm]
    
    %rollMoment = m_tot*aVY*g*h_rc/1000;

    rollAngleChassis = (rollMoment_f + rollMoment_r) / roll_stiffness; % 2000 Platzhalter!
end