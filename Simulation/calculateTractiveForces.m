%% calculateTractiveForces.m
% Function to calculate the tractive forces on the wheels
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [FVX_fl, FVX_fr, FVX_rl, FVX_rr, FVX_r, FVX_f, TC_f, TC_r, t_x] = calculateTractiveForces(Mi, num_motors, i_G, gr, Rdyn_fl, Rdyn_fr, Rdyn_rl, Rdyn_rr, t_x, gear, FWXmax_f, FWXmax_r, t_shift)
    % Calculation of tractive forces 
    if t_x == 0 || t_x >= t_shift
        if num_motors == 4
            FVX_fl = Mi/num_motors*i_G*gr(gear)/Rdyn_fl;    % [N] Tractive Force on front left wheel (AWD) 
            FVX_fr = Mi/num_motors*i_G*gr(gear)/Rdyn_fr;    % [N] Tractive Force on front right wheel (AWD) 
            FVX_f = FVX_fl + FVX_fr;                        % [N] Traction on rear axle (Zugkraft an der Hinterachse)

            FVX_rl = Mi/num_motors*i_G*gr(gear)/Rdyn_rl;    % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
            FVX_rr = Mi/num_motors*i_G*gr(gear)/Rdyn_rr;    % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
            FVX_r = FVX_rl + FVX_rr;                        % [N] Traction on rear axle (Zugkraft an der Hinterachse)

            if FVX_f > FWXmax_f                             % Limiting the tractive force to the traction limit front axle
                FVX_f = FWXmax_f;
                TC_f = 1;                                   % Traction control "on" (Traktionskontrolle "an")
            else
                TC_f = 0;
            end 
        else
            FVX_rl = Mi/2*i_G*gr(gear)/Rdyn_rl;             % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
            FVX_rr = Mi/2*i_G*gr(gear)/Rdyn_rr;             % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
            FVX_r = FVX_rl + FVX_rr;                        % [N] Traction on rear axle (Zugkraft an der Hinterachse)
            
            FVX_fl = 0;
            FVX_fr = 0;
            FVX_f = 0; 
            
            TC_f = 0;
          
        end

        if FVX_r > FWXmax_r                                 % Limiting the tractive force to the traction limit rear axle
            FVX_r = FWXmax_r;
            TC_r = 1;                                       % Traction control "on" (Traktionskontrolle "an")
        else
            TC_r = 0;
        end
        
        if t_x >= t_shift                                   % Reset t_x if higher than shift time
            t_x = 0;
        end
    else
        FVX_fl = 0;                                         % [N] Tractive Force on front left wheel (AWD) 
        FVX_fr = 0;                                         % [N] Tractive Force on front right wheel (AWD) 
        FVX_f = 0;                                          % [N] Tractive Force on rear axle (Zugkraft an der Hinterachse)

        TC_f = 0;
        
        FVX_rl = 0;                                         % [N] Tractive Force on rear left wheel (Zugkraft an linkem Hinterrad)
        FVX_rr = 0;                                         % [N] Tractive Force on rear right wheel (Zugkraft an rechtem Hinterrad)
        FVX_r = 0;                                          % [N] Tractive Force on rear axle (Zugkraft an der Hinterachse)
        
        TC_r = 0;
    end 
end