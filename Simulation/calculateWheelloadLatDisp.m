function [dFWZfl_y, dFWZfr_y, dFWZrl_y, dFWZrr_y] = calculateWheelloadLatDisp(h_COG, track_f, track_r, lr, lf, wheelbase, FVY) 
%% Dynamic wheel load displacement in lateral direction
%   For calculation of wheel load transfer

        dFWZfl_y = -h_COG/track_f*lr/wheelbase*FVY;   % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
        dFWZfr_y = h_COG/track_f*lr/wheelbase*FVY;    % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
        dFWZrl_y = -h_COG/track_r*lf/wheelbase*FVY;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
        dFWZrr_y = h_COG/track_r*lf/wheelbase*FVY;    % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
end