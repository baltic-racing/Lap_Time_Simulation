function [dFWZfl_y, dFWZfr_y, dFWZrl_y, dFWZrr_y] = calculateWheelloadLatDisp(h_COG, B, lh, lv, l, FVY) 
%% Dynamic wheel load displacement in lateral direction
%   For calculation of wheel load transfer

        dFWZfl_y = -h_COG/B*lh/l*FVY;   % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
        dFWZfr_y = h_COG/B*lh/l*FVY;    % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
        dFWZrl_y = -h_COG/B*lv/l*FVY;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
        dFWZrr_y = h_COG/B*lv/l*FVY;    % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
end