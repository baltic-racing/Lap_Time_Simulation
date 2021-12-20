%% calculateWheelloadLatDisp.m
% Function to calculate the dynamic load transfer to left and right respectively
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [dFWZfl_y, dFWZfr_y, dFWZrl_y, dFWZrr_y, dFWZfl_geometric, dFWZfr_geometric, dFWZrl_geometric, dFWZrr_geometric, dFWZfl_elastic, dFWZfr_elastic, dFWZrl_elastic, dFWZrr_elastic] = calculateWheelloadLatDisp(h_COG, track_f, track_r, lr, lf, wheelbase, FVY, h_rc_f, h_rc_r) 
%% Dynamic wheel load displacement in lateral direction
%   For calculation of wheel load transfer

        %{
            Calculate total wheelload transfer:
            deltaW_total = (m * ay * h_com) / t
        %}

        %% Geometric wheel load transfer
        dFWZfl_geometric = -lr/wheelbase*FVY*h_rc_f/track_f;        % [N] Geometric wheel load transfer to front left wheel
        dFWZfr_geometric = lr/wheelbase*FVY*h_rc_f/track_f;         % [N] Geometric wheel load transfer to front right wheel
        dFWZrl_geometric = -lf/wheelbase*FVY*h_rc_r/track_r;        % [N] Geometric wheel load transfer to rear left wheel
        dFWZrr_geometric = lf/wheelbase*FVY*h_rc_r/track_r;         % [N] Geometric wheel load transfer to rear right wheel

        %% Elastic wheel load transfer
        dFWZfl_elastic = -lr/wheelbase*FVY*(h_COG-h_rc_f)/track_f;  % [N] Elastic wheel load transfer to front left wheel
        dFWZfr_elastic = lr/wheelbase*FVY*(h_COG-h_rc_f)/track_f;   % [N] Elastic wheel load transfer to front right wheel
        dFWZrl_elastic = -lf/wheelbase*FVY*(h_COG-h_rc_r)/track_r;  % [N] Elastic wheel load transfer to rear left wheel
        dFWZrr_elastic = lf/wheelbase*FVY*(h_COG-h_rc_r)/track_r;   % [N] Elastic wheel load transfer to rear right wheel

        %% Dynamic wheel load transfer
        dFWZfl_y = dFWZfl_geometric + dFWZfl_elastic;               % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
        dFWZfr_y = dFWZfr_geometric + dFWZfr_elastic;               % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
        dFWZrl_y = dFWZrl_geometric + dFWZrl_elastic;               % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
        dFWZrr_y = dFWZrr_geometric + dFWZrr_elastic;               % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
end