%% calculateTractiveForces.m
% Function to calculate the dynamic load transfer to front or rear respectively
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [dFWZfl_x, dFWZrr_x, dFWZrl_x, dFWZfr_x] = calculateWheelloadLongDisp(h_COG, m_ges, aVX, l)
%% Dynamic wheel load displacement in longitudinal direction
%   For calculation of wheel load transfer

        dFWZfl_x = -m_ges*aVX*h_COG/l/2;  % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
        dFWZfr_x = -m_ges*aVX*h_COG/l/2;  % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
        dFWZrl_x = m_ges*aVX*h_COG/l/2;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
        dFWZrr_x = m_ges*aVX*h_COG/l/2;   % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
        
end