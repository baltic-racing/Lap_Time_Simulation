%% calculateAxleloads.m
% Calculates the axle loads of the race car (front and rear axle and sum). 
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [FWZr, FWZf, FWZges] = calculateAxleloads(FWZ_rl, FWZ_rr, FWZ_fl, FWZ_fr)
%% Axle loads calculation
        FWZr = FWZ_rl + FWZ_rr;  % [N] Rear axle load (Achslast hinten)
        FWZf = FWZ_fl + FWZ_fr;  % [N] Front axle load (Achslast vorne)
        FWZges = FWZr + FWZf;    % [N] Total axle and wheel load (Gesamtachs- und Radlast)       
end