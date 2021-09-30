%% calculateBrakeTorque.m
% Function to calculate the braking torque and other brake parameters.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [BT_fl, BT_fr, BT_rl, BT_rr] = calculateBrakeTorque(vV, n_wheel, BrakePressure, boreFront, boreRear)
    
    R_m = (R_o + R_i) / 2;

    friction_s = 0.35;
    friciton_k = 0.35;
    
    if n_wheel(i) > 0 
        BT_fl = (friciton_k * BrakePressure * pi * boreFront^2 * R_m * 2) / 4;
        BT_fr = (friciton_k * BrakePressure * pi * boreFront^2 * R_m * 2) / 4;
        BT_rl = (friciton_k * BrakePressure * pi * boreRear^2 * R_m * 2) / 4;
        BT_rr = (friciton_k * BrakePressure * pi * boreRear^2 * R_m * 2) / 4;
    else
        BT_fl = (friciton_s)/4
    end
end