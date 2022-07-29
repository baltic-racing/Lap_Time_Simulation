%% plotSimulationReplayDetailedPlot.m
% Plots a single plot of the active variable from the adjacent DropDown
% menu in the simulation replay GUI.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2022, Baltic Racing, all rights reserved.

function plotSimulationReplayDetailedPlot(app, runNumber, result, dropDownValue) 
    figure("Name",'Detailed Plot',"NumberTitle","off")
    plot(result(runNumber).Track(1:end-1,4), dropDownValue)
    xlabel("Distance [m]")
    grid on
end