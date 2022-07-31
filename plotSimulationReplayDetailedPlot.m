%% plotSimulationReplayDetailedPlot.m
% Plots a single plot of the active variable from the adjacent DropDown
% menu in the simulation replay GUI.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2022, Baltic Racing, all rights reserved.

function plotSimulationReplayDetailedPlot(app, runNumber, result, dropDownValue) 

    % If dropDownValue is a char return (No save file has been loaded)
    if isa(dropDownValue,"char")
        return
    end
    
    % Plot dropDownValue data with annotation
    figure("Name",'Detailed Plot',"NumberTitle","off")
    plot(result(runNumber).Track(1:end-1,4), dropDownValue)
    xlabel("Distance [m]")
    annotation('textbox',[.2 .8 .1 .1],'String',['min: ' num2str(min(dropDownValue)) ' max: ' num2str(max(dropDownValue))])
    grid on
end