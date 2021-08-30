%% checkRunID.m
% Checks if the selected runID is valid.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function error = checkRunID(app, runID, saveFile)
    % check if a valid runID has been used
    if (runID > size(saveFile.t,2))     
        uialert(app.LaptimeSimulationUIFigure,'runID is greater than the number of runs in the save file','runID Error')
        error = 1;
    else
        error = 0;
    end         
end