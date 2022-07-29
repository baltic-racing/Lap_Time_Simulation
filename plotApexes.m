%% updateReplayData.m
% plot the apexes on the track in the simulation replay
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2022, Baltic Racing, all rights reserved.

function plotApexes(app, runNumber, result)
    delete(app.Apexes);
            
    apex = (zeros(1,length(result(runNumber).ApexIndexes)));
    
    % Creates an array of every apex on the track
    k = 1;
    for i = 1:length(result(runNumber).Track(:,1))
        if ismember(i,result(runNumber).ApexIndexes)
            apex(k,1) = result(runNumber).Track(i,1);
            apex(k,2) = result(runNumber).Track(i,2);
            isApex(i) = 1;
            k = k + 1;
        else
            isApex(i) = 0;
        end
    end
    
    if app.DrawApexesCheckBox.Value % Only draw apexes if drawApexes checkbox is checked
        hold(app.UIAxes3,'on');
        app.Apexes = scatter(app.UIAxes3,apex(1:end,1),apex(1:end,2),150,'r'); % Draws points for every apex on the track
    end
end