function updateSimulationReplayGUI(app, resultNumbers)
    value = app.Slider_2.Value;                                             % Value of slider (current time)
    
    result = app.results.(['result' num2str(resultNumbers(1))]);
    result2 = [];

    if length(resultNumbers) > 1
        result2 = app.results.(['result' num2str(resultNumbers(2))]);
        numDataPoints2 = length(result2.t);
    end

    %% DEBUG!
    runNumber = 1;
    % ToDo: dynamic borders for runNumber based on .mat size
    
    if isempty(result)
        return
    end
                                                
    t = result(runNumber).t;                                                % Time variable of result    
    ind = interp1(t,1:length(t),value,'nearest');                           % Find index of the time nearest to the time of the slider
    
    if ind == length(t)
        ind = ind-1;
    end
    
    %% Load Variables into table
    fieldNames = fieldnames(result);
    numVariables = numel(fieldNames);
    numDataPoints1 = length(result.t);
    
    varNames = {};
    
    % Get run variables
    for i = 1:numVariables
        fieldName = fieldNames{i};  % Get the current field name
    
        if length(result.(fieldName)) >= numDataPoints1-2
            varNames{(end+1), 1} = fieldName;
            varNames{(end), 2} = num2str(result.(fieldName)(ind));
            varNames{(end), 3} = [];
        end

        if ~isempty(result2) && length(result2.(fieldName)) >= numDataPoints2-2
            varNames{(end), 3} = num2str(result2.(fieldName)(ind));
        end
    end
    
    % Update table
    set(app.UITable, 'Data', varNames, 'ColumnName', {'Name', 'Run 1', 'Run 2'})

    app.DRSStatusLamp.Enable = result.DRS_status(ind);                      % Change status of drs led dependend on DRS_status
    app.RPMGauge.Value = result.vV(ind).*3.6;                               % ToDo: add context menu option for m/s, km/h and mph
    app.RPMGauge_2.Value = result.ni(ind).*3.6;                             % ToDo: add context menu option for m/s, km/h and mph
    
    app.GearLabel.Text = num2str(result.gearSelection(ind));
    app.SpeedLabel.Text = "Speed: " + result.vV(ind) + " [m/s]";
    app.RPMLabel.Text = "RPM: " + round(result.ni(ind));
    app.TimeLabel.Text = "Time: " + result.t(ind) + " [s]";
    app.DistanceLabel.Text = "Distance: " + round(result.s(ind), 2) + " [m]";
    app.RadiusLabel.Text = "Radius: " + round(result.R(ind), 2) + " [m]";

    %xline = xline(app.UIAxes_3, ind, Color='red')
    delete(app.xline);
    app.xline = xline(app.UIAxes_3,ind,'r');

    %rotatedImage = imrotate(app.steeringImage, app.Slider_2.Value);
    %app.Image4_2.ImageSource(rotatedImage);
    
    
    %imshow(rotatedImage,'Parent',app.UIAxes_3);
end