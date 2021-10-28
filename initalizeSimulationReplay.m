%% initalizeSimulationReplay.m
% Initalizes the Simulation Replay GUI / Loads save file and sets inital
% conditions for the replay.
%
% By Eric Dornieden / Sarah Jaeschke, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [runNumber, result] = initalizeSimulationReplay(app)
    %try
        [file_run_3, path] = uigetfile('*.mat');                                    % Load the file name of the selected save file.

        result = load(path+"/"+file_run_3, '-mat');                               % open the .mat file with the given name.
        
        try        
            size(result(1).Track);
            app.RunNumberSpinner.Enable = 'off';
            app.RunNumberSpinner.Value = 1;
        catch 
            result = result.result;
            app.RunNumberSpinner.Limits = [1 size(result, 2)];   
            app.RunNumberSpinner.Enable = 'on';
        end
              
        runNumber = app.RunNumberSpinner.Value; 

        if app. DarkModeCheckBox.Value                                               % Check if Dark Mode is active, if so plot line in white else plot line in black.
            plot(app.UIAxes3,result(runNumber).Track(:,1),result(runNumber).Track(:,2),'w')       % Plot Track in white.
        else
            plot(app.UIAxes3,result(runNumber).Track(:,1),result(runNumber).Track(:,2),'k')       % Plot Track in black.
        end

        title(app.UIAxes3,'','FontSize',12);                                        % Set the title of the UIAxes to an empty String
        app.UIAxes2.LineWidth = 1.5;                                                % Set the Line width of the plot

        app.Slider.Limits = [0 result(runNumber).t(end)];                                    % Set maximum of the time slider to the overall laptime.
        app.Slider.Value = 0;                                                       % Reset Slider Value (Start Replay at the beginning of the lap)
        
        var1 = result(runNumber).vV/max(result(runNumber).vV);
        
        plot(app.UIAxes4,var1,'w')       % Plot Track in white.
        
        %app.SpeedLabel.Text = "Speed: " + num2str(result(runNumber).vV(ceil(app.Slider.Value)));
    %catch error
        %writeToLogfile(error.message);                                             % Write error message to log file.
    %end
end

