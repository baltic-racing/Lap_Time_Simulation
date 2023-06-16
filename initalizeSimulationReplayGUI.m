%% initalizeSimulationReplay.m
% Initalizes the Simulation Replay GUI / Loads save file and sets inital
% conditions for the replay.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [result] = initalizeSimulationReplayGUI(app)
    %try
        [file_run, path] = uigetfile('*.mat');                                    % Load the file name of the selected save file.

        result = load(path+"/"+file_run, '-mat');                               % open the .mat file with the given name.

        % Get the existing data from the uitable
        existingData = get(app.UITable_2, 'Data');
        
        [~, name, ~] = fileparts(file_run);

        % Create a new entry as a cell array or a row vector
        newEntry = {name, 11, 12, 14};
        
        % Append the new entry to the existing data and create a new table
        updatedData = [existingData; newEntry];

        % Update table
        set(app.UITable_2, 'Data', updatedData)
        
        % Get the number of rows in the table
        numRows = size(app.UITable_2.Data, 1);  
        
     
        %% set slider limits
        % ToDo add limit
        runNumber = app.RunNumberSpinner.Value;
        
        [~, columns] = size(result(runNumber).Track);
        %% Add track variables as single parameters
        if (columns == 4)
            result.x_Track = result(runNumber).Track(:,1);           % [m] X-Koordinate der Strecke
            result.y_Track = result(runNumber).Track(:,2);           % [m] Y-Koordinate der Strecke
            result.z_Track = [];                                      % [m] Z-Koordinate der Strecke
            result.s = result(runNumber).Track(:,3);                 % [m] Verlauf der Streckenlänge
            result.R = result(runNumber).Track(:,4);                 % [m] Kurvenradien
        else
            result.x_Track = result(runNumber).Track(:,1);           % [m] X-Koordinate der Strecke
            result.y_Track = result(runNumber).Track(:,2);           % [m] Y-Koordinate der Strecke
            result.z_Track = result(runNumber).Track(:,3);           % [m] Z-Koordinate der Strecke
            result.s = result(runNumber).Track(:,4);                 % [m] Verlauf der Streckenlänge
            result.R = result(runNumber).Track(:,5);                 % [m] Kurvenradien
        end   

        if ~isempty(app.activeResult)
            return
        end

        % Reset style for all rows 
        removeStyle(app.UITable_2)

        % Set style for selected field
        s = uistyle('BackgroundColor','red');
        addStyle(app.UITable_2,s,'row',size(app.UITable_2.Data, 1))

        app.Slider_2.Limits = [0 result(runNumber).t(end)];                           % Set maximum of the time slider to the overall laptime.
        app.Slider_2.Value = 0;                                                       % Reset Slider Value (Start Replay at the beginning of the lap) 

        %% Update setup variables
        % Get Setup variables
        %setupData = cell(numVariables, 2);
        setupData = {};

        fieldNames = fieldnames(result);
        numVariables = numel(fieldNames);

        for i = 1:numVariables
            fieldName = fieldNames{i};  % Get the current field name
    
            if length(result.(fieldName)) == 1
                % Extract the value of the variable
                variableValue = num2str(result.(fieldName));
                
                % Assign the variable name and value to the data cell array
                setupData{(end+1), 1} = fieldName;
                setupData{(end), 2} = variableValue;
            end
        end

        % Update the table with the data and column names
        set(app.UITable2, 'Data', setupData);
       
end
