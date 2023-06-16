function generateStartingParameters(app, exportApp)
% Resets the Simulation status
            
            
            [carDatafile, path] = uigetfile('*.mat');  %open a mat file
            
            % Gets the current value of the track and discipline
            % DropDown
            startingParameters.disciplineID = str2double(app.DisciplineDropDown.Value);    
            startingParameters.carDatafile = carDatafile;
            startingParameters.path = path;
            startingParameters.TrackFileName = app.TrackFileName;
            startingParameters.brakeFunction = app.UseoldbrakesystemCheckBox.Value;
            startingParameters.logCellData = app.LogCellDataCheckBox.Value;
            startingParameters.debugMessages = app.EnableDebugMessagestimeintensivCheckBox.Value;
            startingParameters.startingSpeed = app.StartingSpeed;
            startingParameters.numOfLaps = app.numOfLaps;
            startingParameters.Debug = 0;
            
            % Do not export app when startingParameters need to be saved
            if exportApp
                startingParameters.processDataButtonHandle = app.StartSimulationButton;
                startingParameters.textAreaHandle = app.TextArea;
            else
                startingParameters.processDataButtonHandle = 0;
                startingParameters.textAreaHandle = 0;
            end
            
            % Check if Sensitivity Analysis is selected
            if (app.sensitivityCheckBox.Value)
                startingParameters.sensitivityID = app.VariableDropDown.Value;
                startingParameters.minValue = app.MinEditField.Value;
                startingParameters.stepSize = app.StepsizeEditField.Value;
                startingParameters.numSteps = app.NumStepsEditField.Value;
                
                if (app.secondVariableCheckBox.Value)
                    startingParameters.sensitivityID2 = app.VariableDropDown_2.Value;
                    startingParameters.minValue2 = app.MinEditField_2.Value;
                    startingParameters.stepSize2 = app.StepsizeEditField_2.Value;
                else
                    startingParameters.sensitivityID2 = '0';
                    startingParameters.minValue2 = 0;
                    startingParameters.stepSize2 = 0;
                end
            else
                startingParameters.sensitivityID = '0';
                startingParameters.minValue = 0;
                startingParameters.stepSize = 0;
                startingParameters.numSteps = 1;
                startingParameters.sensitivityID2 = '0';
                startingParameters.minValue2 = 0;
                startingParameters.stepSize2 = 0;
            end
            
            simulationManager(startingParameters);   

end