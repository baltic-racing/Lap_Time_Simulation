function updateColorScheme(app)
    %% Updates the Color Scheme of the Apps (Dark and White Mode)
    % Used for Dark Mode Design:
    % https://material.io/design/color/dark-theme.html#ui-application   
    %
    % By Eric Dornieden, Baltic Racing
    % Copyright (C) 2021, Baltic Racing, all rights reserved.
    %
    % Checks if the setting "ColorScheme" in the app is set to 1 or if the
    % DarkMode Checkbox is checked -> Dark Mode, else White Mode
    if app.settings.ColorScheme || app.DarkModeCheckBox.Value
        %% Dark Mode Colors
        % colors are set to Darkmode colors
        % Primary #BB86FC -> on Primary #000000
        app.PrimaryColor = app.ColorsDark(1);
        app.OnPrimaryColor = app.ColorsDark(2);

        % Primary Variant #3700B3 -> on Primary #000000
        app.PrimaryVariantColor = app.ColorsDark(3);
        app.OnPrimaryVariantColor = app.ColorsDark(4);

        % Secondary 
        app.SecondaryColor = app.ColorsDark(5);
        app.OnSecondaryColor = app.ColorsDark(6);

        % Background
        app.BackgroundColor = app.ColorsDark(7);
        app.OnBackgroundColor = app.ColorsDark(8);

        % Surface 
        app.SurfaceColor = app.ColorsDark(9);
        app.OnSurfaceColor = app.ColorsDark(10);

        % Error
        app.ErrorColor = app.ColorsDark(11);
        app.OnErrorColor = app.ColorsDark(12);

        app.BorderType = app.ColorsDark(13);

        app.Host_Logo = app.ColorsDark(14);
    else
        %% White Mode Colors
        % Colors are set to Whitemode colors
        % Primary #BB86FC -> on Primary #000000
        app.PrimaryColor = app.ColorsBright(1);
        app.OnPrimaryColor = app.ColorsBright(2);

        % Primary Variant #3700B3 -> on Primary #000000
        app.PrimaryVariantColor = app.ColorsBright(3);
        app.OnPrimaryVariantColor = app.ColorsBright(4);

        % Secondary 
        app.SecondaryColor = app.ColorsBright(5);
        app.OnSecondaryColor = app.ColorsBright(6);

        % Background
        app.BackgroundColor = app.ColorsBright(7);
        app.OnBackgroundColor = app.ColorsBright(8);

        % Surface 
        app.SurfaceColor = app.ColorsBright(9);
        app.OnSurfaceColor = app.ColorsBright(10);

        % Error
        app.ErrorColor = app.ColorsBright(11);
        app.OnErrorColor = app.ColorsBright(12);

        app.BorderType = app.ColorsBright(13);             

        app.Host_Logo = app.ColorsBright(14);
    end

    % Update all elements and the other apps if opened (buttons,labels,...)
    
    %% Suspension App
    if (~app.SuspensionButton.Enable)
        app.SuspensionApp.changeColorScheme();
    end

    %% Environment App
    if (~app.EnvironmentButton.Enable)
        app.EnvironmentApp.changeColorScheme();
    end

    %% Aerodynamik App
    if ~app.AerodynamikButton.Enable
        app.AerodynamikApp.changeColorScheme();
    end

    %% Chassis App
    if ~app.ChassisButton.Enable
        app.ChassisApp.changeColorScheme();
    end
    
    %% Customization App
    if ~app.TrackManagerButton.Enable
        app.TrackManagerApp.changeColorScheme();
    end
    
    %% Accumulator App
    if ~app.AccumulatorButton.Enable
        app.AccumulatorApp.changeColorScheme();
    end
    
    %% Drivetrain App
    if ~app.DrivetrainButton.Enable
        app.DrivetrainApp.changeColorScheme();
    end

    %% Track Manager App

    %%% Laptime Simulation GUI
    if ~app.StyleSettingsButton.Enable
        app.CustomizationApp.changeColorScheme();
    end
    app.LaptimeSimulationUIFigure.Color = app.BackgroundColor;           

    %% How-To Page
    app.HowtoTab.BackgroundColor = app.BackgroundColor;  
    %app.HowtoTab.Children

    app.TextArea.BackgroundColor = app.BackgroundColor;
    app.TextArea.FontColor = app.OnSurfaceColor;

    app.HowToLaptimeSimulationPanel.BackgroundColor = app.SurfaceColor;
    app.HowToLaptimeSimulationPanel.ForegroundColor = app.OnSurfaceColor;
    app.HowToLaptimeSimulationPanel.BorderType = app.BorderType;

    app.EditthecarsetupintheCarSetuptab2SaveCarSetupLabel.FontColor = app.OnBackgroundColor;
    app.HowtoTab.BackgroundColor = app.BackgroundColor;
    app.DarkModeCheckBox.FontColor = app.OnSurfaceColor;          
    app.Copyright2021BalticRacingbyEricDorniedenLabel.FontColor = app.OnBackgroundColor;

    app.PatchNotesPanel.BackgroundColor = app.SurfaceColor;
    app.PatchNotesPanel.ForegroundColor = app.OnSurfaceColor;
    app.PatchNotesPanel.BorderType = app.BorderType;

    app.DailyUpdateFeatureRequestsButton.BackgroundColor = app.PrimaryVariantColor;
    app.DailyUpdateFeatureRequestsButton.FontColor = app.OnPrimaryColor;

    app.LaptimeSimulationGitHubButton.BackgroundColor = app.PrimaryVariantColor;
    app.LaptimeSimulationGitHubButton.FontColor = app.OnPrimaryColor; 

    app.WikiBalticRacingButton.BackgroundColor = app.PrimaryVariantColor;
    app.WikiBalticRacingButton.FontColor = app.OnPrimaryColor;

    app.ContactmeButton.BackgroundColor = app.PrimaryVariantColor;
    app.ContactmeButton.FontColor = app.OnPrimaryColor;        
    
    app.StyleSettingsButton.BackgroundColor = app.PrimaryVariantColor;
    app.StyleSettingsButton.FontColor = app.OnPrimaryColor;  
    
    app.FSGRulesButton.BackgroundColor = app.PrimaryVariantColor;
    app.FSGRulesButton.FontColor = app.OnPrimaryColor; 

    app.Image.ImageSource = app.Host_Logo;

    %% Car Setup Page
    app.CarSetupTab.BackgroundColor = app.BackgroundColor; 

    app.EditSetupPanel.BackgroundColor = app.SurfaceColor;
    app.EditSetupPanel.BorderType = app.BorderType;
    app.EditSetupPanel.ForegroundColor = app.OnBackgroundColor;

    app.Copyright2021BalticRacingbyEricDorniedenLabel_2.FontColor = app.OnBackgroundColor;

    app.CarParametersPanel.BackgroundColor = app.SurfaceColor;
    app.CarParametersPanel.ForegroundColor = app.OnSurfaceColor;

    app.SetupPresetsPanel.BackgroundColor = app.SurfaceColor;
    app.SetupPresetsPanel.ForegroundColor = app.OnSurfaceColor;

    app.EnviromentPanel.BackgroundColor = app.SurfaceColor;
    app.EnviromentPanel.ForegroundColor = app.OnSurfaceColor;

    app.TY19Panel.BackgroundColor = app.SurfaceColor;
    app.TY19Panel.ForegroundColor = app.OnSurfaceColor;

    app.TY20Panel.BackgroundColor = app.SurfaceColor;
    app.TY20Panel.ForegroundColor = app.OnSurfaceColor;

    app.TY22Panel.BackgroundColor = app.SurfaceColor;
    app.TY22Panel.ForegroundColor = app.OnSurfaceColor;

    app.LoadButton.BackgroundColor = app.PrimaryVariantColor;
    app.LoadButton.FontColor = app.OnPrimaryColor; 

    app.SaveButton.BackgroundColor = app.PrimaryVariantColor;
    app.SaveButton.FontColor = app.OnPrimaryColor; 

    app.VariableListButton.BackgroundColor = app.PrimaryVariantColor;
    app.VariableListButton.FontColor = app.OnPrimaryColor;

    app.SuspensionButton.BackgroundColor = app.PrimaryVariantColor;
    app.SuspensionButton.FontColor = app.OnPrimaryColor;

    app.ChassisButton.BackgroundColor = app.PrimaryVariantColor;
    app.ChassisButton.FontColor = app.OnPrimaryColor;

    app.DrivetrainButton.BackgroundColor = app.PrimaryVariantColor;
    app.DrivetrainButton.FontColor = app.OnPrimaryColor;

    app.AerodynamikButton.BackgroundColor = app.PrimaryVariantColor;
    app.AerodynamikButton.FontColor = app.OnPrimaryColor;

    app.AccumulatorButton.BackgroundColor = app.PrimaryVariantColor;
    app.AccumulatorButton.FontColor = app.OnPrimaryColor;

    app.EnvironmentButton.BackgroundColor = app.PrimaryVariantColor;
    app.EnvironmentButton.FontColor = app.OnPrimaryColor;

    app.NoDownforceButton.BackgroundColor = app.PrimaryVariantColor;
    app.NoDownforceButton.FontColor = app.OnPrimaryColor;
    app.NoDownforceButton_2.BackgroundColor = app.PrimaryVariantColor;
    app.NoDownforceButton_2.FontColor = app.OnPrimaryColor;
    app.NoDownforceButton_3.BackgroundColor = app.PrimaryVariantColor;
    app.NoDownforceButton_3.FontColor = app.OnPrimaryColor;
    app.LowDownforceButton.BackgroundColor = app.PrimaryVariantColor;
    app.LowDownforceButton.FontColor = app.OnPrimaryColor;
    app.LowDownforceButton_2.BackgroundColor = app.PrimaryVariantColor;
    app.LowDownforceButton_2.FontColor = app.OnPrimaryColor;
    app.LowDownforceButton_3.BackgroundColor = app.PrimaryVariantColor;
    app.LowDownforceButton_3.FontColor = app.OnPrimaryColor;
    app.MedDownforceButton.BackgroundColor = app.PrimaryVariantColor;
    app.MedDownforceButton.FontColor = app.OnPrimaryColor;
    app.MedDownforceButton_2.BackgroundColor = app.PrimaryVariantColor;
    app.MedDownforceButton_2.FontColor = app.OnPrimaryColor;
    app.MedDownforceButton_3.BackgroundColor = app.PrimaryVariantColor;
    app.MedDownforceButton_3.FontColor = app.OnPrimaryColor;
    app.HighDownforceButton.BackgroundColor = app.PrimaryVariantColor;
    app.HighDownforceButton.FontColor = app.OnPrimaryColor;
    app.HighDownforceButton_2.BackgroundColor = app.PrimaryVariantColor;
    app.HighDownforceButton_2.FontColor = app.OnPrimaryColor;
    app.HighDownforceButton_3.BackgroundColor = app.PrimaryVariantColor;
    app.HighDownforceButton_3.FontColor = app.OnPrimaryColor;

    app.Image_2.ImageSource = app.Host_Logo;

    %% Simulation Setup Tab
    app.SimulationSetupTab.BackgroundColor = app.BackgroundColor;
    app.SimulationTextArea.BackgroundColor = app.BackgroundColor;

    app.Copyright2021BalticRacingbyEricDorniedenLabel_3.FontColor = app.OnBackgroundColor;

    app.ScenarioPanel.BackgroundColor = app.SurfaceColor;
    app.ScenarioPanel.ForegroundColor = app.OnSurfaceColor;
    app.ScenarioPanel.BorderType = app.BorderType;

    app.SimulationStatusPanel.BackgroundColor = app.SurfaceColor;
    app.SimulationStatusPanel.ForegroundColor = app.OnSurfaceColor;
    app.SimulationStatusPanel.BorderType = app.BorderType;

    app.SensitivityAnalysisPanel.BackgroundColor = app.SurfaceColor;
    app.SensitivityAnalysisPanel.ForegroundColor = app.OnSurfaceColor;
    app.SensitivityAnalysisPanel.BorderType = app.BorderType;

    app.TrackDropDownLabel.FontColor = app.OnPrimaryColor;
    app.TrackDropDown.BackgroundColor = app.PrimaryColor;
    app.TrackDropDown.FontColor = app.OnPrimaryColor;

    app.DisciplineDropDownLabel.FontColor = app.OnPrimaryColor;
    app.DisciplineDropDown.BackgroundColor = app.PrimaryColor;
    app.DisciplineDropDown.FontColor = app.OnPrimaryColor;

    app.sensitivityCheckBox.FontColor = app.OnPrimaryColor;

    app.VariableDropDownLabel.FontColor = app.OnPrimaryColor;
    app.VariableDropDown.BackgroundColor = app.SurfaceColor; 
    app.VariableDropDown.FontColor = app.OnSurfaceColor;

    app.secondVariableCheckBox.FontColor = app.OnPrimaryColor;
    app.VariableDropDown_2.BackgroundColor = app.SurfaceColor; 
    app.VariableDropDown_2.FontColor = app.OnSurfaceColor;

    app.MinEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.MinEditField.BackgroundColor = app.SurfaceColor; 
    app.MinEditField.FontColor = app.OnSurfaceColor;
    app.MinEditField_2.BackgroundColor = app.SurfaceColor; 
    app.MinEditField_2.FontColor = app.OnSurfaceColor;

    app.StepsizeEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.StepsizeEditField.BackgroundColor = app.SurfaceColor;
    app.StepsizeEditField.FontColor = app.OnSurfaceColor;
    app.StepsizeEditField_2.BackgroundColor = app.SurfaceColor;
    app.StepsizeEditField_2.FontColor = app.OnSurfaceColor;

    app.NumStepsEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.NumStepsEditField.BackgroundColor = app.SurfaceColor;
    app.NumStepsEditField.FontColor = app.OnSurfaceColor;

    app.TrackManagerButton.BackgroundColor = app.PrimaryVariantColor;
    app.TrackManagerButton.FontColor = app.OnPrimaryColor;

    app.UpdateTracksButton.BackgroundColor = app.PrimaryVariantColor;
    app.UpdateTracksButton.FontColor = app.OnPrimaryColor;

    app.ClearButton.BackgroundColor = app.PrimaryVariantColor;
    app.ClearButton.FontColor = app.OnPrimaryColor;

    app.LogCellDataCheckBox.FontColor = app.OnPrimaryColor;
    
    app.TrackDataOptimizationSpinnerLabel.FontColor = app.OnPrimaryColor;
    app.TrackDataOptimizationSpinner.FontColor = app.OnPrimaryColor;
    app.TrackDataOptimizationSpinner.BackgroundColor = app.PrimaryColor;

    app.UIAxes2.Title.Color = app.OnPrimaryColor;
    app.UIAxes2.XColor = app.OnPrimaryColor;
    app.UIAxes2.YColor = app.OnPrimaryColor;
    app.UIAxes2.Color = app.SurfaceColor;
    app.LengthLabel.FontColor = app.OnPrimaryColor;
    
    app.UseoldbrakesystemCheckBox.FontColor = app.OnPrimaryColor;

    %% Simulation Results
    app.SimulationResultsTab.BackgroundColor = app.BackgroundColor;

    app.Results1RunPanel.BackgroundColor = app.SurfaceColor;
    app.Results1RunPanel.ForegroundColor = app.OnSurfaceColor;
    app.Results1RunPanel.BorderType = app.BorderType;

    app.Results2RunPanel.BackgroundColor = app.SurfaceColor;
    app.Results2RunPanel.ForegroundColor = app.OnSurfaceColor;
    app.Results2RunPanel.BorderType = app.BorderType;

    app.DeltaSecondRunFirstRunPanel.BackgroundColor = app.SurfaceColor;
    app.DeltaSecondRunFirstRunPanel.ForegroundColor = app.OnSurfaceColor;
    app.DeltaSecondRunFirstRunPanel.BorderType = app.BorderType;

    app.PlotResultsPanel.BackgroundColor = app.SurfaceColor;
    app.PlotResultsPanel.ForegroundColor = app.OnSurfaceColor;
    app.PlotResultsPanel.BorderType = app.BorderType;

    app.SensitivityAnalysisPanel_2.BackgroundColor = app.SurfaceColor;
    app.SensitivityAnalysisPanel_2.ForegroundColor = app.OnSurfaceColor;
    app.SensitivityAnalysisPanel_2.BorderType = app.BorderType;

    app.CompareRunsPanel.BackgroundColor = app.SurfaceColor;
    app.CompareRunsPanel.ForegroundColor = app.OnSurfaceColor;
    app.CompareRunsPanel.BorderType = app.BorderType;

    app.DataInspectorPanel.BackgroundColor = app.SurfaceColor;
    app.DataInspectorPanel.ForegroundColor = app.OnSurfaceColor;
    app.DataInspectorPanel.BorderType = app.BorderType;

    app.Panel.BackgroundColor = app.SurfaceColor;
    app.Panel.ForegroundColor = app.OnSurfaceColor;
    %app.Panel.BorderType = BorderType;

    app.Load1RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.Load1RunButton.FontColor = app.OnPrimaryColor;
    app.Load2RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.Load2RunButton.FontColor = app.OnPrimaryColor;

    app.TelemetryconverterButton.BackgroundColor = app.PrimaryVariantColor;
    app.TelemetryconverterButton.FontColor = app.OnPrimaryColor;
    app.OpenTrackAnalyzerButton.BackgroundColor = app.PrimaryVariantColor;
    app.OpenTrackAnalyzerButton.FontColor = app.OnPrimaryColor;

    app.DrawPlot1RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.DrawPlot1RunButton.FontColor = app.OnPrimaryColor;
    app.DrawPlot2RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.DrawPlot2RunButton.FontColor = app.OnPrimaryColor;
    app.DrawSensitivityPlotButton.BackgroundColor = app.PrimaryVariantColor;
    app.DrawSensitivityPlotButton.FontColor = app.OnPrimaryColor;
    app.CompareRunPlotsButton.BackgroundColor = app.PrimaryVariantColor;
    app.CompareRunPlotsButton.FontColor = app.OnPrimaryColor;

    app.t_lapLabel.FontColor = app.OnBackgroundColor;
    app.t_lapLabel_2.FontColor = app.OnBackgroundColor;

    app.LapDistanceLabel.FontColor = app.OnBackgroundColor;
    app.LapDistanceLabel_2.FontColor = app.OnBackgroundColor;

    app.EndurancetimeLabel.FontColor = app.OnBackgroundColor;
    app.EndurancetimeLabel_2.FontColor = app.OnBackgroundColor;

    app.EnduranceDistanceLabel.FontColor = app.OnBackgroundColor;
    app.EnduranceDistanceLabel_2.FontColor = app.OnBackgroundColor;

    app.EnergyconsumptionLabel.FontColor = app.OnBackgroundColor;
    app.EnergyconsumptionLabel_2.FontColor = app.OnBackgroundColor;

    app.withrecuperationLabel.FontColor = app.OnBackgroundColor;
    app.withrecuperationLabel_2.FontColor = app.OnBackgroundColor;

    app.PowerLimitLabel.FontColor = app.OnBackgroundColor;
    app.PowerLimitLabel_2.FontColor = app.OnBackgroundColor;

    app.ComputationtimeLabel.FontColor = app.OnBackgroundColor;
    app.ComputationtimeLabel_2.FontColor = app.OnBackgroundColor;
    
    app.SkidPadTimeLabel.FontColor = app.OnBackgroundColor;
    app.SkidPadTimeLabel_2.FontColor = app.OnBackgroundColor;
    
    app.SkidPadSpeedLabel.FontColor = app.OnBackgroundColor;
    app.SkidPadSpeedLabel_2.FontColor = app.OnBackgroundColor;


    if (app.LapDistanceLabel_3.Text == "Lap Distance =")
        app.EnergyconsumptionWOrecuperationLabel_3.FontColor = app.OnBackgroundColor;
        app.WrecuperationLabel_3.FontColor = app.OnBackgroundColor;
        app.t_lapLabel_3.FontColor = app.OnBackgroundColor;
        app.LapDistanceLabel_3.FontColor = app.OnBackgroundColor;
        app.EndurancetimeLabel_3.FontColor = app.OnBackgroundColor;
    end

    app.CreateData1RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.CreateData1RunButton.FontColor = app.OnPrimaryColor;
    app.CreateData2RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.CreateData2RunButton.FontColor = app.OnPrimaryColor;

    app.LoadSetup1RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.LoadSetup1RunButton.FontColor = app.OnPrimaryColor;
    app.LoadSetup2RunButton.BackgroundColor = app.PrimaryVariantColor;
    app.LoadSetup2RunButton.FontColor = app.OnPrimaryColor;

    app.OpenDataInspectorButton.BackgroundColor = app.PrimaryVariantColor;
    app.OpenDataInspectorButton.FontColor = app.OnPrimaryColor;

    app.SelectLapRunnumberEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.SelectLapRunnumberEditField.BackgroundColor = app.SurfaceColor;
    app.SelectLapRunnumberEditField.FontColor = app.OnSurfaceColor;

    app.RunnumberEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.RunnumberEditField.BackgroundColor = app.SurfaceColor;
    app.RunnumberEditField.FontColor = app.OnSurfaceColor;

    app.ResultfilenameEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.ResultfilenameEditField.BackgroundColor = app.SurfaceColor;
    app.ResultfilenameEditField.FontColor = app.OnSurfaceColor;

    app.SetupviewerTextAreaLabel.FontColor = app.OnPrimaryColor;
    app.SetupviewerTextArea.BackgroundColor = app.SurfaceColor;
    app.SetupviewerTextArea.FontColor = app.OnSurfaceColor;

    app.DeletecsvfileCheckBox.FontColor = app.OnPrimaryColor;
    app.ExportasdistanceCheckBox.FontColor = app.OnPrimaryColor;

    app.SelectPlotDropDownLabel.FontColor = app.OnPrimaryColor;
    app.SelectPlotDropDown.BackgroundColor = app.SurfaceColor;
    app.SelectPlotDropDown.FontColor = app.OnSurfaceColor;

    app.xaxisDropDown_2Label.FontColor = app.OnPrimaryColor;
    app.xaxisDropDown_2.BackgroundColor = app.SurfaceColor;
    app.xaxisDropDown_2.FontColor = app.OnSurfaceColor;

    app.yaxisDropDown_2Label.FontColor = app.OnPrimaryColor;
    app.yaxisDropDown_2.BackgroundColor = app.SurfaceColor;
    app.yaxisDropDown_2.FontColor = app.OnSurfaceColor;

    app.zaxisDropDown_2Label.FontColor = app.OnPrimaryColor;
    app.zaxisDropDown_2.BackgroundColor = app.SurfaceColor;
    app.zaxisDropDown_2.FontColor = app.OnSurfaceColor;

    app.ResolutionxEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.ResolutionxEditField.BackgroundColor = app.SurfaceColor;
    app.ResolutionxEditField.FontColor = app.OnSurfaceColor;

    app.yEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.ResolutionyEditField.BackgroundColor = app.SurfaceColor;
    app.ResolutionyEditField.FontColor = app.OnSurfaceColor;

    app.DrawDotsCheckBox.FontColor = app.OnPrimaryColor;

    app.SelectPlotDropDown_3Label.FontColor = app.OnPrimaryColor;
    app.SelectPlotDropDown_2.BackgroundColor = app.SurfaceColor;
    app.SelectPlotDropDown_2.FontColor = app.OnSurfaceColor;

    %% Simulation Replay
    app.SimulationReplayWIPTab.BackgroundColor = app.BackgroundColor;

    app.Copyright2021BalticRacingbyEricDorniedenLabel_4.FontColor = app.OnBackgroundColor;

    app.Panel_2.BackgroundColor = app.SurfaceColor;
    app.Panel_2.ForegroundColor = app.OnSurfaceColor;
    app.Panel_2.BorderType = app.BorderType;
    
    % Variable Dropdowns
    app.DropDown.BackgroundColor = app.PrimaryColor;
    %app.DropDown.FontColor = app.OnSurfaceColor;
    app.DropDown_2.BackgroundColor = app.PrimaryColor;
    %app.DropDown_2.FontColor = app.OnSurfaceColor;
    app.DropDown_3.BackgroundColor = app.PrimaryColor;
    %app.DropDown_3.FontColor = app.OnSurfaceColor;
    app.DropDown_4.BackgroundColor = app.PrimaryColor;
    %app.DropDown_4.FontColor = app.OnSurfaceColor;
    app.DropDown_5.BackgroundColor = app.PrimaryColor;
    %app.DropDown_5.FontColor = app.OnSurfaceColor;
    app.DropDown_6.BackgroundColor = app.PrimaryColor;
    %app.DropDown_6.FontColor  = app.OnSurfaceColor;
    app.DropDown_8.BackgroundColor = app.PrimaryColor;
    app.DropDown_8.FontColor = app.OnSurfaceColor;
    app.DropDown_9.BackgroundColor = app.PrimaryColor;
    app.DropDown_9.FontColor = app.OnSurfaceColor;
    app.DropDown_10.BackgroundColor = app.PrimaryColor;
    app.DropDown_10.FontColor = app.OnSurfaceColor;
    app.DropDown_11.BackgroundColor = app.PrimaryColor;
    app.DropDown_11.FontColor = app.OnSurfaceColor;
    app.DropDown_12.BackgroundColor = app.PrimaryColor;
    app.DropDown_12.FontColor = app.OnSurfaceColor;
    app.DropDown_13.BackgroundColor = app.PrimaryColor;
    app.DropDown_13.FontColor = app.OnSurfaceColor;
    app.DropDown_14.BackgroundColor = app.PrimaryColor;
    app.DropDown_14.FontColor = app.OnSurfaceColor;


    app.LoadReplayButton.BackgroundColor = app.PrimaryVariantColor;
    app.LoadReplayButton.FontColor = app.OnPrimaryColor;
    
    app.PlayButton.BackgroundColor = app.PrimaryVariantColor;
    app.PlayButton.FontColor = app.OnPrimaryColor;
    
    app.StopButton.BackgroundColor = app.PrimaryVariantColor;
    app.StopButton.FontColor = app.OnPrimaryColor;
    
    app.BrakeLampLabel.FontColor = app.OnPrimaryColor;
    app.DRSStatusLampLabel.FontColor = app.OnPrimaryColor;
    
    app.RunNumberSpinnerLabel.FontColor = app.OnPrimaryColor;
    app.RunNumberSpinner.BackgroundColor = app.SurfaceColor;
    app.RunNumberSpinner.FontColor = app.OnSurfaceColor;

    app.UIAxes3.Title.Color = app.OnPrimaryColor;
    app.UIAxes3.XColor = app.OnPrimaryColor;
    app.UIAxes3.YColor = app.OnPrimaryColor;
    app.UIAxes3.Color = app.SurfaceColor;
    
    app.UIAxes4.Title.Color = app.OnPrimaryColor;
    app.UIAxes4.XColor = app.OnPrimaryColor;
    app.UIAxes4.YColor = app.OnPrimaryColor;
    app.UIAxes4.Color = app.SurfaceColor;

%     app.var1Label.FontColor = app.OnPrimaryColor;
%     app.var2Label.FontColor = app.OnPrimaryColor;
%     app.var3Label.FontColor = app.OnPrimaryColor;
%     app.var4Label.FontColor = app.OnPrimaryColor;
%     app.var5Label.FontColor = app.OnPrimaryColor;
%     app.var6Label.FontColor = app.OnPrimaryColor;
    app.var7Label.FontColor = app.OnPrimaryColor;
    app.var8Label.FontColor = app.OnPrimaryColor;
    app.var9Label.FontColor = app.OnPrimaryColor;
    app.var10Label.FontColor = app.OnPrimaryColor;
    app.var11Label.FontColor = app.OnPrimaryColor;
    app.var12Label.FontColor = app.OnPrimaryColor;
    
    app.SpeedLabel.FontColor = app.OnPrimaryColor;
    app.GearLabel.FontColor = app.OnPrimaryColor;
    app.RPMLabel.FontColor = app.OnPrimaryColor;
    app.RPMLabel_2.FontColor = app.OnPrimaryColor;
    app.TimeLabel.FontColor = app.OnPrimaryColor;
    app.DistanceLabel.FontColor = app.OnPrimaryColor;
    app.RadiusLabel.FontColor = app.OnPrimaryColor;
    
    app.Slider.FontColor = app.OnPrimaryColor;
    app.Slider_2.FontColor = app.OnPrimaryColor;

    app.LowerLimitEditFieldLabel.FontColor = app.OnPrimaryColor;
    app.LowerLimitEditField.BackgroundColor = app.SurfaceColor;
    app.LowerLimitEditField.FontColor = app.OnSurfaceColor;
    
    app.UpperLimitEditField_2Label.FontColor = app.OnPrimaryColor;
    app.UpperLimitEditField.BackgroundColor = app.SurfaceColor;
    app.UpperLimitEditField.FontColor = app.OnSurfaceColor;

    app.AutoCheckBox.FontColor = app.OnSurfaceColor;

    app.DrawApexesCheckBox.FontColor = app.OnSurfaceColor;

    app.wheelloadLabel.FontColor = app.OnSurfaceColor;
    
    app.RPMGauge.FontColor = app.OnPrimaryColor;
    app.RPMGauge.BackgroundColor = app.BackgroundColor;

    app.ExportasCSVButton.BackgroundColor = app.PrimaryVariantColor;
    app.ExportasCSVButton.FontColor = app.OnPrimaryColor;
    
    app.Copyright2021BalticRacingbyEricDorniedenLabel_5.FontColor = app.OnSurfaceColor;
end