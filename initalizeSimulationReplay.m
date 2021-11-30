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

        cla(app.UIAxes3);                                                   % Clear UIAxes when loading a new save file
        cla(app.UIAxes4);
        
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

%         for i=1:length(result(runNumber).Track(:,1))
%             M(i) = recordPedalPlot(app, i, runNumber, result);
%         end

        drawPedalPlots(app)

        if app. DarkModeCheckBox.Value                                               % Check if Dark Mode is active, if so plot line in white else plot line in black.
            plot(app.UIAxes3,result(runNumber).Track(:,1),result(runNumber).Track(:,2),'w')       % Plot Track in white.
        else
            plot(app.UIAxes3,result(runNumber).Track(:,1),result(runNumber).Track(:,2),'k')       % Plot Track in black.
        end

        hold(app.UIAxes3,'on');

        title(app.UIAxes3,'','FontSize',12);                                        % Set the title of the UIAxes to an empty String
        app.UIAxes2.LineWidth = 1.5;                                                % Set the Line width of the plot

        app.Slider.Limits = [0 result(runNumber).t(end)];                           % Set maximum of the time slider to the overall laptime.
        app.Slider.Value = 0;                                                       % Reset Slider Value (Start Replay at the beginning of the lap)       

        runID = runNumber;
        saveFile = result;



        saveFileData = {'A_Accu_cell' saveFile(runID).A_accu_cell(1:end-1); 'DRS_status' saveFile(runID).DRS_status(1:end-1);...
            'aVX' saveFile(runID).aVX(:)./9.81; 'BPPsignal' saveFile(runID).BPPsignal(:); 'E_Accu' saveFile(runID).E_Accu(1:end-1);...
            'E_Akku_Recu' saveFile(runID).E_Accu_Recu(1:end-1); 'E_heat' saveFile(runID).E_heat(1:end-1); 'E_res' saveFile(runID).E_res(1:end-1);...
            'vV' saveFile(runID).vV(1:end-1).*3.6; 'psi1' saveFile(runID).psi1(1:end-1); 'gear' saveFile(runID).gearSelection(1:end);...
            'alpha_f' saveFile(runID).alpha_f(1:end-1); 'alpha_r' saveFile(runID).alpha_r(1:end-1); 'aRev' saveFile(runID).aRev(:);...
            'aVY' saveFile(runID).aVY(:)./9.81; 'beta' saveFile(runID).beta(1:end-1).*(180/pi); 'FL' saveFile(runID).FL(:); 'FR' saveFile(runID).FR(:);...
            'vVYmax' saveFile(runID).vVYmax(1:end-1).*3.6; 'FVX' saveFile(runID).FVX(:); 'FVX_rl' saveFile(runID).FVX_rl(:); 'FVX_rr' saveFile(runID).FVX_rr(:);...
            'FVY' saveFile(runID).FVY(:); 'FWXmax_r' saveFile(runID).FWXmax_r(1:end-1); 'FVX_fl' saveFile(runID).FVX_fl(:); 'FVX_fr' saveFile(runID).FVX_fr(:);...
            'FWXmax_rl' saveFile(runID).FWXmax_rl(1:end-1); 'FWXmax_rr' saveFile(runID).FWXmax_rr(1:end-1); 'FWXmax_f' saveFile(runID).FWXmax_f(1:end-1);...
            'FWXmax_fl' saveFile(runID).FWXmax_fl(1:end-1); 'FWXmax_fr' saveFile(runID).FWXmax_fr(1:end-1); 'FWYf' saveFile(runID).FWYf(1:end);...
            'FWYr' saveFile(runID).FWYr(1:end); 'FWYmax_r' saveFile(runID).FWYmax_r(1:end-1);...
            'FWYmax_rl' saveFile(runID).FWYmax_rl(1:end-1); 'FWYmax_rr' saveFile(runID).FWYmax_rr(1:end-1); 'FWYmax_f' saveFile(runID).FWYmax_f(1:end-1);...
            'FWYmax_fl' saveFile(runID).FWYmax_fl(1:end-1); 'FWYmax_fr' saveFile(runID).FWYmax_fr(1:end-1); 'FWZr' saveFile(runID).FWZr(1:end-1);...
            'FWZf' saveFile(runID).FWZf(1:end-1); 'Faero' saveFile(runID).Faero(:); 'Fdrag' saveFile(runID).Fdr(:); 'Mi' saveFile(runID).Mi(:);...
            'P_M' saveFile(runID).P_M(:); 'P_el' saveFile(runID).P_el(:); 'TC_front' saveFile(runID).TC_front(1:end-1);...
            'Rdyn_rl' saveFile(runID).Rdyn_rl(1:end-1); 'Rdyn_rr' saveFile(runID).Rdyn_rr(1:end-1);...
            'Rdyn_fl' saveFile(runID).Rdyn_fl(1:end-1); 'Rdyn_fr' saveFile(runID).Rdyn_fr(1:end-1); 'slipY_r' saveFile(runID).slipY_r(1:end-1);...
            'slipY_f' saveFile(runID).slipY_f(1:end-1); 'TC' saveFile(runID).TC(1:end-1); 'Tirelimit' saveFile(runID).Tirelimit(1:end-1);...
            'cZ_rl' saveFile(runID).cZ_rl(1:end-1); 'cZ_rr' saveFile(runID).cZ_rr(1:end-1);...
            'cZ_fl' saveFile(runID).cZ_fl(1:end-1); 'cZ_fr' saveFile(runID).cZ_fr(1:end-1);...
            'dFWZrl_aero' saveFile(runID).dFWZrl_aero(:); 'dFWZrl_x' saveFile(runID).dFWZrl_x(:); 'dFWZrl_y' saveFile(runID).dFWZrl_y(:);...
            'dFWZrr_aero' saveFile(runID).dFWZrr_aero(:); 'dFWZrr_x' saveFile(runID).dFWZrr_x(:); 'dFWZrr_y' saveFile(runID).dFWZrr_y(:);...
            'dFWZfl_aero' saveFile(runID).dFWZfl_aero(:); 'dFWZfl_x' saveFile(runID).dFWZfl_x(:); 'dFWZfl_y' saveFile(runID).dFWZfl_y(:);...
            'dFWZfr_aero' saveFile(runID).dFWZfr_aero(:); 'dFWZfr_x' saveFile(runID).dFWZfr_x(:); 'dFWZfr_y' saveFile(runID).dFWZfr_y(:);...
            'kappa_fl' saveFile(runID).kappa_fl(1:end-1); 'kappa_fr' saveFile(runID).kappa_fr(1:end-1);...
            'kappa_rl' saveFile(runID).kappa_rl(1:end-1); 'kappa_rr' saveFile(runID).kappa_rr(1:end-1);...
            'l_contact_patch_fl' saveFile(runID).l_contact_patch_fl(1:end-1); 'l_contact_patch_fr' saveFile(runID).l_contact_patch_fr(1:end-1);...
            'l_contact_patch_rl' saveFile(runID).l_contact_patch_rl(1:end-1); 'l_contact_patch_rr' saveFile(runID).l_contact_patch_rr(1:end-1);...
            'motor_eff' saveFile(runID).motor_eff(1:end-1); 'ni' saveFile(runID).ni(:); 'vRev' saveFile(runID).vRev(:).*3.6;...
            'x-Coordinate' saveFile(runID).Track(1:end-1,1); 'y-Coordinate' saveFile(runID).Track(1:end-1,2); 'z-Coordinate' saveFile(runID).Track(1:end-1,3);...
            'Distance' saveFile(runID).Track(1:end-1,4);'Corner Radius' saveFile(runID).Track(1:end-1,5); 'vWoBrake' saveFile(runID).vWoBrake(1:end-1).*3.6;...
            'alpha_rl' saveFile(runID).alpha_rl(1:end-1); 'alpha_rr' saveFile(runID).alpha_rr(1:end-1); 'alpha_fr' saveFile(runID).alpha_fr(1:end-1); 'alpha_fl' saveFile(runID).alpha_fl(1:end-1);...
            'ABS status' saveFile(runID).ABS(1:end-1); 'FWZtot' saveFile(runID).FWZtot(1:end-1); 'M_tractive' saveFile(runID).M_tractive(1:end);...
            'FB' saveFile(runID).FB(1:end); 'FB_fl' saveFile(runID).FB_fl(1:end); 'FB_fr' saveFile(runID).FB_fr(1:end); 'FB_rl' saveFile(runID).FB_rl(1:end); 'FB_rr' saveFile(runID).FB_rr(1:end);...
            'FU_fl' saveFile(runID).FU_fl(1:end); 'FU_fr' saveFile(runID).FU_fr(1:end); 'FU_rl' saveFile(runID).FU_rl(1:end); 'FU_rr' saveFile(runID).FU_rr(1:end);...
            'FWZ_fl' saveFile(runID).FWZ_fl(1:end-1); 'FWZ_fr' saveFile(runID).FWZ_fr(1:end-1); 'FWZ_rl' saveFile(runID).FWZ_rl(1:end-1); 'FWZ_rr' saveFile(runID).FWZ_rr(1:end-1);...
            'P_Bh' saveFile(runID).P_Bh(1:end-1); 'P_Mloss' saveFile(runID).P_Mloss(1:end-1); 'P_tractive' saveFile(runID).P_tractive(1:end); 't' saveFile(runID).t(1:end-1);...
            'delta' saveFile(runID).delta(1:end-1).*(180/pi); 'delta_fl' saveFile(runID).delta_fl(1:end-1).*(180/pi); 'delta_fr' saveFile(runID).delta_fr(1:end-1).*(180/pi);
            'rollMoment_f' saveFile(runID).rollMoment_f(:); 'rollMoment_r' saveFile(runID).rollMoment_r(:); 'ackermann' saveFile(runID).ackermann(1:end-1); 'ackermannPercent' saveFile(runID).ackermannPercent(1:end-1);
            'rollAngleChassis' saveFile(runID).rollAngleChassis(:);};
    
        app.DropDown_14.Items = saveFileData(:,1);
        app.DropDown_14.ItemsData = saveFileData(:,2);

        app.DropDown_13.Items = saveFileData(:,1);
        app.DropDown_13.ItemsData = saveFileData(:,2);

        app.DropDown_12.Items = saveFileData(:,1);
        app.DropDown_12.ItemsData = saveFileData(:,2);

        app.DropDown_11.Items = saveFileData(:,1);
        app.DropDown_11.ItemsData = saveFileData(:,2);

        app.DropDown_10.Items = saveFileData(:,1);
        app.DropDown_10.ItemsData = saveFileData(:,2);

        app.DropDown_9.Items = saveFileData(:,1);
        app.DropDown_9.ItemsData = saveFileData(:,2);

        app.DropDown_8.Items = saveFileData(:,1);
        app.DropDown_8.ItemsData = saveFileData(:,2);

        app.DropDown_6.Items = saveFileData(:,1);
        app.DropDown_6.ItemsData = saveFileData(:,2);

        app.DropDown_5.Items = saveFileData(:,1);
        app.DropDown_5.ItemsData = saveFileData(:,2);

        app.DropDown_4.Items = saveFileData(:,1);
        app.DropDown_4.ItemsData = saveFileData(:,2);

        app.DropDown_3.Items = saveFileData(:,1);
        app.DropDown_3.ItemsData = saveFileData(:,2);

        app.DropDown_2.Items = saveFileData(:,1);
        app.DropDown_2.ItemsData = saveFileData(:,2);

        app.DropDown.Items = saveFileData(:,1);
        app.DropDown.ItemsData = saveFileData(:,2);
        
        %app.SpeedLabel.Text = "Speed: " + num2str(result(runNumber).vV(ceil(app.Slider.Value)));
    %catch error
        %writeToLogfile(error.message);                                             % Write error message to log file.
    %end
    
%         axes(handles.axes1);
%         y = imread('SteeringWheel1.png');
%         imshow(app.UIAxes6)
end

function drawPedalPlots(app)
    X = categorical({'Brake','Throttle'});
    
    brake = 0;
    throttle = 0;
    pedalInputs = [brake throttle];
    
    c = bar(app.UIAxes5, X,[100,100],'FaceColor','flat');
    c(1).CData = [1 1 1; 1 1 1];
    
    hold(app.UIAxes5,'on');
    
    pedalPlot = bar(app.UIAxes5,X, pedalInputs,'FaceColor','flat');
    pedalPlot.YDataSource = 'pedalInputs';
    
    pedalPlot(1).CData = [1 0 0; 0 1 0];
    ylim([0 100])
    ax = gca;
    disableDefaultInteractivity(ax)
    axis off

    xtips2 = pedalPlot(1).XEndPoints;
    ytips2 = [10 10];
    labels2 = string(pedalPlot(1).YData);
    delete(app.pedalLabel);
    app.pedalLabel = text(app.UIAxes5, xtips2,ytips2,labels2,'HorizontalAlignment','center','FontSize',10);
end

function [frame] = recordPedalPlot(app, i, runNumber, result)
    X = categorical({'Brake','Throttle'});
    
    throttle = result(runNumber).P_M/max(result(runNumber).P_M);
    brake = result(runNumber).FB/max(result(runNumber).FB);

    pedalInputs = [round(brake(i)*100,0) round(throttle(i)*100,0)];
    
    c = bar(app.UIAxes5, X,[100,100],'FaceColor','flat');
    c(1).CData = [1 1 1; 1 1 1];
    
    hold(app.UIAxes5,'on');
    
    pedalPlot = bar(app.UIAxes5,X, pedalInputs,'FaceColor','flat');
    pedalPlot.YDataSource = 'pedalInputs';
    
    pedalPlot(1).CData = [1 0 0; 0 1 0];
    ylim([0 100])
    ax = gca;
    disableDefaultInteractivity(ax)
    axis off

%     xtips2 = pedalPlot(1).XEndPoints;
%     ytips2 = [10 10];
%     labels2 = string(pedalPlot(1).YData);
%     delete(app.pedalLabel);
%     app.pedalLabel = text(app.UIAxes5, xtips2,ytips2,labels2,'HorizontalAlignment','center','FontSize',10);

    frame = getframe;
end

