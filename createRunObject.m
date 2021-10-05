%% createRunObject.m
% Creates an object for the Simulink Data Inspector from an result .mat
% file.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function createRunObject(app, saveFile, saveID)
            
    filename = app.ResultfilenameEditField.Value;

    runID = app.SelectLapRunnumberEditField.Value;          

    % checks if the selected runID is a valid run
    if checkRunID(app, runID, saveFile)
        return;
    end

    % Creates Name for the Run if no name was defined
    if app.ResultfilenameEditField.Value == ""
        runName = ['Result ' num2str(saveID)];
    else
        runName = filename;
    end

    saveFileData = {'A_Akkuzelle' saveFile.A_Akkuzelle(1:end-1,runID); 'DRS_status' saveFile.DRS_status(1:end-1,runID);...
        'aVX' saveFile.aVX(:,runID)./9.81; 'BPPsignal' saveFile.BPPsignal(:,runID); 'E_Akku' saveFile.E_Akku(1:end-1,runID);...
        'E_Akku_Reku' saveFile.E_Akku_Reku(1:end-1,runID); 'E_Waerme' saveFile.E_Waerme(1:end-1,runID); 'E_ges' saveFile.E_ges(1:end-1,runID);...
        'FB' saveFile.FB(1:end-1,runID); 'vV' saveFile.vV(1:end-1,runID).*3.6; 'psi1' saveFile.psi1(1:end-1,runID); 'gear' saveFile.gear(:,runID);...
        'alpha_f' saveFile.alpha_f(1:end-1,runID); 'alpha_r' saveFile.alpha_r(1:end-1,runID); 'aRev' saveFile.aRev(:,runID);...
        'aVY' saveFile.aVY(:,runID)./9.81; 'beta' saveFile.beta(1:end-1,runID); 'FL' saveFile.FL(:,runID); 'FR' saveFile.FR(:,runID);...
        'vVYmax' saveFile.vVYmax(1:end-1,runID).*3.6; 'FVX' saveFile.FVX(:,runID); 'FVX_rl' saveFile.FVX_hl(:,runID); 'FVX_rr' saveFile.FVX_hr(:,runID);...
        'FVXid' saveFile.FVXid(:,runID); 'FVXre' saveFile.FVXre(:,runID); 'FVY' saveFile.FVY(:,runID); 'FWXmax_r' saveFile.FWXmax_h(1:end-1,runID);...
        'FWXmax_rl' saveFile.FWXmax_hl(1:end-1,runID); 'FWXmax_rr' saveFile.FWXmax_hr(1:end-1,runID); 'FWXmax_f' saveFile.FWXmax_v(1:end-1,runID);...
        'FWXmax_vl' saveFile.FWXmax_vl(1:end-1,runID); 'FWXmax_fr' saveFile.FWXmax_vr(1:end-1,runID); 'FWYf' saveFile.FWYv(:,runID);...
        'FWYr' saveFile.FWYh(:,runID); 'FWYmax_r' saveFile.FWYmax_h(1:end-1,runID);...
        'FWYmax_rl' saveFile.FWYmax_hl(1:end-1,runID); 'FWYmax_rr' saveFile.FWYmax_hr(1:end-1,runID); 'FWYmax_f' saveFile.FWYmax_v(1:end-1,runID);...
        'FWYmax_vl' saveFile.FWYmax_vl(1:end-1,runID); 'FWYmax_fr' saveFile.FWYmax_vr(1:end-1,runID); 'FWZr' saveFile.FWZh(1:end-1,runID);...
        'FWZf' saveFile.FWZv(1:end-1,runID); 'Faero' saveFile.Faero(:,runID); 'Fdrag' saveFile.Fdr(:,runID); 'Mi' saveFile.Mi(:,runID);...
        'P_M' saveFile.P_M(:,runID); 'P_el' saveFile.P_el(:,runID);...
        'Rdyn_rl' saveFile.Rdyn_hl(1:end-1,runID); 'Rdyn_rr' saveFile.Rdyn_hr(1:end-1,runID);...
        'Rdyn_fl' saveFile.Rdyn_vl(1:end-1,runID); 'Rdyn_fr' saveFile.Rdyn_vr(1:end-1,runID); 'RutschenY_r' saveFile.RutschenY_h(1:end-1,runID);...
        'RutschenY_f' saveFile.RutschenY_v(1:end-1,runID); 'TC' saveFile.TC(1:end-1,runID); 'Tirelimit' saveFile.Tirelimit(1:end-1,runID);...
        'cZ_rl' saveFile.cZ_hl(1:end-1,runID); 'cZ_rr' saveFile.cZ_hr(1:end-1,runID);...
        'cZ_fl' saveFile.cZ_vl(1:end-1,runID); 'cZ_fr' saveFile.cZ_vr(1:end-1,runID);...
        'dFWZrl_aero' saveFile.dFWZhl_aero(:,runID); 'dFWZrl_x' saveFile.dFWZhl_x(:,runID); 'dFWZrl_y' saveFile.dFWZhl_y(:,runID);...
        'dFWZrr_aero' saveFile.dFWZhr_aero(:,runID); 'dFWZrr_x' saveFile.dFWZhr_x(:,runID); 'dFWZrr_y' saveFile.dFWZhr_y(:,runID);...
        'dFWZfl_aero' saveFile.dFWZvl_aero(:,runID); 'dFWZfl_x' saveFile.dFWZvl_x(:,runID); 'dFWZfl_y' saveFile.dFWZvl_y(:,runID);...
        'dFWZfr_aero' saveFile.dFWZvr_aero(:,runID); 'dFWZfr_x' saveFile.dFWZvr_x(:,runID); 'dFWZfr_y' saveFile.dFWZvr_y(:,runID);...
        'kappa_fl' saveFile.kappa_fl(1:end-1,runID); 'kappa_fr' saveFile.kappa_fr(1:end-1,runID);...
        'kappa_rl' saveFile.kappa_rl(1:end-1,runID); 'kappa_rr' saveFile.kappa_rr(1:end-1,runID);...
        'l_contact_patch_fl' saveFile.l_contact_patch_fl(1:end-1,runID); 'l_contact_patch_fr' saveFile.l_contact_patch_fr(1:end-1,runID);...
        'l_contact_patch_rl' saveFile.l_contact_patch_rl(1:end-1,runID); 'l_contact_patch_rr' saveFile.l_contact_patch_rr(1:end-1,runID);...
        'motor_eff' saveFile.motor_eff(1:end-1,runID); 'ni' saveFile.ni(:,runID); 'vRev' saveFile.vRev(:,runID).*3.6;...
        'x-Coordinate' saveFile.Track(1:end-1,1); 'y-Coordinate' saveFile.Track(1:end-1,2); 'z-Coordinate' saveFile.Track(1:end-1,3);...
        'Distance' saveFile.Track(1:end-1,4);'Corner Radius' saveFile.Track(1:end-1,5); 'vWoBrake' saveFile.vWoBrake(1:end-1,runID).*3.6};

    names = saveFileData;
    names(:,1) = lower(names(:,1));

    [~, idx] = sortrows(names, 1); 

    % Sorts the table alphabeticaly case insensitiv
    saveFileSorted = saveFileData(idx, :);

    % Create textheader by inserting commas and adding time (also
    % needed for distance because the Data Inspector uses only time
    % for x-axis
    textHeader = strcat('time,',strjoin(saveFileSorted(:,1), ','));

    % Checks if the data should be exported as time or distance
    % based
    if app.ExportasdistanceCheckBox.Value
        % add distance variable
        Data = saveFile.Track(1:end-1,4)./1000;
    else               
        % add time variable
        Data = saveFile.t(1:end-1,runID);
    end

    % Write empty csv file
    writematrix('', 'import.csv');


    % Write data to an array
    for i = 1:length(saveFileData)
        Data(:,i+1) = cell2mat(saveFileSorted(i,2));
    end

    %write header to file
    fid = fopen('import.csv','w'); 
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    %write data to end of file
    dlmwrite('import.csv',Data,'-append');
    %writetable(Data,'import.csv','WriteMode','append','Delimiter','comma')

    %% Creating Data for Simulink Data Inspector
    csvRunID = Simulink.sdi.createRun(runName,'file','import.csv');

    if app.DeletecsvfileCheckBox.Value
        delete('import.csv');
    end    
end