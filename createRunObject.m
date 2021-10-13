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

    saveFileData = {'A_Akkuzelle' saveFile(runID).A_Akkuzelle(1:end); 'DRS_status' saveFile(runID).DRS_status(1:end-1);...
        'aVX' saveFile(runID).aVX(:)./9.81; 'BPPsignal' saveFile(runID).BPPsignal(:); 'E_Akku' saveFile(runID).E_Akku(1:end-1);...
        'E_Akku_Reku' saveFile(runID).E_Akku_Reku(1:end-1); 'E_Waerme' saveFile(runID).E_Waerme(1:end-1); 'E_ges' saveFile(runID).E_ges(1:end-1);...
        'vV' saveFile(runID).vV(1:end-1).*3.6; 'psi1' saveFile(runID).psi1(1:end-1); 'gear' saveFile(runID).gear(1:end-1);...
        'alpha_f' saveFile(runID).alpha_f(1:end-1); 'alpha_r' saveFile(runID).alpha_r(1:end-1); 'aRev' saveFile(runID).aRev(:);...
        'aVY' saveFile(runID).aVY(:)./9.81; 'beta' saveFile(runID).beta(1:end-1); 'FL' saveFile(runID).FL(:); 'FR' saveFile(runID).FR(:);...
        'vVYmax' saveFile(runID).vVYmax(1:end-1).*3.6; 'FVX' saveFile(runID).FVX(:); 'FVX_rl' saveFile(runID).FVX_hl(:); 'FVX_rr' saveFile(runID).FVX_hr(:);...
        'FVY' saveFile(runID).FVY(:); 'FWXmax_r' saveFile(runID).FWXmax_h(1:end-1);...
        'FWXmax_rl' saveFile(runID).FWXmax_hl(1:end-1); 'FWXmax_rr' saveFile(runID).FWXmax_hr(1:end-1); 'FWXmax_f' saveFile(runID).FWXmax_v(1:end-1);...
        'FWXmax_vl' saveFile(runID).FWXmax_vl(1:end-1); 'FWXmax_fr' saveFile(runID).FWXmax_vr(1:end-1); 'FWYf' saveFile(runID).FWYv(1:end-1);...
        'FWYr' saveFile(runID).FWYh(1:end-1); 'FWYmax_r' saveFile(runID).FWYmax_h(1:end-1);...
        'FWYmax_rl' saveFile(runID).FWYmax_hl(1:end-1); 'FWYmax_rr' saveFile(runID).FWYmax_hr(1:end-1); 'FWYmax_f' saveFile(runID).FWYmax_v(1:end-1);...
        'FWYmax_vl' saveFile(runID).FWYmax_vl(1:end-1); 'FWYmax_fr' saveFile(runID).FWYmax_vr(1:end-1); 'FWZr' saveFile(runID).FWZh(1:end-1);...
        'FWZf' saveFile(runID).FWZv(1:end-1); 'Faero' saveFile(runID).Faero(:); 'Fdrag' saveFile(runID).Fdr(:); 'Mi' saveFile(runID).Mi(:);...
        'P_M' saveFile(runID).P_M(:); 'P_el' saveFile(runID).P_el(:);...
        'Rdyn_rl' saveFile(runID).Rdyn_hl(1:end-1); 'Rdyn_rr' saveFile(runID).Rdyn_hr(1:end-1);...
        'Rdyn_fl' saveFile(runID).Rdyn_vl(1:end-1); 'Rdyn_fr' saveFile(runID).Rdyn_vr(1:end-1); 'RutschenY_r' saveFile(runID).RutschenY_h(1:end-1);...
        'RutschenY_f' saveFile(runID).RutschenY_v(1:end-1); 'TC' saveFile(runID).TC(1:end-1); 'Tirelimit' saveFile(runID).Tirelimit(1:end-1);...
        'cZ_rl' saveFile(runID).cZ_hl(1:end-1); 'cZ_rr' saveFile(runID).cZ_hr(1:end-1);...
        'cZ_fl' saveFile(runID).cZ_vl(1:end-1); 'cZ_fr' saveFile(runID).cZ_vr(1:end-1);...
        'dFWZrl_aero' saveFile(runID).dFWZhl_aero(:); 'dFWZrl_x' saveFile(runID).dFWZhl_x(:); 'dFWZrl_y' saveFile(runID).dFWZhl_y(:);...
        'dFWZrr_aero' saveFile(runID).dFWZhr_aero(:); 'dFWZrr_x' saveFile(runID).dFWZhr_x(:); 'dFWZrr_y' saveFile(runID).dFWZhr_y(:);...
        'dFWZfl_aero' saveFile(runID).dFWZvl_aero(:); 'dFWZfl_x' saveFile(runID).dFWZvl_x(:); 'dFWZfl_y' saveFile(runID).dFWZvl_y(:);...
        'dFWZfr_aero' saveFile(runID).dFWZvr_aero(:); 'dFWZfr_x' saveFile(runID).dFWZvr_x(:); 'dFWZfr_y' saveFile(runID).dFWZvr_y(:);...
        'kappa_fl' saveFile(runID).kappa_fl(1:end-1); 'kappa_fr' saveFile(runID).kappa_fr(1:end-1);...
        'kappa_rl' saveFile(runID).kappa_rl(1:end-1); 'kappa_rr' saveFile(runID).kappa_rr(1:end-1);...
        'l_contact_patch_fl' saveFile(runID).l_contact_patch_fl(1:end-1); 'l_contact_patch_fr' saveFile(runID).l_contact_patch_fr(1:end);...
        'l_contact_patch_rl' saveFile(runID).l_contact_patch_rl(1:end-1); 'l_contact_patch_rr' saveFile(runID).l_contact_patch_rr(1:end-1);...
        'motor_eff' saveFile(runID).motor_eff(1:end-1); 'ni' saveFile(runID).ni(:); 'vRev' saveFile(runID).vRev(:).*3.6;...
        'x-Coordinate' saveFile(runID).Track(1:end-1,1); 'y-Coordinate' saveFile(runID).Track(1:end-1,2); 'z-Coordinate' saveFile(runID).Track(1:end-1,3);...
        'Distance' saveFile(runID).Track(1:end-1,4);'Corner Radius' saveFile(runID).Track(1:end-1,5); 'vWoBrake' saveFile(runID).vWoBrake(1:end-1).*3.6};

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
        Data = saveFile(runID).Track(1:end-1,4)./1000;
    else               
        % add time variable
        Data = saveFile(runID).t(1:end-1);
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