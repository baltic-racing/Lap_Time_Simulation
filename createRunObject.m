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

    saveFileData = {'A_Accu_cell' saveFile(runID).A_accu_cell(1:end-1); 'DRS_status' saveFile(runID).DRS_status(1:end-1);...
        'aVX' saveFile(runID).aVX(:)./9.81; 'BPPsignal' saveFile(runID).BPPsignal(:); 'E_Accu' saveFile(runID).E_Accu(1:end-1);...
        'E_Akku_Reku' saveFile(runID).E_Accu_Recu(1:end-1); 'E_heat' saveFile(runID).E_heat(1:end-1); 'E_res' saveFile(runID).E_res(1:end-1);...
        'vV' saveFile(runID).vV(1:end-1).*3.6; 'psi1' saveFile(runID).psi1(1:end-1); 'gear' saveFile(runID).gearSelection(1:end);...
        'alpha_f' saveFile(runID).alpha_f(1:end-1); 'alpha_r' saveFile(runID).alpha_r(1:end-1); 'aRev' saveFile(runID).aRev(:);...
        'aVY' saveFile(runID).aVY(:)./9.81; 'beta' saveFile(runID).beta(1:end-1); 'FL' saveFile(runID).FL(:); 'FR' saveFile(runID).FR(:);...
        'vVYmax' saveFile(runID).vVYmax(1:end-1).*3.6; 'FVX' saveFile(runID).FVX(:); 'FVX_rl' saveFile(runID).FVX_rl(:); 'FVX_rr' saveFile(runID).FVX_rr(:);...
        'FVY' saveFile(runID).FVY(:); 'FWXmax_r' saveFile(runID).FWXmax_r(1:end-1);...
        'FWXmax_rl' saveFile(runID).FWXmax_fl(1:end-1); 'FWXmax_rr' saveFile(runID).FWXmax_rr(1:end-1); 'FWXmax_f' saveFile(runID).FWXmax_f(1:end-1);...
        'FWXmax_fl' saveFile(runID).FWXmax_fl(1:end-1); 'FWXmax_fr' saveFile(runID).FWXmax_fr(1:end-1); 'FWYf' saveFile(runID).FWYf(1:end);...
        'FWYr' saveFile(runID).FWYr(1:end); 'FWYmax_r' saveFile(runID).FWYmax_r(1:end-1);...
        'FWYmax_rl' saveFile(runID).FWYmax_rl(1:end-1); 'FWYmax_rr' saveFile(runID).FWYmax_rr(1:end-1); 'FWYmax_f' saveFile(runID).FWYmax_f(1:end-1);...
        'FWYmax_fl' saveFile(runID).FWYmax_fl(1:end-1); 'FWYmax_fr' saveFile(runID).FWYmax_fr(1:end-1); 'FWZr' saveFile(runID).FWZr(1:end-1);...
        'FWZf' saveFile(runID).FWZf(1:end-1); 'Faero' saveFile(runID).Faero(:); 'Fdrag' saveFile(runID).Fdr(:); 'Mi' saveFile(runID).Mi(:);...
        'P_M' saveFile(runID).P_M(:); 'P_el' saveFile(runID).P_el(:);...
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
        'alpha_rl' saveFile(runID).alpha_rl(1:end);'alpha_rr' saveFile(runID).alpha_rr(1:end); 'alpha_fr' saveFile(runID).alpha_fr(1:end); 'alpha_fl' saveFile(runID).alpha_fl(1:end)};
    
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