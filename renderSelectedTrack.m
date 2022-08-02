%% renderSelectedTrack.m
% Plots the selected racetrack from the scenario panel track dropdown
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function renderSelectedTrack(app)
            
    app.TrackFileName = app.Tracks.Data(app.TrackDropDown.Value,7);             % Get Track filename with the track name from the track dropdown.

    load(app.TrackFileName,'Track');                                            % Load track with given TrackFileName
    x_Track_preview = Track(:,1);                                               % [m] X-Koordinate der Strecke
    y_Track_preview = Track(:,2);                                               % [m] Y-Koordinate der Strecke
    s_preview = Track(:,4);                                                     % [m] Verlauf der Streckenlänge
    TrackName = app.Tracks.Data(app.TrackDropDown.Value,1);                     % Name of the race track
    app.StartingSpeed = str2double(app.Tracks.Data(app.TrackDropDown.Value,6));
    app.numOfLaps = str2double(app.Tracks.Data(app.TrackDropDown.Value,5));

    % Enables the selection dropdown for the discipline
    if str2double(app.Tracks.Data(app.TrackDropDown.Value,4))
        app.DisciplineDropDown.Enable = true;
    else
        app.DisciplineDropDown.Enable = false;
    end

    app.UIAxes2.reset; 

    if app.DarkModeCheckBox.Value
        plot(app.UIAxes2,x_Track_preview,y_Track_preview,'w')
    else
        plot(app.UIAxes2,x_Track_preview,y_Track_preview,'k')   
    end
    
    if size(Track,2) > 5
        x_inner = Track(:,6);
        y_inner = Track(:,7);
        
        x_outer = Track(:,8);
        y_outer = Track(:,9);
        
        hold(app.UIAxes2,'on')
        plot(app.UIAxes2,x_inner,y_inner,'b');
        plot(app.UIAxes2,x_outer,y_outer,'b');
    end    
    
    title(app.UIAxes2,TrackName,'FontSize',12);
    app.UIAxes2.LineWidth = 1.5;

    % Sets Color to Dark or white Mode
    app.UIAxes2.Title.Color = app.OnPrimaryColor;
    app.UIAxes2.XColor = app.OnPrimaryColor;
    app.UIAxes2.YColor = app.OnPrimaryColor;
    app.UIAxes2.Color = app.SurfaceColor;       

    app.LengthLabel.Text = "Track Length: " + max(s_preview) + " m";

    drawApexes = 1; % ToDo: add as option in GUI
    if drawApexes
        hold(app.UIAxes2,"on");
        ApexIndizes = Apexes(Track(:,5));
        scatter(app.UIAxes2,x_Track_preview(ApexIndizes),y_Track_preview(ApexIndizes),20,'r','filled')
    end
end            