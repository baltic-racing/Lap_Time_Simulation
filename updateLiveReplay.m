%% updateLiveReplayData.m
% Update the replay data of the Simulation Replay. 
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function updateLiveReplay(app, runNumber, result, ind)

    if isempty(result)
        return
    end
   
    delete(app.pointHandle2)                 % Deletes previous pointHandle
              
    x = result.Track(ind,1);
    y = result.Track(ind,2);
    app.pointHandle2 = plot(app.UIAxes3, x, y, 'Marker', 'o', 'LineWidth',5,'MarkerSize',11,'Color','r');   % Creates pointHandle (marks selected point from UITable in the UIAxes)

    if result(runNumber).ni(ind) >= result(runNumber).n_max
        app.RPMGauge.FontColor = 'r';
        app.RPMLabel_2.FontColor = 'r';
    else
        app.RPMGauge.FontColor = 'w';
        app.RPMLabel_2.FontColor = 'w';
    end

    app.SpeedLabel.Text = "Speed: " + num2str(round(result(runNumber).vV(ind).*3.6,2)) + " km/h";       % Update Current Speed Label
    app.GearLabel.Text = "Gear: " + num2str(result(runNumber).gearSelection(ind));                      % Update Current Gear Label
    app.RPMLabel.Text = "RPM: " + num2str(round(result(runNumber).ni(ind)));                            % Update Current RPM Label
    app.TimeLabel.Text = "Time: " + num2str(result(runNumber).t(ind)) + " s";                           % Update Current Time Label
    app.DistanceLabel.Text = "Distance: " + num2str(round(result(runNumber).Track(ind,4),2)) + " m";    % Update Current Distance Label
    app.RadiusLabel.Text = "Radius: " + num2str(round(result(runNumber).Track(ind,5),2)) + " m";        % Update Current Radius Label

    app.RPMGauge.Value = result(runNumber).ni(ind);

    app.var1Label.Text = num2str(app.DropDown.Value(ind));
    app.var2Label.Text = num2str(app.DropDown_2.Value(ind));
    app.var3Label.Text = num2str(app.DropDown_3.Value(ind));
    app.var4Label.Text = num2str(app.DropDown_4.Value(ind));
    app.var5Label.Text = num2str(app.DropDown_5.Value(ind));
    app.var6Label.Text = num2str(app.DropDown_6.Value(ind));
    app.var7Label.Text = num2str(app.DropDown_9.Value(ind));
    app.var8Label.Text = num2str(app.DropDown_10.Value(ind));
    app.var9Label.Text = num2str(app.DropDown_11.Value(ind));
    app.var10Label.Text = num2str(app.DropDown_12.Value(ind));
    app.var11Label.Text = num2str(app.DropDown_13.Value(ind));
    app.var12Label.Text = num2str(app.DropDown_14.Value(ind));

    gearing = 45;
    value = result(runNumber).beta(ind)*180/pi * gearing;
    angle = (abs(value)/360-fix(abs(value)/360)) * 360 * abs(value)/value;
    angle = round(angle/5)*5;
    if angle > 180 
        angle = angle - 180;
    elseif angle < -180
        angle = angle + 180;
    elseif isnan(angle)
        angle = 0;
    end

    app.Image4.ImageSource = "SteeringWheel" + angle + ".png";

    if result(runNumber).DRS_status(ind)
        app.DRSStatusLamp.Color = 'green';
    else
        app.DRSStatusLamp.Color = 'red';
    end

    if result(runNumber).BPPsignal(ind)
        app.BrakeLamp.Color = 'green';
    else
        app.BrakeLamp.Color = 'red';
    end
    
    delete(app.xline);
    app.xline = xline(app.UIAxes4,ind,'r');

    throttle = result(runNumber).P_M/max(result(runNumber).P_M);
    brake = result(runNumber).FB/max(result(runNumber).FB);

    if (ind > 1 && (round(brake(ind)*100,0) ~= round(brake(ind-1)*100,0) || round(throttle(ind)*100,0) ~= round(throttle(ind-1)*100,0)))
        pedalInputs = [round(brake(ind)*100,0) round(throttle(ind)*100,0)];
        delete(app.pedalLabel);
        app.pedalLabel = text(app.UIAxes5, [1 2],[10 10],strsplit(num2str(pedalInputs)),'HorizontalAlignment','center','FontSize',10);
        refreshdata(app.UIAxes5,'caller');      
    end

    drawnow limitrate;
end