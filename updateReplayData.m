function updateReplayData(app, runNumber, result)
    
    val = app.Slider.Value;                                                 % Value of slider (current time)
    
    t = result(runNumber).t;                                                % Time variable of result

    ind = interp1(t,1:length(t),val,'nearest');                             % Find index of the time nearest to the time of the slider
    
    if ind == length(t)
        ind = ind-1;
    end

    delete(app.pointHandle2)                 % Deletes previous pointHandle

    if result(runNumber).ni(ind) >= result(runNumber).n_max
        app.RPMGauge.FontColor = 'r';
        app.RPMLabel_2.FontColor = 'r';
    else
        app.RPMGauge.FontColor = 'w';
        app.RPMLabel_2.FontColor = 'w';
    end
            
    hold(app.UIAxes3,'on')
    x = result.Track(ind,1);
    y = result.Track(ind,2);
    app.pointHandle2 = plot(app.UIAxes3, x, y, 'Marker', 'o', 'LineWidth',5,'MarkerSize',11,'Color','r');   % Creates pointHandle (marks selected point from UITable in the UIAxes)         

    app.SpeedLabel.Text = "Speed: " + num2str(round(result(runNumber).vV(ind).*3.6,2)) + " km/h";   % Update Current Speed Label
    app.GearLabel.Text = "Gear: " + num2str(result(runNumber).gearSelection(ind));   % Update Current Speed Label
    app.RPMLabel.Text = "RPM: " + num2str(round(result(runNumber).ni(ind)));   % Update Current Speed Label
    app.TimeLabel.Text = "Time: " + num2str(result(runNumber).t(ind)) + " s";   % Update Current Speed Label
    app.DistanceLabel.Text = "Distance: " + num2str(round(result(runNumber).Track(ind,4),2)) + " m";   % Update Current Speed Label
    app.RadiusLabel.Text = "Radius: " + num2str(round(result(runNumber).Track(ind,5),2)) + " m";   % Update Current Speed Label

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
   
    %% Display wheel load in percent
    FWZ_fl = result(runNumber).FWZ_fl(ind);
    FWZ_fr = result(runNumber).FWZ_fr(ind);
    FWZ_rl = result(runNumber).FWZ_rl(ind);
    FWZ_rr = result(runNumber).FWZ_rr(ind);

    wheelload = FWZ_fl + FWZ_fr + FWZ_rl + FWZ_rr;

    percentfl = num2str(round(FWZ_fl/wheelload*100));
    percentfr = num2str(round(FWZ_fr/wheelload*100));
    percentrl = num2str(round(FWZ_rl/wheelload*100));
    percentrr = num2str(round(FWZ_rr/wheelload*100));

    StringW = [percentfl '      ' percentfr newline percentrl '      ' percentrr];

    app.wheelloadLabel.Text = StringW;


    %% Display DRS and Brake Status
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

    pedalInputs = [round(brake(ind)*100,0) round(throttle(ind)*100,0)];
    delete(app.pedalLabel);
    app.pedalLabel = text(app.UIAxes5, [1 2],[10 10],strsplit(num2str(pedalInputs)),'HorizontalAlignment','center','FontSize',10);
    
    refreshdata(app.UIAxes5,'caller');
    drawnow;

end