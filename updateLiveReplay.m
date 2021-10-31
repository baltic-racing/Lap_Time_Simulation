function updateLiveReplay(app, runNumber, result, ind)

    delete(app.pointHandle2)                 % Deletes previous pointHandle
            
    hold(app.UIAxes3,'on')
    x = result.Track(ind,1);
    y = result.Track(ind,2);
    app.pointHandle2 = plot(app.UIAxes3, x, y, 'Marker', 'o', 'LineWidth',5,'MarkerSize',11,'Color','r');   % Creates pointHandle (marks selected point from UITable in the UIAxes)

    %comet(app.UIAxes3, x, y, 0);

    app.SpeedLabel.Text = "Speed: " + num2str(round(result(runNumber).vV(ind).*3.6,2)) + " km/h";   % Update Current Speed Label
    app.GearLabel.Text = "Gear: " + num2str(result(runNumber).gearSelection(ind));   % Update Current Speed Label
    app.RPMLabel.Text = "RPM: " + num2str(round(result(runNumber).ni(ind)));   % Update Current Speed Label
    app.TimeLabel.Text = "Time: " + num2str(result(runNumber).t(ind)) + " s";   % Update Current Speed Label
    app.DistanceLabel.Text = "Distance: " + num2str(result(runNumber).Track(ind,4)) + " m";   % Update Current Speed Label
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
    
    drawnow limitrate;
end