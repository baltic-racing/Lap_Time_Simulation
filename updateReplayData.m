function updateReplayData(app, runNumber, result)
    
    val = app.Slider.Value;                                                 % Value of slider (current time)
    
    t = result(runNumber).t;                                                % Time variable of result

    ind = interp1(t,1:length(t),val,'nearest');                             % Find index of the time nearest to the time of the slider
    
    if ind == length(t)
        ind = ind-1;
    end

    app.SpeedLabel.Text = "Speed: " + num2str(result(runNumber).vV(ind));   % Update Current Speed Label
    
    app.RPMGauge.Value = result(runNumber).ni(ind);
    
    delete(app.xline);
    app.xline = xline(app.UIAxes4,ind,'r');
end