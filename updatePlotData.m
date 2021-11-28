function updatePlotData(app, runNumber, result)

    hold(app.UIAxes4,'on');

    if app.CheckBox.Value
        delete(app.var1plot);
        var1 = app.DropDown.Value;
        var1norm = var1/max(var1);
        app.var1plot = plot(app.UIAxes4,var1norm,'color',app.PlotColors(1));       % Plot Track in white.
    else
        hold(app.UIAxes4,'off');
        delete(app.var1plot);
        hold(app.UIAxes4,'on');
    end

    if app.CheckBox_2.Value
        delete(app.var2plot);
        var2 = app.DropDown_2.Value;
        var2norm = var2/max(var2);
        app.var2plot = plot(app.UIAxes4,var2norm,'color',app.PlotColors(2));       % Plot Track in white.
    else
        hold(app.UIAxes4,'off');
        delete(app.var2plot);
        hold(app.UIAxes4,'on');
    end

    if app.CheckBox_3.Value
        delete(app.var3plot);
        var3 = app.DropDown_3.Value;
        var3norm = var3/max(var3);
        app.var3plot = plot(app.UIAxes4,var3norm,'color',app.PlotColors(3));       % Plot Track in white.
    else
        hold(app.UIAxes4,'off');
        delete(app.var3plot);
        hold(app.UIAxes4,'on');
    end

    if app.CheckBox_4.Value
        delete(app.var4plot);
        var4 = app.DropDown_4.Value;
        var4norm = var4/max(var4);
        app.var4plot = plot(app.UIAxes4,var4norm,'color',app.PlotColors(4));       % Plot Track in white.
    else
        hold(app.UIAxes4,'off');
        delete(app.var4plot);
        hold(app.UIAxes4,'on');
    end

    if app.CheckBox_5.Value
        delete(app.var5plot);
        var5 = app.DropDown_5.Value;
        var5norm = var5/max(var5);
        app.var5plot = plot(app.UIAxes4,var5norm,'color',app.PlotColors(5));       % Plot Track in white.
    else
        hold(app.UIAxes4,'off');
        delete(app.var5plot);
        hold(app.UIAxes4,'on');
    end

    if app.CheckBox_6.Value
        delete(app.var6plot);
        var6 = app.DropDown_6.Value;
        var6norm = var6/max(var6);
        app.var6plot = plot(app.UIAxes4,var6norm,'color',app.PlotColors(6));       % Plot Track in white.
    else
        hold(app.UIAxes4,'off');
        delete(app.var6plot);
        hold(app.UIAxes4,'on');
    end
    
    updateReplayData(app, runNumber, result); 
    hold(app.UIAxes4,'off');
end