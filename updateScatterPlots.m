function updateScatterPlots(app, runNumber, result)
    

    scatterVar = app.DropDown_8.Value;

    scatter(app.UIAxes3,result(runNumber).Track(1:end-1,1),result(runNumber).Track(1:end-1,2),40,scatterVar,'filled');
    
    if ~app.AutoCheckBox.Value
        lowerLimit = app.LowerLimitEditField.Value;
        upperLimit = app.UpperLimitEditField.Value; 

        if lowerLimit >= upperLimit
            lowerLimit = 0;
            upperLimit = 1;
        end
    else
        lowerLimit = min(scatterVar);
        upperLimit = max(scatterVar); 

        if lowerLimit >= upperLimit
            lowerLimit = 0;
            upperLimit = 1;
        end
    end

    set(app.UIAxes3,'clim',[lowerLimit upperLimit]);
end