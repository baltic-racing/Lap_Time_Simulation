function plotDataAnalysisPlot(app, event)
    displayRow = event.InteractionInformation.DisplayRow;
    displayColumn = event.InteractionInformation.DisplayColumn;
    
    % Return if column headline have been clicked
    if isempty(displayRow)
        return
    end

    varName = app.UITable.Data(displayRow,1);
    result1 = app.results.(['result' num2str(app.activeResult(1))]);

    if event.EventName == "Clicked"
        

        if app.HoldonCheckBox.Value
            hold(app.UIAxes_3, "on")

            % Save old label text to append new text
            oldText = app.plotNameLabel.Value;
            app.plotNameLabel.Value = string([oldText; varName]);
        else
            hold(app.UIAxes_3, "off")
            app.plotNameLabel.Value = string(varName);
        end

        plot(app.UIAxes_3, result1.(string(varName)))
        hold(app.UIAxes_3, "on")

        if length(app.activeResult) > 1
            result2 = app.results.(['result' num2str(app.activeResult(2))]);
            plot(app.UIAxes_3, result2.(string(varName)))
        end

        hold(app.UIAxes_3, "off")
    else
        plot(app.results.result1.t, result1.(string(varName)))
        hold on

        if length(app.activeResult) > 1
            result2 = app.results.(['result' num2str(app.activeResult(2))]);
            plot(app.results.result1.t, result2.(string(varName)))
        end
        title(varName)

        hold off
    end
end