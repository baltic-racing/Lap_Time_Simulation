function SimulationTest()
    %% Test if setup is valid
    setup1 = 'TY19_LowDownforce.mat';
    path1 = 'E:\Laptime Simulation\Lap_Time_Simulation\Presets';
    result1 = Vehiclesim_Endurance_GUI_Version(setup1, path1, 'TrackAccel.mat', 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1);
    
    setup2 = 'TY19_HighDownforce.mat';
    path2 = 'E:\Laptime Simulation\Lap_Time_Simulation\Presets';
    result2 = Vehiclesim_Endurance_GUI_Version(setup2, path2, 'TrackAccel.mat', 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1);
    
    numSteps = 10;
    
    parfor(i = 1:numSteps)
        result(:,i) = Vehiclesim_Endurance_GUI_Version(setup2, path2, 'TrackAccel.mat', 1, 6, i/2, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 1, 1);
    end
    
    plot(result(6).ni)
    
    %% Usage of results
    % x=result(4).ni -> result(runNumber).ni
    
%     % combine results
%     for i = 1:numSteps
%         %result2(:,i) = result(i); 
%         result2(i,:) = result(i); 
%     end

    %% Test 1
    assert(result1.t(end)>result2.t(end),'First run slower!')

    %% Test 2
    
end

