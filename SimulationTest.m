function [assertResult1] = SimulationTest()
    %% Test if setup is valide
    setup1 = 'TY19_LowDownforce.mat';
    path1 = 'C:\Users\Eric\Desktop\Simulation_Git\Lap_Time_Simulation\Presets';
    result1 = Vehiclesim_Endurance_GUI_Version(setup1, path1, 'TrackAccel.mat', 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1);
    
    setup2 = 'TY19_HighDownforce.mat';
    path2 = 'C:\Users\Eric\Desktop\Simulation_Git\Lap_Time_Simulation\Presets';
    result2 = Vehiclesim_Endurance_GUI_Version(setup2, path2, 'TrackAccel.mat', 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 1, 1);
    
    assert(result1.t(end)>result2.t(end),'First run slower!')
    
    x1 = result1.t(end)
    x2 = result2.t(end)

    %% Test single class
    exp = 'single';
    act = ones('single');
    assert(isa(act,exp))
end

