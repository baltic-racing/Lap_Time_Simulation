classdef SimulationTests < matlab.unittest.TestCase

    methods(TestClassSetup)
        % Shared setup for the entire test class     
    end

    methods(TestMethodSetup)
        % Setup for each test
    end

    methods(Test)
        % Test methods
       
        function accelTest(testCase)
            startingParameters = load('Tests/Setups/startingParametersAccel.mat');
            startingParameters = startingParameters.startingParameters;

            Vehiclesim_Endurance_GUI_Version(startingParameters)
        end

        function sensitivityAccelTest(testCase)
            startingParameters = load('Tests/Setups/startingParametersAccelSensitivity.mat');
            startingParameters = startingParameters.startingParameters;

            result = Vehiclesim_Endurance_GUI_Version(startingParameters);
            % ToDo: assert that time decreases
            % assert(size(result, 2) == 3); % ToDo: implement correctly     
        end
        
        function sensitivityAccelTestTwoParams(testCase)
            startingParameters = load('Tests/Setups/startingParametersAccelSensitivityTwoParams.mat');
            startingParameters = startingParameters.startingParameters;

            Vehiclesim_Endurance_GUI_Version(startingParameters);
            % assert(size(result, 2) == 3); % ToDo: implement correctly     
        end
        
        function ottobianoEndurance(testCase)
            startingParameters = load('Tests/Setups/startingParametersOttobianoEndu.mat');
            startingParameters = startingParameters.startingParameters;

            Vehiclesim_Endurance_GUI_Version(startingParameters);
            % assert(size(result, 2) == 3); % ToDo: implement correctly     
        end

        % function testCustomPlots(testCase)
        %     startingParameters = load('Tests/Setups/startingParametersAccelSensitivity.mat');
        %     startingParameters = startingParameters.startingParameters;
        % 
        %     result = Vehiclesim_Endurance_GUI_Version(startingParameters);
        % 
        %     for plotID = 1:43
        %         PlotResults(result, plotID, 1, 1)
        %     end
        % end

        % ADD GUI test cases and t
    end

end