%% loadTrack.m
% Loads the track given by the TrackFileName and returns all of its data.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [x_Track, y_Track, z_Track, s, R, Track, ApexIndexes, lapLength] = loadTrack(TrackFileName, disciplineID, numOfLaps, expandData)
%% Loads the track for the simulation    

    % DEBUG!
    expandData = 1;
        
    % Loads all the track variables from the Trackfile
    load(TrackFileName,'Track')
    
    % Length of the track before adding laps for endurance
    lapLength = length(Track);

    % Checks if AutoX or Endurance is selected (1==AutoX,2==Endurance)
    if disciplineID == 1    % AutoX
        
        x_Track = Track(:,1);           % [m] X-Coordinate of the Track
        y_Track = Track(:,2);           % [m] Y-Coordinate of the Track
        z_Track = Track(:,3);           % [m] Z-Coordinate of the Track
        s = Track(:,4);                 % [m] Track Pathway (Verlauf der Streckenlänge)
        R = Track(:,5);                 % [m] Radius of curves (Kurvenradien)

        % Calls the .m file 'Apexes' to calculate the Apexes of the given track
        ApexIndexes = Apexes(abs(R));  
        
    elseif disciplineID == 2    % Endurance
        x_Track = [];
        y_Track = [];
        z_Track = [];
        s = [];
        R = [];
        
        s_max = max(Track(:,4));
        
        % Expand Track to Endurance Distance
        for i = 1:numOfLaps
            x_Track = [x_Track; Track(:,1)];     % [m] X-Coordinate of the Track
            y_Track = [y_Track; Track(:,2)];     % [m] Y-Coordinate of the Track
            z_Track = [z_Track; Track(:,3)];     % [m] Z-Coordinate of the Track
            s = [s; Track(:,4)+s_max*(i-1)];                 % [m] Track Pathway (Verlauf der Streckenlänge)
            R = [R; Track(:,5)];                 % [m] Radius of curves (Kurvenradien)

            % Calls the .m file 'Apexes' to calculate the Apexes of the given track
            ApexIndexes = Apexes(abs(R));  
        end
        
        % Complete Endurance Track (All Laps)
        Track = [x_Track, y_Track, z_Track, s, R];  

    end   
    
    if expandData
        trackLength = length(Track);

        s(1) = 0;

        for i = 1:trackLength-1
            x(i*3-2) = x_Track(i);
            x(i*3-1) = x_Track(i);
            x(i*3) = x_Track(i);

            y(i*3-2) = y_Track(i);
            y(i*3-1) = y_Track(i);
            y(i*3) = y_Track(i);

            z(i*3-2) = z_Track(i);
            z(i*3-1) = z_Track(i);
            z(i*3) = z_Track(i);

            if i ~= 1
                s_step = (s(i+1)-s(i))/3;
                s_new(i*3-2) = s(i);
                s_new(i*3-1) = s(i)+s_step;
                s_new(i*3) = s(i)+s_step*2;
            else
                s_step = (s(2)-s(1))/3;
                s_new(i*3-2) = 0;
                s_new(i*3-1) = 0+s_step;
                s_new(i*3) = 0+s_step*2;
            end

            R_new(i*3-2) = R(i);
            R_new(i*3-1) = R(i);
            R_new(i*3) = R(i);
        end

        Track = [x', y', z', s_new', R_new']; 

        x_Track = Track(:,1);           % [m] X-Coordinate of the Track
        y_Track = Track(:,2);           % [m] Y-Coordinate of the Track
        z_Track = Track(:,3);           % [m] Z-Coordinate of the Track
        s = Track(:,4);                 % [m] Track Pathway (Verlauf der Streckenlänge)
        R = Track(:,5);                 % [m] Radius of curves (Kurvenradien)
        ApexIndexes = ApexIndexes.*3;
    end
end