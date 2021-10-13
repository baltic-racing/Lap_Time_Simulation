%% loadTrack.m
% Loads the track given by the TrackFileName and returns all of its data.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [x_Track, y_Track, z_Track, s, R, Track, ApexIndexes, lapLength] = loadTrack(TrackFileName, disciplineID, numOfLaps)
%% Loads the track for the simulation    

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


        trackLength = length(Track);

        for i = 1:trackLength
            x(i*3-2) = x_Track(i);
            x(i*3-1) = x_Track(i);
            x(i*3) = x_Track(i);

            y(i*3-2) = y_Track(i);
            y(i*3-1) = y_Track(i);
            y(i*3) = y_Track(i);

            z(i*3-2) = z_Track(i);
            z(i*3-1) = z_Track(i);
            z(i*3) = z_Track(i);

            s(i*3-2) = s(i);
            s(i*3-1) = s(i);
            s(i*3) = s(i);

            R(i*3-2) = R(i);
            R(i*3-1) = R(i);
            R(i*3) = R(i);
        end

        Track = [x, y, z, s, R];  
    end   
end