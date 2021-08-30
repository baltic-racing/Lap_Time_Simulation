%% writeToLogfile.m
% Writes error and status messages to log file and to Simulation Status
% Textbox, if Debug is true and textAreaHandle is used as parameter.
%
% use writeToLogfile(text) or writeToLogfile(text, Debug) to only write to
% log file.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function writeToLogfile(text, Debug, textAreaHandle)
    %UNTITLED Summary of this function goes here
    %   Detailed explanation goes here
    fid = fopen('ErrorLog.txt', 'a');
    if fid == -1
      error('Cannot open log file.');
    end
    fprintf(fid, '%s: %s\n', datestr(now, 0), text);
    fclose(fid);
    
    if nargin == 3 && Debug
        textAreaHandle.Value{end+1} = 'loaded Track!';
    end
end

