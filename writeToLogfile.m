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

