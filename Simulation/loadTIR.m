function originalTyre = loadTIR(FileNameLocation)
%LOADTIR Load a user specified TIR file into a structure.
% 
% Syntax: [CurrentTyre] = loadTIR(FileNameLocation)
% where FileNameLocation is a string including full path and name

fileID = fopen(FileNameLocation,'r'); % Open file filename for reading ('r')

counter = 0;
LineNumber = 0;

%While not at end of file
while feof(fileID) == 0;
    
    currentLine = fgetl(fileID);
    LineNumber = LineNumber + 1;
    
    if ~isempty(currentLine)
        %If current line is not a title, comment or description
        if (currentLine(1) ~= '[' && currentLine(1) ~= '$' && currentLine(1) ~= '!')
            %Append parameter counter
            counter = counter + 1;
            
            %Find location of = and $ in each row
            index1 = strfind(currentLine, '=') + 1;
            index2 = strfind(currentLine, '$') - 1;
            
            %Parameter Name
            value = currentLine(1:index1-2);
            MFParamName{counter,1} = strtrim(value);
            %Parameter Value
            if (isempty(index2) == 1)
                value = currentLine(index1:end);
                MFParamValue{counter,1} = strtrim(value);
            else
                value = currentLine(index1:index2);
                MFParamValue{counter,1} = strtrim(value);
            end
            %Parameter Description
            value = currentLine(index2+1:end);
            MFParamDescription{counter,1} = strtrim(value);
            
            
        end
    end
    
end

fclose(fileID);

clearvars -except MFParamName MFParamValue MFParamDescription

% Mass is used as a variable name twice in the MF6.1 TIR
% format. Once in the Units section and once in the Inertia
% section. Rename the inertia section as MASS1. Undo this when
% writing the TIR file.
indSecondMass = find(strcmp(MFParamName,'MASS'),1,'last');
MFParamName(indSecondMass) = {'MASS1'};


%generate structure containing data of CurrentTyre
for ii = 1:length(MFParamName)
    if cellfun(@length,MFParamName(ii)) > 1
        if any(isstrprop(MFParamValue{ii},'digit')) > 0
            originalTyre.(MFParamName{ii}) = str2double(MFParamValue{ii});
        else
            originalTyre.(MFParamName{ii}) = MFParamValue{ii};
        end
    end
end
% inlcude the description information in a seperate field of
% the structure
for jj = 1:length(MFParamName)
    if cellfun(@length,MFParamName(jj)) > 1
        originalTyre.TIRDescription.(MFParamName{jj}) = MFParamDescription{jj};
    end
end

end
