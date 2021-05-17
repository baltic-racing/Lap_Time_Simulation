%% Excel Einlesen
% Setup the Import Options
opts = spreadsheetImportOptions("NumVariables", 2);

% Specify sheet and range
opts.Sheet = "Tabelle1";
opts.DataRange =  "A2:B8" ;

% Specify column names and types
opts.VariableTypes = ["double", "double"];

% Import the data
Tabelle = readtable("D:\Users\Sven Weishaupt\Desktop\E-Auto Getriebe\Matlab Motor Auslegung\EMRAX Motor\Drehmoment_Drehzahl_Tabelle.xlsx", opts, "UseExcel", false);

% Clear temporary variables
clear opts

% Tabellen Excelwert als array 
Tabelle = table2array(Tabelle);
% Erste Variable
n = Tabelle(:,1);
% Zweite Variable 
M = Tabelle(:,2);

%% Mat Datei aus Daten erstellen
save('EMRAX_208_TorqueData_Peak.mat','n','M');