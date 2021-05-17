clear all
close all

load('Aero Downforce Daten.mat');    % Aero Daten
% Variable für Pfad des Tir-Files
FileNameLocation = ('C:\Users\Sven Weishaupt\Dropbox\MATLAB\_MATLAB_\C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

% Tir-File in Struktur laden
TIRparam = loadTIR(FileNameLocation);
R = 10;
m = 280;
FZ_stat = 800;
k = 0;

for i = 1:8
    
    v1 = 5;
    v2 = 100;
    FZ_stat = 800 + (i-1)*200;

    while v1 < v2
        
        v1 = v1 + k;
        
        Faero = interp1(v,FA,v1*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
        
        FZ(i) = FZ_stat + Faero/2;
        
        v2 = sqrt(max(MF52_Fy_cs(0:-0.1:-12,FZ(i),0,0,TIRparam))*R/m);
        
        k = k+0.0001;
        
    end
    
    aQmax(i) = v2^2/R/9.81;

end

scatter(FZ,aQmax)
