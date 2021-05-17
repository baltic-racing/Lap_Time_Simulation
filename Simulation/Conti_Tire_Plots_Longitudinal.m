%% Conti-Reifen Plots - Longitudinal
clear all
close all
clc

i = 0;

while i <= 12
    
    ALPHA = i;        % [°] Schräglaufwinkel
    GAMMA = 0;          % [°] Sturz
    KAPPA = -0.2:0.01:0.2;  % [-] Reifenschlupfmatrix, von 0% Schlupf bis 20% Schlupf
    
    % Variable für Pfad des Tir-Files
    FileNameLocation = ('C:\Users\Sven Weishaupt\Dropbox\MATLAB\_MATLAB_\C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

    % Tir-File in Struktur laden
    TIRparam = loadTIR(FileNameLocation);
    
    TIRparam.LMUX = 1;              % [-] Skalierungsfaktor der Reifen (0.75 für optimale Reifentemperatur, 0.6 für niedrige Reifentemperatur)
    TIRparam.RHX1 = 0;      % Einfluss ALPHA auf KAPPA
    TIRparam.LXAL = 1.2;
    TIRparam.RBX3 = 0;
    
    Fz_r = 1400;   % [N] Radlast hinten im Stand(1) -> statisch
    Fx = MF52_Fx_cs(ALPHA,Fz_r,GAMMA,KAPPA,TIRparam); % [N] Übertragbare Radkraft
    
    figure(1)
    hold on
    plot(KAPPA,Fx,'LineWidth',1)
    title('Übertragbare Längskraft in Abhängigkeit des Längsschlupfs und Schräglaufwinkels (FZ = 1400N)','FontSize',12)
    xlabel('Reifenschlupf [-]','FontSize',10)
    ylabel('Übertragbare Längskraft [N]','FontSize',10)
    xlim([-0.2 0.2])
    grid on
    box on
    
    i = i + 2;

end

legend('0° Schräglaufwinkel','2° Schräglaufwinkel','4° Schräglaufwinkel',...
    '6° Schräglaufwinkel','8° Schräglaufwinkel','10° Schräglaufwinkel',...
    '12° Schräglaufwinkel','FontSize',12,'Location','northwest')

        