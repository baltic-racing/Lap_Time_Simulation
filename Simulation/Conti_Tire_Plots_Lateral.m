%% Conti-Reifen Plots - Lateral
clear all
close all
clc

i = 0;

while i <= 0.2
    
    ALPHA = -12:0.1:12;   % [°] Schräglaufwinkel
    GAMMA = 0;            % [°] Sturz
    KAPPA = i;            % [-] Reifenschlupf
    
    % Variable für Pfad des Tir-Files
    FileNameLocation = ('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

    % Tir-File in Struktur laden
    TIRparam = loadTIR(FileNameLocation);
    
    TIRparam.LMUY = 1;              % [-] Skalierungsfaktor der Reifen (0.75 für optimale Reifentemperatur, 0.6 für niedrige Reifentemperatur)
    TIRparam.LYKA = 1.2;
    
    Fz_r = 1400;   % [N] Radlast hinten im Stand(1) -> statisch
    Fy = MF52_Fy_cs(ALPHA,Fz_r,GAMMA,KAPPA,TIRparam); % [N] Übertragbare Radkraft
    
    figure(1)
    hold on
    plot(ALPHA,Fy,'LineWidth',1)
    title('Übertragbare Querkraft in Abhängigkeit des Schräglaufwinkels und Längsschlupfs','FontSize',12)
    xlabel('Schräglaufwinkel [°]','FontSize',10)
    ylabel('Übertragbare Querkraft [N]','FontSize',10)
    xlim([-12 12])
    grid on
    box on
    
    i = i + 0.05;

end

legend('0% Schlupf','5% Schlupf','10% Schlupf','15% Schlupf','20% Schlupf','FontSize',12)
        