%% Acceleration Simulation EMRAX Motor
clear all
close all
clc

%% Matfiles einlesen
load('Aero Downforce Daten.mat');    % Aero Daten
load('EMRAX_208_TorqueData_Peak.mat');    % E-Motor Daten

%% Fahrzeugdaten   
m_ges = 265;    % [kg] Fahrzeuggesamtmasse inkl. Fahrer
g = 9.81;       % [m/s²] Erdbeschleunigung
t_L = 25;       % [°C] Umgebungslufttemperatur
p_L = 1.013;    % [bar] Umgebungsluftdruck
R_L = 287;      % [J/(kg*K)] Gaskonstante Luft
l = 1525;       % [mm] Radstand
h_COG = 280;    % [mm] Höhe Fahrzeugschwerpunkt
m_ph = 55;      % [%] Prozentualer Radlastanteil der Hinterachse
k_R = 0.01;     % [-] Rollwiderstandsbeiwert
c_w = 1.23;    % [-] cw-Wert des Fahrzeugs
A_S = 1.086;      % [m²] Stirnfläche des Wagens
x_vA = 470;     % [mm] x-Koordinate der Vorderachse in CAD
x_COP = 1600;   % [mm] x-Koordinate Center of Pressure in CAD

%% Motor- und Getriebedaten
i_ges = 3.88;          % [-] Gesamtübersetung Motor zu Rad
P_Mmax = 68000*2;     % [W] Maximalleistung des Motors

%% Reifenmodell
% Variable für Pfad des Tir-Files
FileNameLocation = ('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

% Tir-File in Struktur laden
TIRparam = loadTIR(FileNameLocation);

ALPHA = 0;          % [°] Schräglaufwinkel
GAMMA = 0;          % [°] Sturz 
KAPPA = 0:0.01:0.2; % [-] Längsschlupfmatrix für Reifenmodell

% Vertikal
FZ0   = TIRparam.FNOMIN;          % [N] Nominelle Radlast
R0    = TIRparam.UNLOADED_RADIUS; % [m] Fertigungsradius des Reifens
TIRparam.LMUX = 0.6;              % [-] Skalierungsfaktor der Reifen (0.75 für optimale Reifentemperatur, 0.6 für niedrige Reifentemperatur)

% Weitere Reifenparameter
c_tire = TIRparam.VERTICAL_STIFFNESS;          % [N/m] Vertikale Steifigkeit des Reifens
Fz_r(1) = m_ph/100*m_ges * g/2;                % [N] Radlast hinten im Stand(1) -> statisch
r_dyn(1) = R0 - Fz_r(1) /c_tire;               % [m] Dynamischer Reifenradius im Stand(1) -> statisch
Fx_Achse(1) = 2*max(MF52_Fx_ps(ALPHA,Fz_r,GAMMA,KAPPA,TIRparam)); % [N] Übertragbare Achskraft, insbesondere abhängig vom Schlupf

%% Vorbereitungen zur Rechnung
delta_t = 0.001;                        % [s] Zeitschritt für Simulation
rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m³] Luftdichte

% Startwerte
v_V(1) = 0;             % [m/s] Anfangsgeschwindigkeit
s(1) = 0;               % [m] Anfangsweg
t(1) = 0;               % [s] Anfangszeit
k = 1;                  % [-] Zähler
a = 0;                  % [-] Hilfsgröße zur Schnittpunktbestimmung Traktionsgrenze/reale Zugkrafthyperbel

%% Start der Rechnung
while s(k) < 75 % [m] Abbruchbedingung, Acceleration Distanz
    
    % Motordrehzahl
    ni(k) = v_V(k)*30/pi*i_ges/r_dyn(k);  % [1/min]
    
    % Motormoment für beide Motoren interpoliert
    Mi(k) = 2*interp1(n,M,ni(k),'linear','extrap'); % [Nm]
    
    % Aero-Abtriebskraft interpoliert
    % F_aero(k) = interp1(v,FA,v_V(k)*3.6,'linear','extrap'); % [N]
    F_aero(k) =  2.96 * 1.086 * 1.1643 * v_V(k)^2 / 2;
    
    % Reale Zugkrafthyperbel
    F_VXre(k) = P_Mmax/v_V(k);   % [N]
    
    % Motorleistung & Begrenzung auf 80 kW aus FS-Rules
    P_M(k) = Mi(k)*ni(k)/60*2*pi;     % [W] Gesamt-Motorleistung
    if P_M(k) > 80000
       P_M(k) = 80000;                % [W] Begrenzte Leistung
       Mi(k) = P_M(k)*60/ni(k)/2/pi;  % [Nm] Begrenzen des Moments
    end
    
    % Zugkraft an der Hinterachse
    F_VX(k)= Mi(k)*i_ges/r_dyn(k);  % [N]
    
    % Begrenzen der Zugkraft auf die Traktionsgrenze
    if F_VX(k) >= Fx_Achse(k) 
        F_VX(k) = Fx_Achse(k);
    end
        
    % Fahrwiderstände
    F_R = k_R*m_ges*g;                  % [N] Rollwiderstand
    F_L(k) = rho_L*v_V(k)^2/2*c_w*A_S;  % [N] Luftwiderstand
    F_dr(k) = F_R+F_L(k);               % [N] Gesamtwiderstand 
            
    % Bewegungsgleichungen des Fahrzeugs
    a_X(k) =(F_VX(k)-F_dr(k))/(m_ges);   % [m/s²] Beschleunigung
    v_V(k+1)= v_V(k)+a_X(k)*delta_t;     % [m/s] Geschwindigkeit
    s(k+1)= s(k)+v_V(k)*delta_t;         % [m] Weg
    
    % Dynamische Radlastverlagerung
    deltaF_Rad(k) = (m_ges*a_X(k)*h_COG)/l/2;  % [N] 
    
    % Radlast eines Hinterreifens
    Fz_r(k+1) = Fz_r(1) + deltaF_Rad(k)+F_aero(k)/2*(l-(l-(x_COP-x_vA)))/l; % [N]
       
    % Reifenmodell für übertragbare Achskraft hinten
    Fx_Achse(k+1) = 2*max(MF52_Fx_ps(ALPHA,Fz_r(k+1),GAMMA,KAPPA,TIRparam)); % [N] 

    % Schnittpunkt Traktionsgrenze/reale Zugkraft für Plot
    if a == 0 && k > 5 && Fx_Achse(k+1) >= F_VXre (k) 
        a = k;
    end
   
    % Dynamischer Reifenradius
    r_dyn(k+1) = R0 - Fz_r(k+1)/c_tire;  % [m]
    
    k = k + 1;          % Erhöhen des Zählers 
    t = t + delta_t;    % [s] Zeitintervall vergrößern
    
end

F_VXmax = max(F_VX); % [N] Maximale Zugkraft

%% Vektoren zum Plotten erzeugen
t = 0:delta_t:(k-2)*delta_t;
Fx_Achse = Fx_Achse(1:1:end-1);
Fz_r = Fz_r(1:1:end-1);
r_dyn = r_dyn(1:1:end-1);
s = s(1:1:end-1);
v_V = v_V(1:1:end-1);

%% Plotten der Ergebnisse
% Geschwindigkeitsverlauf
figure(1) 
plot(s,v_V*3.6,'b','LineWidth',1.5) 
title('Geschwindigkeitsverlauf','FontSize',12) 
xlabel('Strecke [m]','FontSize',10) 
ylabel('Geschwindigkeit [km/h]','FontSize',10)
ylim([0 140])
xlim([0 75])
grid on
box on

% Zeitverlauf
figure(2) 
plot(s,t,'b','LineWidth',1.5) 
title('Zeit','FontSize',12) 
xlabel('Strecke [m]','FontSize',10) 
ylabel('Zeit [s]','FontSize',10)
ylim([0 6])
xlim([0 75])
grid on
box on

% Drehzahlverlauf
figure(3)
plot(s,ni,'b','LineWidth',1.5) 
title('Motordrehzahl','FontSize',12) 
xlabel('Strecke [m]','FontSize',10) 
ylabel('Drehzahl [min^-^1]','FontSize',10)
ylim([0 6000])
xlim([0 75])
grid on 
box on

% Zugkraftdiagramm
figure(4)
hold on
plot(v_V*3.6,F_VX,'LineWidth',1.5) 
plot(v_V*3.6,F_dr,'LineWidth',1.5)  
plot(v_V(1:a)*3.6,Fx_Achse(1:a),'LineWidth',1.5) 
plot(v_V(3:end)*3.6,F_VXre(3:end),'LineWidth',1.5)
hold off
title('Zugkraftdiagramm','FontSize',12)
xlabel('Geschwindigkeit [km/h]','FontSize',10)
ylabel('Zugkraft [N]','FontSize',10)
legend('Zugkraft','Fahrwiderstand',...
   'Traktionsgrenze','Zugkrafthyperbel','Location','best','FontSize',15)
ylim([0 4500])
xlim([0 v_V(end)*3.6]) 
grid on
box on

% Maximale Beschleunigungen 
figure(5)
plot(v_V*3.6,a_X,'b','LineWidth',1.5) 
title('Maximale Beschleunigung','FontSize',12)
xlabel('Geschwindigkeit [km/h]','FontSize',10)
ylabel('Beschleunigung [m/s²]','FontSize',10)
xlim([0 140])
ylim([0 20]) 
grid on
box on

% Getriebeplan 
figure(6)
plot(v_V*3.6,ni)
title('Drehzahl und Geschwindigkeit','FontSize',12)
xlabel('Geschwindigkeit [km/h]','FontSize',10)
ylabel('Drehzahl [1/min]','FontSize',10)
xlim([0 140])
ylim([0 6000]) 
grid on
box on

% Leistungsverlauf
figure(7)
plot(s,P_M/1000)
title('Gesamt-Motorleistung','FontSize',12)
xlabel('Weg [m]','FontSize',10)
ylabel('Motorleistung [kW]','FontSize',10)
xlim([0 75])
ylim([0 100]) 
grid on
box on 

% Momentenverlauf am Rad
figure(8)
plot(s,Mi/2*i_ges)
title('Gesamt-Radmoment','FontSize',12)
xlabel('Weg [m]','FontSize',10)
ylabel('Moment [Nm]','FontSize',10)
xlim([0 75])
ylim([0 500]) 
grid on
box on

% Verlauf der Radlast hinten
figure(9)
plot(s,Fz_r)
title('Radlastverlauf hinten','FontSize',12)
xlabel('Weg [m]','FontSize',10)
ylabel('Radlast [N]','FontSize',10)
xlim([0 75])
ylim([0 2000]) 
grid on
box on

% Ausgabe der Werte
disp(['Die Acceleration Zeit für 75 m beträgt ' num2str(t(end)) ' s'])
disp(['Die Endgeschwindigkeit beträgt ' num2str(v_V(end)) ' m/s = ' num2str(v_V(end)*3.6) ' km/h'])
disp(['Die Motordrehzahl nach 75 m beträgt ' num2str(ni(end)) ' 1/min'])