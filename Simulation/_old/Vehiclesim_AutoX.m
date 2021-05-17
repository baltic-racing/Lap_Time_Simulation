%% Fahrzeugsimulation
clear all; clc; close all;
tic
addpath(genpath('D:\Users\Sven Weishaupt\Desktop\Vehicle Simulation'))

%% Noch zu bearbeiten
% - Schräglaufwiderstand als Fahrwiderstand einbauen (Handbuch
% Rennwagentechnik, Band Antrieb -> Interpolieren)
% - Giergeschwindigkeit + Massenträgheitsmoment des Fahrzeugs rein:
% Gierbeschleunigung und -geschwindigkeit ist sehr gering, da die
% Querkräfte zeitgleich aufgebaut werden!
% - Plots in anderen Code schreiben + GUI erstellen (APP schon erstellt)
% - Brake Bias einbauen

%% Matfiles einlesen
load('Dynojet Momentenverlauf.mat');        % Momentenverlauf von dynamischem Rollenprüfstand
load('Aero Downforce Daten.mat');           % Aero Daten
load('AutoXTrack.mat');                     % Streckendaten
load('Schräglaufwiderstand.mat');           % Schräglaufwiderstandsbeiwerte
 
x_Track = Track(:,1);           % [m] X-Koordinate der Strecke
y_Track = Track(:,2);           % [m] Y-Koordinate der Strecke
s = Track(:,3);                 % [m] Verlauf der Streckenlänge
R = Track(:,4);                 % [m] Kurvenradien
ApexIndizes = Apexes(abs(R));   % [-] Indizes der Apexes

%% Fahrzeugdaten
m_ges = 275;    % [kg] Fahrzeuggesamtmasse inkl. Fahrer
g = 9.81;       % [m/s²] Erdbeschleunigung
FG = m_ges*g;   % [N] Gewichtskraft des Fahrzeugs

l = 1525;       % [mm] Radstand
B = 1500;       % [mm] Spurweite für vorne und hinten
h_COG = 246.37; % [mm] Höhe Fahrzeugschwerpunkt
% h_ROLL = 150;   % [mm] Höhe Rollzentrum
% h_NICK = 200;   % [mm] Höhe Nickzentrum
x_COG = 1400;   % [mm] x-Koordinate des Fahrzeugschwerpunkts in CAD
x_vA = 470;     % [mm] x-Koordinate der Vorderachse in CAD
x_COP = 1600;   % [mm] x-Koordinate Center of Pressure in CAD
aero_ph = (l-(l-(x_COP-x_vA)))/l;   % [-] Anteil Aerokraft auf Hinterachse
aero_pv = 1-aero_ph;                % [-] Anteil Aerokraft auf Vorderachse

thetaV_X = 29.3828;  % [kg*m²] Trägheitsmoment Fahrzeug um X-Achse
thetaV_Y = 104.837;  % [kg*m²] Trägheitsmoment Fahrzeug um Y-Achse
thetaV_Z = 119.2005; % [kg*m²] Trägheitsmoment Fahrzeug um Z-Achse

m_ph = 66;              % [%] Prozentualer Radlastanteil Hinterachse
lv = l*m_ph/100;        % [mm] Abstand Vorderachse zu Fahrzeugschwerpunkt
lh = l-lv;              % [mm] Abstand Hinterachse zu Fahrzeugschwerpunkt

FB = 7000;      % [N] Maximale Bremskraft

k_R = 0.01;     % [-] Rollwiderstandsbeiwert
c_w = 1.463;    % [-] cw-Wert
A_S = 1.3;      % [m²] Stirnfläche

% Umgebungsbedingungen
t_L = 25;       % [°C] Umgebungslufttemperatur
p_L = 1.013;    % [bar] Umgebungsluftdruck
R_L = 287;      % [J/(kg*K)] Gaskonstante Luft
rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m³] Luftdichte

% Motordaten
P = M.*n*2*pi/60;   % [W] Leistungsmatrix einlesen
P_Mmax = max(P);    % [W] Ermittlung der Maximalleistung
n_Mmax = 10500;     % [1/min] Maximaldrehzahl des Motors
n_shift = 10500;    % [1/min] Motordrehzahl beim Schalten
eta_ges = 0.90;     % [-] Gesamtwirkungsgrad des Antriebsstrangs
t_shift = 0.035;    % [s] Schaltzeit (Zugkraftunterbrechung)

% Getriebedaten
i_Primaer = 1.848;      % [-] Primärübersetzung
z_Kettenblatt = 37;     % [-] Zähnezahl des Kettenblatts
z_Ritzel = 11;          % [-] Zähnezahl des Ritzels
i_Sekundaer = z_Kettenblatt / z_Ritzel; % [-] Sekundärübersetzung

i = [2.579 1.720 1.429 1.238]; % [-] Übersetzungen der einzelnen Gänge
% i = [1.857 1.565 1.35 1.238 1.136]; % [-] Übersetzungen der einzelnen Gänge
ig = i.*i_Primaer*i_Sekundaer;      % [-] Gesamtübersetzungen
Gang_max = length(i);               % [-] Höchster Gang

%% Initialisieren aller Variablen
alpha_v = zeros(1,length(Track)-1);
alpha_vl = zeros(1,length(Track)-1);
alpha_vr = zeros(1,length(Track)-1);
alpha_h = zeros(1,length(Track)-1);
alpha_hl = zeros(1,length(Track)-1);
alpha_hr = zeros(1,length(Track)-1);
aRev = zeros(1,length(Track)-1);
aVX = zeros(1,length(Track)-1);
aVY = zeros(1,length(Track)-1);
beta = zeros(1,length(Track)-1);
BPPsignal = zeros(1,length(Track)-1);
cY_vl = zeros(1,length(Track));
cY_vr = zeros(1,length(Track));
cY_hl = zeros(1,length(Track));
cY_hr = zeros(1,length(Track));
cZ_vl = zeros(1,length(Track));
cZ_vr = zeros(1,length(Track));
cZ_hl = zeros(1,length(Track));
cZ_hr = zeros(1,length(Track));
delta = zeros(1,length(Track)-1);
dFWZhl_aero = zeros(1,length(Track)-1);
dFWZhr_aero = zeros(1,length(Track)-1);
dFWZvl_aero = zeros(1,length(Track)-1);
dFWZvr_aero = zeros(1,length(Track)-1);
dFWZhl_x = zeros(1,length(Track)-1);
dFWZhr_x = zeros(1,length(Track)-1);
dFWZvl_x = zeros(1,length(Track)-1);
dFWZvr_x = zeros(1,length(Track)-1);
dFWZhl_y = zeros(1,length(Track)-1);
dFWZhr_y = zeros(1,length(Track)-1);
dFWZvl_y = zeros(1,length(Track)-1);
dFWZvr_y = zeros(1,length(Track)-1);
Gangauswahl = zeros(1,length(Track)-1);
Faero = zeros(1,length(Track)-1);
FB = FB*ones(1,length(Track));
FVX = zeros(1,length(Track)-1);
FVX_hl = zeros(1,length(Track)-1);
FVX_hr = zeros(1,length(Track)-1);
FVXre = zeros(1,length(Track)-1);
FVXid = zeros(1,length(Track)-1);
FR = zeros(1,length(Track)-1);
FL = zeros(1,length(Track)-1);
Fdr = zeros(1,length(Track)-1);
FVY = zeros(1,length(Track)-1);
FWXmax_v = zeros(1,length(Track));
FWXmax_h = zeros(1,length(Track));
FWXmax_vl = zeros(1,length(Track));
FWXmax_vr = zeros(1,length(Track));
FWXmax_hl = zeros(1,length(Track));
FWXmax_hr = zeros(1,length(Track));
FWYv = zeros(1,length(Track)-1);
FWYh = zeros(1,length(Track)-1);
FWYmax_v = zeros(1,length(Track));
FWYmax_h = zeros(1,length(Track));
FWYmax_vl = zeros(1,length(Track));
FWYmax_vr = zeros(1,length(Track));
FWYmax_hl = zeros(1,length(Track));
FWYmax_hr = zeros(1,length(Track));
FWZges = zeros(1,length(Track));
FWZ_hl = zeros(1,length(Track));
FWZ_hr = zeros(1,length(Track));
FWZ_vl = zeros(1,length(Track));
FWZ_vr = zeros(1,length(Track));
FWZh = zeros(1,length(Track)-1);
FWZv = zeros(1,length(Track)-1);
Mi = zeros(1,length(Track)-1);
ni = zeros(1,length(Track)-1);
psi1 = zeros(1,length(Track)-1);
Rdyn_vl = zeros(1,length(Track));
Rdyn_vr = zeros(1,length(Track));
Rdyn_hl = zeros(1,length(Track));
Rdyn_hr = zeros(1,length(Track));
RutschenY_v = zeros(1,length(Track));
RutschenY_h = zeros(1,length(Track));
t = zeros(1,length(Track)-1);
TC = zeros(1,length(Track));
Tirelimit = zeros(1,length(Track));
vAPEXmax = zeros(1,length(ApexIndizes));
vV = zeros(1,length(Track)-1);
vRev = zeros(1,length(Track)-1);
vVYmax = zeros(1,length(Track)-1);

%% Reifenmodell - Magic Tire Formula Model 5.2
load('VerticalStiffness_65kPA_IA0.mat');    % Lookup-Table für vertikale Reifensteifigkeit
load('CorneringStiffness_65kPA_IA0.mat');   % Lookup-Table für Schräglaufseitensteifigkeit

% Variable für Pfad des Tir-Files
FileNameLocation = ('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

% Tir-File in Struktur laden
TIRparam = loadTIR(FileNameLocation);

% Reifendaten
J_Tire = 0.2;                  % [kg*m²] Trägheitsmoment des Reifens um Drehachse
R0 = TIRparam.UNLOADED_RADIUS; % [m] Fertigungsradius des Reifens
bW = TIRparam.WIDTH;           % [m] Breite des Reifens (Aufstandsbreite)
p_infl = 65000;                % [Pa] Luftdruck des Reifens
GAMMA = 0;                     % [°] Sturz

% Skalierungsfaktoren für Grip-Niveau
% 0.75 für optimale Reifentemperatur, 0.6 für niedrige Reifentemperatur
TIRparam.LMUX = 0.6;          % [-] Skalierungsfaktor Längsrichtung
TIRparam.LMUY = 0.6;          % [-] Skalierungsfaktor Querrichtung

% Zusätzliche Faktoren für Wechselwirkung Längs/Querschlupf
TIRparam.LXAL = 1.2;    % [-] Einfluss Schräglaufwinkel auf übertragbare Längskraft
TIRparam.LYKA = 1.2;    % [-] Einfluss Längsschlupf auf übertragbare Querkraft
TIRparam.RBX3 = 0;      % [-] Zusätzlicher Faktor für combined slip Fx

% (Statische) Achs- und Radlasten
FWZges(1) = FG;                 % [N] Statische Gesamtachslast 
FWZh(1) = m_ph/100*FWZges(1);   % [N] Statische Achslast hinten
FWZv(1) = FWZges(1)-FWZh(1);    % [N] Statische Achslast vorne
FWZ_vr(1) = FWZv(1)/2;          % [N] Statische Radlast vorne rechts  
FWZ_vl(1) = FWZv(1)/2;          % [N] Statische Radlast vorne links
FWZ_hr(1) = FWZh(1)/2;          % [N] Statische Radlast hinten rechts  
FWZ_hl(1) = FWZh(1)/2;          % [N] Statische Radlast hinten links  

% Vertikale Reifensteifigkeiten interpoliert
cZ_vl(1) = interp1(Fz,cZ_tire,FWZ_vl(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
cZ_vr(1) = interp1(Fz,cZ_tire,FWZ_vr(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
cZ_hl(1) = interp1(Fz,cZ_tire,FWZ_hl(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
cZ_hr(1) = interp1(Fz,cZ_tire,FWZ_hr(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert

% Schräglaufseitensteifigkeiten interpoliert
cY_vl(1) = interp1(Fzy,cY_tire,FWZ_vl(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
cY_vr(1) = interp1(Fzy,cY_tire,FWZ_vr(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
cY_hl(1) = interp1(Fzy,cY_tire,FWZ_hl(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
cY_hr(1) = interp1(Fzy,cY_tire,FWZ_hr(1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert

% Dynamische Reifenradien (im Stand = statisch)
Rdyn_vl(1) = R0 - FWZ_vl(1)/cZ_vl(1);    % [m] Dynamischer Reifenradius 
Rdyn_vr(1) = R0 - FWZ_vr(1)/cZ_vr(1);    % [m] Dynamischer Reifenradius 
Rdyn_hl(1) = R0 - FWZ_hl(1)/cZ_hl(1);    % [m] Dynamischer Reifenradius 
Rdyn_hr(1) = R0 - FWZ_hr(1)/cZ_hr(1);    % [m] Dynamischer Reifenradius  

FWXmax_h(1) = Inf;

%% Berechnen der maximalen Kurvengeschwindigkeiten für alle Apexes (numerisch)

% Statische Reifenlasten für Berechnung
FWZ_vl_stat = FWZ_vl(1);
FWZ_vr_stat = FWZ_vr(1);
FWZ_hl_stat = FWZ_hl(1);
FWZ_hr_stat = FWZ_hr(1);

for i = 1:length(ApexIndizes)
 
    FWYv(i) = 0;            % [N] Startwert Querkraft Vorderachse
    FWYh(i) = 0;            % [N] Startwert Querkraft Hinterachse
    FWYmax_v(i) = 0.1;      % [N] Startwert maximal übertragbare Querkraft Vorderachse
    FWYmax_h(i) = 0.1;      % [N] Startwert maximal übertragbare Querkraft Hinterachse
    vV(i) = 0;              % [m/s] Startwert Fahrzeuggeschwindigkeit
    
    if R(ApexIndizes(i)) > 0
        f = -1;
    else
        f = 1;
    end
    
    while  FWYv(i) < FWYmax_v(i) && FWYh(i) < FWYmax_h(i) && vV(i) < 30
        
        vV(i) = vV(i) + 0.01;   % [m/s] Erhöhen der Fahrzeuggeschwindigkeit
        
        Faero(i) = interp1(v,FA,vV(i)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert    
        
        FVY(i) = m_ges*vV(i)^2/R(ApexIndizes(i));    % [N] Zentrifugalkraft
        
        aVY(i) = vV(i)^2/R(i);  % [m/s²] Querbeschleunigung
        
        % Aufzubringende Querkräfte an Vorder- und Hinterachse
        FWYv(i) = lh/l*abs(FVY(i));   % [N] Aufzubringende Querkraft der Vorderachse
        FWYh(i) = lv/l*abs(FVY(i));   % [N] Aufzubringende Querkraft der Hinterachse
        
        cY_v(i) = interp1(Fzy,cY_tire,(FWZ_vl(i)+FWZ_vr(i))/2,'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
        cY_h(i) = interp1(Fzy,cY_tire,(FWZ_hl(i)+FWZ_hr(i))/2,'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
        
        Rh(i) = sqrt(R(ApexIndizes(i))^2-(lh/1000)^2);
        
        delta_a(i) = atan(l/1000/Rh(i));    % [rad] Ackermannwinkel seitenkraftfreie Kurvenfahrt
        
        % Eigenlenkgradient
        EG(i) = (m_ges*(cY_h(i)*180/pi*lh-cY_v(i)*180/pi*lv))/(l*cY_v(i)*180/pi*cY_h(i)*180/pi);
        
        DELTA(i) = delta_a(i)+EG(i)*abs(aVY(i));
        
        % Lenkwinkel, Schwimmwinkel, Gierrate, Schräglaufwinkel  
        delta(i) = f*atan((l/1000)/sqrt(R(ApexIndizes(i))^2-(lh/1000)^2));   % [rad] Lenkwinkel
        deltam(i) = atan(l/1000/R(ApexIndizes(i)));
        beta(i) = f*atan((lh/1000)/sqrt(R(ApexIndizes(i))^2-(lh/1000)^2));  % [rad] Schwimmwinkel
% beta(i) = 0;
        psi1(i) = vV(i)/R(ApexIndizes(i));                               % [rad/s] Gierrate
        alpha_v(i) = 180/pi*(delta(i)-(lv/1000)/vV(i)*psi1(i)-beta(i));  % [°] Schräglaufwinkel vorne
        alpha_h(i) = 180/pi*((lh/1000)/vV(i)*psi1(i)-beta(i));           % [°] Schräglaufwinkel hinten
        
        alphaVR(i) = delta(i)-atan((vV(i)*sin(beta(i))+lv/1000*psi1(i))/...
            (vV(i)*cos(beta(i))+0.5*B/1000*psi1(i)));
        alphaVL(i) = delta(i)-atan((vV(i)*sin(beta(i))+lv/1000*psi1(i))/...
            (vV(i)*cos(beta(i))-0.5*B/1000*psi1(i)));
        alphaHR(i) = -atan((vV(i)*sin(beta(i))+lh/1000*psi1(i))/...
            (vV(i)*cos(beta(i))+0.5*B/1000*psi1(i)));
        alphaHL(i) = -atan((vV(i)*sin(beta(i))-lh/1000*psi1(i))/...
            (vV(i)*cos(beta(i))-0.5*B/1000*psi1(i)));
        
        % Radlastverlagerung in Folge von Aerokräften
        dFWZhl_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf linkes Hinterrad
        dFWZhr_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf rechtes Hinterrad
        dFWZvl_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf linkes Vorderrad
        dFWZvr_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf rechtes Vorderrad
        
        % Dynamische Radlastverlagerung in Längsrichtung = 0 angenommen
        dFWZhl_x(i) = 0;   % [N] Dynamische Radlastverlagerung linkes Hinterrad
        dFWZhr_x(i) = 0;   % [N] Dynamische Radlastverlagerung rechtes Hinterrad
        dFWZvl_x(i) = 0;   % [N] Dynamische Radlastverlagerung linkes Vorderrad
        dFWZvr_x(i) = 0;   % [N] Dynamische Radlastverlagerung rechtes Vorderrad
        
        % Dynamische Radlastverlagerung in Querrichtung
        dFWZvl_y(i) = -h_COG/B*lh/l*FVY(i);  % [N] Dynamische Radlastverlagerung linkes Vorderrad
        dFWZvr_y(i) = h_COG/B*lh/l*FVY(i);  % [N] Dynamische Radlastverlagerung rechtes Vorderrad
        dFWZhl_y(i) = -h_COG/B*lv/l*FVY(i);  % [N] Dynamische Radlastverlagerung linkes Hinterrad
        dFWZhr_y(i) = h_COG/B*lv/l*FVY(i);  % [N] Dynamische Radlastverlagerung rechtes Hinterrad
        
        % Radlasten
        FWZ_hl(i) = FWZ_hl_stat + dFWZhl_aero(i) + dFWZhl_x(i) + dFWZhl_y(i); % [N] Radlast hinten links
        FWZ_hr(i) = FWZ_hr_stat + dFWZhr_aero(i) + dFWZhr_x(i) + dFWZhr_y(i); % [N] Radlast hinten rechts
        FWZ_vl(i) = FWZ_vl_stat + dFWZvl_aero(i) + dFWZvl_x(i) + dFWZvl_y(i); % [N] Radlast vorne links
        FWZ_vr(i) = FWZ_vr_stat + dFWZvr_aero(i) + dFWZvr_x(i) + dFWZvr_y(i); % [N] Radlast vorne rechts
        
        alphaREAR(i) = MF52_ALPHA_cs(FWYh(i)/2,(FWZ_hl(i)+FWZ_hr(i))/2,0,TIRparam);
        
        % Maximal übertragbare Reifenkräfte in Querrichtung
        FWYmax_vl(i) = abs(MF52_Fy_cs(alpha_v(i),FWZ_vl(i),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
        FWYmax_vr(i) = abs(MF52_Fy_cs(alpha_v(i),FWZ_vr(i),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
        FWYmax_hl(i) = abs(MF52_Fy_cs(alpha_h(i),FWZ_hl(i),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
        FWYmax_hr(i) = abs(MF52_Fy_cs(alpha_h(i),FWZ_hr(i),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
       
        FWYmax_v(i) = FWYmax_vl(i) + FWYmax_vr(i);    % [N] Maximal übertragbare Achskraft
        FWYmax_h(i) = FWYmax_hl(i) + FWYmax_hr(i);    % [N] Maximal übertragbare Achskraft
                
    end
    
    vAPEXmax(i) = vV(i);   % [m/s] Maximalgeschwindigkeit für jede Apex
end

% figure(2)
% plot(1:65,delta(1:65),1:65,deltam,1:65,DELTA,1:65,delta_a)
% legend('alt','neu','neu2','Ackermann')
% 
% figure(3)
% plot(1:65,delta(1:65)*180/pi)

figure(1)
plot(1:65,alphaREAR*180/pi)

%% Startwerte für Simulation
Gang = 1;      % [-] Fahrzeug startet im ersten Gang
vV(1) = 0;     % [m/s] Geschwindigkeit
t(1) = 0;      % [s] Zeit

% Hilfsgrößen
t_z = zeros(1,Gang_max);    % [s] Bestimmung der Zeit pro Gang  
t_x = 0;                    % [s] Berücksichtigung der Zugkraftunterbechung während des Schaltens
z = 1;                      % [-] Bestimmung der anstehenden Apex

%% Simulation
for i = 1:length(Track)-1
    
    % Bestimmen von Motordrehzahl und Gang
    Gangauswahl(i) = Gang; % [-] Speichern des Gangverlaufs
    ni(i) = vV(i)*30/pi*ig(Gang)/Rdyn_hl(i); % [1/min] Aktuelle Motordrehzahl ermitteln
    
    if ni(i) < 6000 && Gang > 1 % Runterschalten
       Gang = Gang-1;    
       ni(i) = vV(i)*30/pi*ig(Gang)/Rdyn_hl(i); % [1/min] Drehzahl nach Gangwechsel ermitteln
    end
    
    if ni(i) >= n_shift % Gang erhöhen, wenn Schaltdrehzahl erreicht
        Gang = Gang + 1;
        t_x = t(i)-t(i-1); % Schrittweite der Hilfsgröße festlegen
            if Gang >= Gang_max  % Begrenzen der höchsten Gangzahl
                Gang = Gang_max;
            end 
        ni(i)= vV(i)*30/pi*ig(Gang)/Rdyn_hl(i); % [1/min] Drehzahl nach Gangwechsel ermitteln      
            if ni(i)>= n_Mmax % Drehzahlbegrenzer
                ni(i) = n_Mmax;
            end      
    end
    
    % Bestimmen von Aero-Kräften und Motormoment
    Faero(i) = interp1(v,FA,vV(i)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
    Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motormoment interpoliert
    
    % Berechnen der Zugkraft
    if t_x == 0
        FVX_hl(i) = Mi(i)/2*ig(Gang)/Rdyn_hl(i); % [N] Zugkraft an linkem Hinterrad
        FVX_hr(i) = Mi(i)/2*ig(Gang)/Rdyn_hr(i); % [N] Zugkraft an rechtem Hinterrad  
    elseif t_x >=0 && t_shift > t_x
        FVX_hl(i) = 0;
        FVX_hr(i) = 0;
        t_x = t_x+t(i)-t(i-1);
    elseif t_x >= t_shift
        t_x = 0;
        FVX_hl(i) = Mi(i)/2*ig(Gang)/Rdyn_hl(i); % [N] Zugkraft an linkem Hinterrad
        FVX_hr(i) = Mi(i)/2*ig(Gang)/Rdyn_hr(i); % [N] Zugkraft an rechtem Hinterrad
    end
    
    FVX(i) = FVX_hl(i) + FVX_hr(i);  % [N] Zugkraft an der Hinterachse
    
    if FVX(i) > FWXmax_h(i)   % Begrenzen der Zugkraft auf Traktionsgrenze
        FVX(i) = FWXmax_h(i);
    end
    
    % Lenkwinkel, Schwimmwinkel, Gierrate, Schräglaufwinkel
    if R(i) > 0
        f = -1;     % [-] Hilfsgröße für Richtung der Winkel
    else
        f = 1;
    end
    
    delta(i) = f*atan((l/1000)/sqrt(R(i)^2-(lh/1000)^2));   % [rad] Lenkwinkel
    beta(i) = f*atan((lh/1000)/sqrt(R(i)^2-(lh/1000)^2));   % [rad] Schwimmwinkel
    psi1(i) = vV(i)/R(i);                                   % [rad/s] Gierrate
    alpha_v(i) = 180/pi*(delta(i)-(lv/1000)/vV(i)*psi1(i)-beta(i));  % [°] Schräglaufwinkel vorne
    alpha_h(i) = 180/pi*((lh/1000)/vV(i)*psi1(i)-beta(i));           % [°] Schräglaufwinkel hinten
    
    % Fahrwiderstände
    k_av(i) = abs(interp1(ALPHA,k_alpha,alpha_v(i),'linear','extrap')); % [-] Schräglaufwiderstandsbeiwert vorne
    k_ah(i) = abs(interp1(ALPHA,k_alpha,alpha_h(i),'linear','extrap')); % [-] Schräglaufwiderstandsbeiwert hinten
    
    F_ALPHA_v(i) = k_av(i)*FWZv(i); % [N] Schräglaufwiderstand vorne
    F_ALPHA_h(i) = k_ah(i)*FWZh(i); % [N] Schräglaufwiderstand hinten
    
    FR(i) = k_R*FWZges(i);            % [N] Rollwiderstand
    FL(i) = rho_L*vV(i)^2/2*c_w*A_S;  % [N] Luftwiderstand
    Fdr(i) = FR(i)+FL(i);             % [N] Gesamtwiderstand
     
    % Fahrzeug
    FVY(i) = m_ges*vV(i)^2/R(i);      % [N] Zentrifugalkraft
    aVX(i) = (FVX(i)-Fdr(i))/m_ges;   % [m/s²] Längsbeschleunigung
    aVY(i) = vV(i)^2/R(i);            % [m/s²] Querbeschleunigung

    if ismember(i,ApexIndizes)
        if vV(i) > vAPEXmax(z)   % Begrenzen auf maximale Kurvengeschwindigkeit in Apexes
            vV(i) = vAPEXmax(z);
        end
        z = z + 1;
    end
    
    vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Gesamt-Fahrzeuggeschwindigkeit
    t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Zeit
    
    % Aufzubringende Querkräfte an Vorder- und Hinterachse
    FWYv(i) = lh/l*FVY(i);   % [N] Aufzubringende Querkraft der Vorderachse
    FWYh(i) = lv/l*FVY(i);   % [N] Aufzubringende Querkraft der Hinterachse
        
    % Radlastverlagerung in Folge von Aerokräften
    dFWZhl_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf linkes Hinterrad
    dFWZhr_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf rechtes Hinterrad
    dFWZvl_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf linkes Vorderrad
    dFWZvr_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf rechtes Vorderrad
     
    % Dynamische Radlastverlagerungen in Längsrichtung
    dFWZhl_x(i) = m_ges*aVX(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung linkes Hinterrad
    dFWZhr_x(i) = m_ges*aVX(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung rechtes Hinterrad
    dFWZvl_x(i) = -m_ges*aVX(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung linkes Vorderrad
    dFWZvr_x(i) = -m_ges*aVX(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung rechtes Vorderrad
    
    % Dynamische Radlastverlagerung in Querrichtung
    dFWZvl_y(i) = -h_COG/B*lh/l*FVY(i);  % [N] Dynamische Radlastverlagerung linkes Vorderrad
    dFWZvr_y(i) = h_COG/B*lh/l*FVY(i);   % [N] Dynamische Radlastverlagerung rechtes Vorderrad
    dFWZhl_y(i) = -h_COG/B*lv/l*FVY(i);  % [N] Dynamische Radlastverlagerung linkes Hinterrad
    dFWZhr_y(i) = h_COG/B*lv/l*FVY(i);   % [N] Dynamische Radlastverlagerung rechtes Hinterrad

    % Radlasten
    FWZ_hl(i+1) = FWZ_hl(1) + dFWZhl_aero(i) + dFWZhl_x(i) + dFWZhl_y(i); % [N] Radlast hinten links
    FWZ_hr(i+1) = FWZ_hr(1) + dFWZhr_aero(i) + dFWZhr_x(i) + dFWZhr_y(i); % [N] Radlast hinten rechts
    FWZ_vl(i+1) = FWZ_vl(1) + dFWZvl_aero(i) + dFWZvl_x(i) + dFWZvl_y(i); % [N] Radlast vorne links
    FWZ_vr(i+1) = FWZ_vr(1) + dFWZvr_aero(i) + dFWZvr_x(i) + dFWZvr_y(i); % [N] Radlast vorne rechts
    
    % Begrenzen der Radlasten auf (quasi) Null
    if FWZ_vl(i+1) < 0
        FWZ_vl(i+1) = 0.001;
    end
    if FWZ_vr(i+1) < 0
        FWZ_vr(i+1) = 0.001;
    end
    if FWZ_hl(i+1) < 0
        FWZ_hl(i+1) = 0.001;
    end
    if FWZ_hr(i+1) < 0
        FWZ_hr(i+1) = 0.001;
    end
    
    % Achslasten
    FWZh(i+1) = FWZ_hl(i+1) + FWZ_hr(i+1);  % [N] Achslast hinten
    FWZv(i+1) = FWZ_vl(i+1) + FWZ_vr(i+1);  % [N] Achslast vorne
    FWZges(i+1) = FWZh(i+1) + FWZv(i+1);    % [N] Gesamtachs- und Radlast
    
    % Vertikale Reifensteifigkeiten
    cZ_vl(i+1) = interp1(Fz,cZ_tire,FWZ_vl(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    cZ_vr(i+1) = interp1(Fz,cZ_tire,FWZ_vr(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    cZ_hl(i+1) = interp1(Fz,cZ_tire,FWZ_hl(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    cZ_hr(i+1) = interp1(Fz,cZ_tire,FWZ_hr(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    
    % Dynamische Reifenradien
    Rdyn_vl(i+1) = R0 - FWZ_vl(i+1)/cZ_vl(i+1);    % [m] Dynamischer Reifenradius
    Rdyn_vr(i+1) = R0 - FWZ_vr(i+1)/cZ_vr(i+1);    % [m] Dynamischer Reifenradius
    Rdyn_hl(i+1) = R0 - FWZ_hl(i+1)/cZ_hl(i+1);    % [m] Dynamischer Reifenradius
    Rdyn_hr(i+1) = R0 - FWZ_hr(i+1)/cZ_hr(i+1);    % [m] Dynamischer Reifenradius
    
    % Maximal übertragbare Reifenkräfte in Längsrichtung
    FWXmax_vl(i+1) = max(abs(MF52_Fx_cs(alpha_v(i),FWZ_vl(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_vr(i+1) = max(abs(MF52_Fx_cs(alpha_v(i),FWZ_vr(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_hl(i+1) = max(abs(MF52_Fx_cs(alpha_h(i),FWZ_hl(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_hr(i+1) = max(abs(MF52_Fx_cs(alpha_h(i),FWZ_hr(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_v(i+1) = FWXmax_vl(i+1) + FWXmax_vr(i+1);    % [N] Maximal übertragbare Achskraft
    FWXmax_h(i+1) = FWXmax_hl(i+1) + FWXmax_hr(i+1);    % [N] Maximal übertragbare Achskraft
    
    % Maximal übertragbare Reifenkräfte in Querrichtung
    FWYmax_vl(i+1) = abs(MF52_Fy_cs(alpha_v(i),FWZ_vl(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_vr(i+1) = abs(MF52_Fy_cs(alpha_v(i),FWZ_vr(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_hl(i+1) = abs(MF52_Fy_cs(alpha_h(i),FWZ_hl(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_hr(i+1) = abs(MF52_Fy_cs(alpha_h(i),FWZ_hr(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_v(i+1) = FWYmax_vl(i+1) + FWYmax_vr(i+1);    % [N] Maximal übertragbare Achskraft
    FWYmax_h(i+1) = FWYmax_hl(i+1) + FWYmax_hr(i+1);    % [N] Maximal übertragbare Achskraft
    
    t_z(Gang) = t_z(Gang) + t(i+1) - t(i); % [s] Zeit pro Gang
end

% Zwischenspeichern von diversen Variablen für Plots
vV_NoBrake = vV;    
Gangauswahl_NoBrake = Gangauswahl;
ApexGang_Soll = Gangauswahl(ApexIndizes);

%% BREMSPUNKTBERECHNUNG
BrakeIndizes = [];

for k = 1:length(ApexIndizes)
    
    Counter = 0;
    j = ApexIndizes(k);
    vRev(j-1) = vAPEXmax(k);
    
    while vRev(j-1) < vV(j-1)
        Faero(j-1) = interp1(v,FA,vRev(j-1)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
        FR(j-1) = k_R*(FG+Faero(j-1));
        FL(j-1) = rho_L*vRev(j-1)^2/2*c_w*A_S;
        Fdr(j-1) = FR(j-1)+FL(j-1);
        aRev(j-1) = (-Fdr(j-1)-FB(j-1))/m_ges;
        vRev(j-2) = sqrt(vRev(j-1)^2-2*aRev(j-1)*(s(j-1)-s(j-2)));
        j = j - 1;
        Counter = Counter + 1;
    end
   
    if Counter > 0
        BrakeIndizes = [BrakeIndizes j-1:ApexIndizes(k)-1];
    end
    
end

%% Startwerte für Simulation MIT BREMSEN
Gang = 1;      % [-] Fahrzeug startet im ersten Gang
vV(1) = 0;     % [m/s] Geschwindigkeit
t(1) = 0;      % [s] Zeit

% Hilfsgrößen
t_z = zeros(1,Gang_max);    % [s] Bestimmung der Zeit pro Gang  
t_x = 0;                    % [s] Berücksichtigung der Zugkraftunterbechung während des Schaltens
z = 1;                      % [-] Bestimmung der anstehenden Apex

%% SIMULATION MIT BREMSEN
for i = 1:length(Track)-1
    
    % Überprüfen, vor welcher Apex das Auto ist
    if ismember(i,ApexIndizes)  
        z = z + 1;
    end
    
    % Bestimmen von Motordrehzahl und Gang
    Gangauswahl(i) = Gang; % [-] Speichern des Gangverlaufs
    
    ni(i) = vV(i)*30/pi*ig(Gang)/Rdyn_hl(i); % [1/min] Aktuelle Motordrehzahl ermitteln
    
    if ni(i) < 6000 && Gang > 1 % Runterschalten
       Gang = Gang-1;    
       ni(i) = vV(i)*30/pi*ig(Gang)/Rdyn_hl(i); % [1/min] Drehzahl nach Gangwechsel ermitteln
    end
    
    if ni(i) >= n_shift % Gang erhöhen, wenn Schaltdrehzahl erreicht
        Gang = Gang + 1;
        t_x = t(i)-t(i-1); % Schrittweite der Hilfsgröße festlegen
            if Gang >= Gang_max  % Begrenzen der höchsten Gangzahl
                Gang = Gang_max;
            end 
        ni(i)= vV(i)*30/pi*ig(Gang)/Rdyn_hl(i); % [1/min] Drehzahl nach Gangwechsel ermitteln      
            if ni(i)>= n_Mmax % Drehzahlbegrenzer
                ni(i) = n_Mmax;
            end      
    end
    
    % Bestimmen von Aero-Kräften und Motormoment
    Faero(i) = interp1(v,FA,vV(i)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
    Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motormoment interpoliert
    
    % Berechnen der Zugkraft
    if t_x == 0
        FVX_hl(i) = Mi(i)/2*ig(Gang)/Rdyn_hl(i); % [N] Zugkraft an linkem Hinterrad
        FVX_hr(i) = Mi(i)/2*ig(Gang)/Rdyn_hr(i); % [N] Zugkraft an rechtem Hinterrad  
    elseif t_x >=0 && t_shift > t_x
        FVX_hl(i) = 0;
        FVX_hr(i) = 0;
        t_x = t_x+t(i)-t(i-1);
    elseif t_x >= t_shift
        t_x = 0;
        FVX_hl(i) = Mi(i)/2*ig(Gang)/Rdyn_hl(i); % [N] Zugkraft an linkem Hinterrad
        FVX_hr(i) = Mi(i)/2*ig(Gang)/Rdyn_hr(i); % [N] Zugkraft an rechtem Hinterrad
    end
    
    % Prüfen, ob gebremst werden muss
    if ismember(i,BrakeIndizes) % Einleiten des Bremsvorgangs     
        BPPsignal(i) = 1;   %  Bremssignal
        FVX_hl(i) = 0;      % [N] Zugkraft an linkem Hinterrad
        FVX_hr(i) = 0;      % [N] Zugkraft an rechtem Hinterrad
    else
        FB(i) = 0;          % [N] Bremskraft
    end
        
    FVX(i) = FVX_hl(i) + FVX_hr(i);           % [N] Zugkraft an der Hinterachse
    
    if FVX(i) > FWXmax_h(i)     % Begrenzen der Zugkraft auf Traktionsgrenze
        FVX(i) = FWXmax_h(i);
        TC(i) = 1;              % Traktionskontrolle "an"
    end
    
    % Fahrwiderstände
    FR(i) = k_R*FWZges(i);            % [N] Rollwiderstand
    FL(i) = rho_L*vV(i)^2/2*c_w*A_S;  % [N] Luftwiderstand
    Fdr(i) = FR(i)+FL(i);             % [N] Gesamtwiderstand
    
    % Fahrzeug
    FVY(i) = m_ges*vV(i)^2/R(i);                    % [N] Zentrifugalkraft 
    aVX(i) = (FVX(i)-Fdr(i)-FB(i))/(m_ges);         % [m/s²] Longitudinal Beschleunigung
    aVY(i) = vV(i)^2/R(i);                          % [m/s²] Querbeschleunigung
    vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Gesamt-Fahrzeuggeschwindigkeit
    
    if ismember(i,BrakeIndizes) && vV(i+1) < vAPEXmax(z) % Begrenzen der Geschwindigkeit auf ApexGeschwindigkeit
        vV(i+1) = vAPEXmax(z);   % [m/s] Gesamt-Fahrzeuggeschwindigkeit
    end
    
    t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Zeit
    
    % Querkräfte an Vorder- und Hinterachse
    FWYv(i) = lh/l*FVY(i);   % [N] Aufzubringende Querkraft der Vorderachse
    FWYh(i) = lv/l*FVY(i);   % [N] Aufzubringende Querkraft der Hinterachse
    
    if FWYv(i) > FWYmax_v(i)
       RutschenY_v(i) = 1;  
    end
    
     if FWYh(i) > FWYmax_h(i)
       RutschenY_h(i) = 1;  
    end
    
     % Radlastverlagerung in Folge von Aerokräften 
     dFWZhl_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf linkes Hinterrad
     dFWZhr_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf rechtes Hinterrad
     dFWZvl_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf linkes Vorderrad
     dFWZvr_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf rechtes Vorderrad
     
     % Dynamische Radlastverlagerungen in Längsrichtung
     dFWZvl_x(i) = -m_ges*aVX(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung linkes Vorderrad
     dFWZvr_x(i) = -m_ges*aVX(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung rechtes Vorderrad
     dFWZhl_x(i) = m_ges*aVX(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung linkes Hinterrad
     dFWZhr_x(i) = m_ges*aVX(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung rechtes Hinterrad
     
     % Dynamische Radlastverlagerung in Querrichtung
     dFWZvl_y(i) = -h_COG/B*lh/l*FVY(i);   % [N] Dynamische Radlastverlagerung linkes Vorderrad
     dFWZvr_y(i) = h_COG/B*lh/l*FVY(i);    % [N] Dynamische Radlastverlagerung rechtes Vorderrad
     dFWZhl_y(i) = -h_COG/B*lv/l*FVY(i);   % [N] Dynamische Radlastverlagerung linkes Hinterrad
     dFWZhr_y(i) = h_COG/B*lv/l*FVY(i);    % [N] Dynamische Radlastverlagerung rechtes Hinterrad

     % Radlasten
     FWZ_vl(i+1) = FWZ_vl(1) + dFWZvl_aero(i) + dFWZvl_x(i) + dFWZvl_y(i); % [N] Radlast vorne links
     FWZ_vr(i+1) = FWZ_vr(1) + dFWZvr_aero(i) + dFWZvr_x(i) + dFWZvr_y(i); % [N] Radlast vorne rechts
     FWZ_hl(i+1) = FWZ_hl(1) + dFWZhl_aero(i) + dFWZhl_x(i) + dFWZhl_y(i); % [N] Radlast hinten links
     FWZ_hr(i+1) = FWZ_hr(1) + dFWZhr_aero(i) + dFWZhr_x(i) + dFWZhr_y(i); % [N] Radlast hinten rechts
     
     % Begrenzen der Radlasten auf (quasi) Null
     if FWZ_vl(i+1) < 0
         FWZ_vl(i+1) = 0.001;
     end
     if FWZ_vr(i+1) < 0
         FWZ_vr(i+1) = 0.001;
     end
     if FWZ_hl(i+1) < 0
         FWZ_hl(i+1) = 0.001;
     end
     if FWZ_hr(i+1) < 0
         FWZ_hr(i+1) = 0.001;
     end
     
     % Achslasten
     FWZh(i+1) = FWZ_hl(i+1) + FWZ_hr(i+1);  % [N] Achslast hinten
     FWZv(i+1) = FWZ_vl(i+1) + FWZ_vr(i+1);  % [N] Achslast vorne
     FWZges(i+1) = FWZh(i+1) + FWZv(i+1);    % [N] Gesamtachs- und Radlast    
           
    % Vertikale Reifensteifigkeiten
    cZ_vl(i+1) = interp1(Fz,cZ_tire,FWZ_vl(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    cZ_vr(i+1) = interp1(Fz,cZ_tire,FWZ_vr(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    cZ_hl(i+1) = interp1(Fz,cZ_tire,FWZ_hl(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert
    cZ_hr(i+1) = interp1(Fz,cZ_tire,FWZ_hr(i+1),'linear','extrap'); % [N/m] Reifensteifigkeit interpoliert

    % Dynamische Reifenradien
    Rdyn_vl(i+1) = R0 - FWZ_vl(i+1)/cZ_vl(i+1);    % [m] Dynamischer Reifenradius 
    Rdyn_vr(i+1) = R0 - FWZ_vr(i+1)/cZ_vr(i+1);    % [m] Dynamischer Reifenradius 
    Rdyn_hl(i+1) = R0 - FWZ_hl(i+1)/cZ_hl(i+1);    % [m] Dynamischer Reifenradius 
    Rdyn_hr(i+1) = R0 - FWZ_hr(i+1)/cZ_hr(i+1);    % [m] Dynamischer Reifenradius
   
    % Lenkwinkel, Schwimmwinkel, Gierrate, Schräglaufwinkel
    if R(i) > 0
        f = -1;
    else
        f = 1;
    end
    
    delta(i) = f*atan((l/1000)/sqrt(R(i)^2-(lh/1000)^2));   % [rad] Lenkwinkel
    beta(i) = f*atan((lh/1000)/sqrt(R(i)^2-(lh/1000)^2));   % [rad] Schwimmwinkel
    psi1(i) = vV(i)/R(i);                                   % [rad/s] Gierrate
    alpha_v(i) = 180/pi*(delta(i)-(lv/1000)/vV(i)*psi1(i)-beta(i));  % [°] Schräglaufwinkel vorne
    alpha_h(i) = 180/pi*((lh/1000)/vV(i)*psi1(i)-beta(i));           % [°] Schräglaufwinkel hinten

    % Maximal übertragbare Reifenkräfte in Längsrichtung
    FWXmax_vl(i+1) = max(abs(MF52_Fx_cs(alpha_v(i),FWZ_vl(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_vr(i+1) = max(abs(MF52_Fx_cs(alpha_v(i),FWZ_vr(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_hl(i+1) = max(abs(MF52_Fx_cs(alpha_h(i),FWZ_hl(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_hr(i+1) = max(abs(MF52_Fx_cs(alpha_h(i),FWZ_hr(i+1),GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximal übertragbare Radkraft
    FWXmax_v(i+1) = FWXmax_vl(i+1) + FWXmax_vr(i+1);    % [N] Maximal übertragbare Achskraft
    FWXmax_h(i+1) = FWXmax_hl(i+1) + FWXmax_hr(i+1);    % [N] Maximal übertragbare Achskraft
    
    % Maximal übertragbare Reifenkräfte in Querrichtung
    FWYmax_vl(i+1) = abs(MF52_Fy_cs(alpha_v(i),FWZ_vl(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_vr(i+1) = abs(MF52_Fy_cs(alpha_v(i),FWZ_vr(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_hl(i+1) = abs(MF52_Fy_cs(alpha_h(i),FWZ_hl(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_hr(i+1) = abs(MF52_Fy_cs(alpha_h(i),FWZ_hr(i+1),GAMMA,0,TIRparam)); % [N] Maximal übertragbare Radkraft
    FWYmax_v(i+1) = FWYmax_vl(i+1) + FWYmax_vr(i+1);    % [N] Maximal übertragbare Achskraft
    FWYmax_h(i+1) = FWYmax_hl(i+1) + FWYmax_hr(i+1);    % [N] Maximal übertragbare Achskraft
    
    t_z(Gang) = t_z(Gang) + t(i+1) - t(i); % [s] Zeit pro Gang
end

%% Plotten der Ergebnisse
% Plotten der Strecke
% for z = 1:32
%    [Val,Index100(z)] = min(abs(z*25-s)); 
% end
% figure(1)
% plot(x_Track,y_Track,'k')
% hold on
% for z = 1:length(ApexIndizes)
% scatter(x_Track(ApexIndizes),y_Track(ApexIndizes),20,'r','filled')
% text(x_Track(ApexIndizes(z)),y_Track(ApexIndizes(z)),[' ' num2str(z)])
% end
% title('AutoX-Strecke Formula Student Germany','FontSize',12) 
% text(x_Track(1),y_Track(1),'Start')
% text(x_Track(end),y_Track(end),'Ziel')
% % for z = 1:32
% %    scatter(x_Track(Index100(z)),y_Track(Index100(z)),20,'b','filled')
% %    hold on
% %    text(x_Track(Index100(z)),y_Track(Index100(z)),['  ' num2str(z*25) 'm'])
% % end
% annotation('textbox',[0.6, 0.8, 0.25, 0.05],'String',...
%     ['Track Length: ' num2str(max(s)) 'm'],'FitBoxToText','on')
% hold off
% box on

% Geschwindigkeit über Strecke
% figure(2)
% hold on
% scatter(s(Gangauswahl==1)',vV(Gangauswahl==1)*3.6,5,'b','filled')
% scatter(s(Gangauswahl==2)',vV(Gangauswahl==2)*3.6,5,'m','filled')
% scatter(s(Gangauswahl==3)',vV(Gangauswahl==3)*3.6,5,'g','filled')
% scatter(s(ApexIndizes),vAPEXmax*3.6,'r','filled')
% hold off
% title('Geschwindigkeitsverlauf','FontSize',12) 
% legend('1.Gang','2.Gang','3.Gang','Apexmaximalgeschwindigkeiten','Location','best')
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Geschwindigkeit [km/h]','FontSize',10)
% grid on
% grid minor
% box on
% 
% % Geschwindigkeitsverlauf auf Strecke
% figure(3)
% scatter(Track(:,1),Track(:,2),40,vV*3.6,'filled')
% title('Geschwindigkeitsverlauf','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% colormap(jet)
% c = colorbar;
% ylabel(c,'Geschwindigkeit [km/h]')
% box on
% axis equal
% 
% % Drehzahlverlauf über Strecke
% figure(4) 
% plot(s(1:end-1),ni,'b','LineWidth',1.5) 
% title('Motordrehzahl','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Drehzahl [min^-^1]','FontSize',10)
% grid on 
% box on
% % 
% % Drehzahlverlauf auf Strecke
% figure(5)
% scatter(Track(1:end-1,1),Track(1:end-1,2),40,ni,'filled')
% title('Drehzahlverlauf','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% colormap(jet)
% c = colorbar;
% ylabel(c,'Motordrehzahl [1/min]')
% box on
% axis equal
% 
% % Beschleunigungsverlauf
% figure(6)
% hold on
% plot(s(1:end-1),aVX./g,'b','LineWidth',0.5)
% plot(s(1:end-1),aVY./g,'r','LineWidth',0.5) 
% hold off
% title('Beschleunigungsverlauf','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Beschleunigung [g]','FontSize',10)
% legend('Längsbeschleunigung [g]','Querbeschleunigung [g]')
% grid on
% box on
% 
% % Längsbeschleunigungsverlauf auf Strecke
% figure(7)
% scatter(Track(1:end-1,1),Track(1:end-1,2),40,aVX/9.81,'filled')
% title('Längsbeschleunigungsverlauf','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% colormap(jet)
% c = colorbar;
% ylabel(c,'Längsbeschleunigung [g]')
% box on
% axis equal

% % Querbeschleunigungsverlauf auf Strecke
% figure(8)
% scatter(Track(1:end-1,1),Track(1:end-1,2),40,aVY/9.81,'filled')
% title('Querbeschleunigungsverlauf','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% colormap([flip(jet);jet])
% c = colorbar;
% ylabel(c,'Querbeschleunigung [g]')
% box on
% axis equal
% 
% % Gangauswahl 
% figure(9)
% plot(s(1:end-1),Gangauswahl,'b','LineWidth',1.5) 
% title('Gangauswahl','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Gang [-]','FontSize',10)
% grid on
% box on
% 
% Gangverlauf auf Strecke
figure(9)
scatter(Track(1:end-1,1),Track(1:end-1,2),40,Gangauswahl,'filled')
title('Gangauswahl','FontSize',12) 
xlabel('X-Koordinate [m]','FontSize',10)
ylabel('Y-Koordinate [m]','FontSize',10)
colormap(jet(Gang_max))
c = colorbar('Ticks',[1:1:Gang_max]);
caxis([1 Gang_max])
ylabel(c,'Gangauswahl')
box on
axis equal
% 
% % Zeit pro Gang
% figure(10)
% b = bar(1:Gang_max,t_z);
% set(b,'FaceColor','b')
% text(b.XEndPoints,b.YEndPoints,num2str(t_z'),'vert','bottom','horiz','center')
% title('Zeit pro Gang','FontSize',12)
% xlabel('Gang [-]','FontSize',10)
% ylabel('Zeit [s]','FontSize',10)
% 
% % Radlasten vorne
% figure(11)
% plot(s,FWZ_vl,s,FWZ_vr)
% title('Radlasten vorne')
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Radlast [N]','FontSize',10)
% legend('Radlast vorne links','Radlast vorne rechts')
% grid on
% box on
% 
% % Radlasten hinten
% figure(12)
% plot(s,FWZ_hl,s,FWZ_hr)
% title('Radlasten hinten')
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Radlast [N]','FontSize',10)
% legend('Radlast hinten links','Radlast hinten rechts')
% grid on
% grid minor
% box on
% 
% % Achslasten
% figure(13)
% plot(s,FWZh,s,FWZv)
% title('Achslasten')
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Achslast [N]','FontSize',10)
% legend('Achslast hinten','Achslast vorne')
% grid on
% box on

% % Dynamische Reifenradien
% figure(14)
% plot(s,Rdyn_vl,s,Rdyn_vr,s,Rdyn_hl,s,Rdyn_hr)
% title('Dynamische Reifenradien','FontSize',12) 
% legend('Vorne links','Vorne rechts','Hinten links','Hinten rechts')
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Dynamischer Reifenradius [m]','FontSize',10)
% grid on
% box on
% grid on

% % Bremssignal
% figure(15)
% scatter(Track(1:end-1,1),Track(1:end-1,2),40,BPPsignal,'filled')
% title('Bremssignalverlauf','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% legend('Bremse nicht betätigt','FontSize',12) 
% colormap([0 1 0;1 0 0])
% box on
% axis equal

% % Schräglaufwinkel
% figure(40)
% plot(s,alpha_vl,s,alpha_vr,s,alpha_hl,s,alpha_hr)
% title('Schräglaufwinkel')
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Schräglaufwinkel [°]','FontSize',10)
% legend('Vorne links','Vorne rechts','Hinten links','Hinten rechts')
% ylim([-20 20])
% grid on
% box on

% 1. Lauf ohne Bremsen vs. 2. Lauf mit Bremsen
figure(20)
plot(s,vV_NoBrake*3.6,'g',s,vV*3.6)
hold on
scatter(s(ApexIndizes),vAPEXmax*3.6,'r','filled')
hold off
title('Geschwindigkeitsverlauf mit/ohne Bremsen','FontSize',12) 
legend('Ohne Bremsen','Mit Bremsen')
xlabel('Strecke [m]','FontSize',10) 
ylabel('Geschwindigkeit [km/h]','FontSize',10)
grid on
box on

% % Maximal übertragbare und tatsächliche Querkräfte Vorderachse
% figure(23)
% plot(s,FWYmax_v,s(1:end-1),abs(FWYv))
% title('Maximal übertragbare und tatsächliche Querkraft Vorderachse','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Querkraft [N]','FontSize',10)
% legend('Max. übertragbare Querkraft','Tatsächliche Querkraft')
% grid on
% grid minor
% box on

% % Maximal übertragbare und tatsächliche Querkräfte Hinterachse
% figure(24)
% plot(s,FWYmax_h,s(1:end-1),abs(FWYh))
% title('Maximal übertragbare und tatsächliche Querkraft Hinterachse','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Querkraft [N]','FontSize',10)
% legend('Max. übertragbare Querkraft','Tatsächliche Querkraft')
% grid on
% grid minor
% box on

% figure(40)
% plot(s(1:end-1),FVX,s,FWXmax_hr,s,FWXmax_hl,s,speicherhr,s,speicherhl)
% legend('Zugkraft','Mit Bremsen HR','Mit Bremsen HL','Ohne Bremsen HR','Ohne Bremsen HL')
% 
% figure(42)
% % plot(s,FWZ_hr,s,FWZ_hl,s,FWZ_vr,s,FWZ_vl,...
% %     s,speicherFWZhr,s,speicherFWZhl,s,speicherFWZvr,s,speicherFWZvl)
% plot(s(BrakeIndizes),FWZ_hr(BrakeIndizes),s(BrakeIndizes),FWZ_hl(BrakeIndizes),...
%     s(BrakeIndizes),speicherFWZhr(BrakeIndizes),s(BrakeIndizes),speicherFWZhl(BrakeIndizes))
% legend('Mit Bremsen HR','Mit Bremsen HL','Ohne Bremsen HR','Ohne Bremsen HL')

% % Begrenzung der Geschwindigkeit auf Strecke
% figure(25)
% hold on
% scatter(Track(Tirelimit==0,1),Track(Tirelimit==0,2),40,'k','filled')
% scatter(Track(Tirelimit==1,1),Track(Tirelimit==1,2),40,'b','filled')
% scatter(Track(Tirelimit==2,1),Track(Tirelimit==2,2),40,'g','filled')
% scatter(Track(Tirelimit==3,1),Track(Tirelimit==3,2),40,'r','filled')
% hold off
% title('Limits des Fahrzeugs','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% legend('Begrenzung durch Leistung','Begrenzung durch Reifen (Beschleunigung)',...
%     'Begrenzung durch Reifen (Bremsung)','Begrenzung durch Reifen quer')
% box on
% axis equal

% % Längsschlupf
% figure(80)
% plot(s,kappa_hl,s,kappa_hr,s,kappa_vl,s,kappa_vr)
% yline(0)
% title('Längsschlupf der Hinterräder')
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Längsschlupf [-]','FontSize',10)
% legend('Hinten links','Hinten rechts','Vorne links','Vorne rechts')
% ylim([-0.25 0.25])
% grid on
% grid minor
% box on

% % Giergeschwindigkeit 
% figure(90)
% plot(s,psi1*180/pi)
% title('Giergeschwindigkeit')
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Giergeschwindigkeit [°/s]','FontSize',10)
% grid on
% grid minor
% box on

% figure(80)
% plot(s,RutschenY_v,s,RutschenY_h)
% hold on
% scatter(s(ApexIndizes),ones(1,length(ApexIndizes)),20,'r','filled')

% % Verlauf der übertragbaren Querkräfte (Conti-Plot)
% figure(100)
% for f = 2000:-200:0
% plot([-12:0.01:12],MF52_Fy_ps([-12:0.01:12],f,GAMMA,0,TIRparam))
% hold on
% end
% legend('2000 N','1800 N','1600 N','1400 N','1200 N','1000 N','800 N','600 N','400 N','200 N','0 N','FontSize',15)
% grid on
% box on

% % Verlauf der Zugkräfte an beiden Hinterrädern 
% figure(70)
% plot(s(1:end-1),FVX_hl,s(1:end-1),FVX_hr)
% grid on
% box on

% % Gierrate
% figure(72)
% hold on
% scatter(s(R>0),psi1neu(R>0)*180/pi,5,'g')
% scatter(s(R(1:end-1)<0),psi1neu(R(1:end-1)<0)*180/pi,5,'r')
% hold off
% title('Gierrate')
% grid on 
% grid minor

% % Lenkwinkel Vorderachse
% figure(73)
% plot(s(1:end-1),delta*180/pi)
% title('Lenkwinkel Vorderachse')
% grid on 
% grid minor

% % Schwimmwinkel
% figure(74)
% plot(s(1:end-1),beta*180/pi)
% title('Schwimmwinkel')
% grid on 
% grid minor

% Schräglaufwinkel
figure(75)
plot(s(1:end-1),alpha_v,s(1:end-1),alpha_h)
legend('Vorne','Hinten')
title('Schräglaufwinkel')
grid on 
grid minor

% % Traktionsgrenze 
% figure(80)
% % plot(s(1:end-1),FVX,s,FWXmax_h)
% % legend('Zugkraft','Traktionsgrenze')
% scatter(Track(:,1),Track(:,2),40,TC,'filled')
% title('Einsatz Traktionskontrolle','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% legend('Traktionskontrolle betätigt','FontSize',12) 
% colormap([0 1 0;1 0 0])
% box on
% axis equal

%%  Ausgabe der Werte
disp(['AutoX-Zeit: ' num2str(t(end)) ' s']);
toc