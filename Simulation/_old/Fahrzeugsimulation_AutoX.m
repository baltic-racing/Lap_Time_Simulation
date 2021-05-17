%% Fahrzeugsimulation
clear all; clc; close all;

%% Matfiles einlesen
load('Dynojet Momentenverlauf.mat'); % Momentenverlauf von dynamischem Rollenprüfstand
load('Aero Downforce Daten.mat');    % Aero Daten
load('AutoXTrack.mat');              % Streckendaten
load('Rdyn_65kPA_IA0.mat');          % Lookup-Table für dynamischen Reifendurchmesser
 
x_Track = Track(:,1);       % [m] X-Koordinate der Strecke
y_Track = Track(:,2);       % [m] Y-Koordinate der Strecke
s = Track(:,3);             % [m] Verlauf der Streckenlänge
R = Track(:,4);             % [m] Kurvenradien
ApexIndizes = Apexes(R);    % [-] Indizes der Apexes

%% Fahrzeugdaten
m_ges = 270;    % [kg] Fahrzeuggesamtmasse inkl. Fahrer
g = 9.81;       % [m/s²] Erdbeschleunigung
FG = m_ges*g;   % [N] Gewichtskraft des Fahrzeugs

l = 1525;       % [mm] Radstand
B = 1500;       % [mm] Spurweite
h_COG = 280;    % [mm] Höhe Fahrzeugschwerpunkt
x_COG = 1400;   % [mm] x-Koordinate des Fahrzeugschwerpunkts in CAD
x_vA = 470;     % [mm] x-Koordinate der Vorderachse in CAD
x_COP = 1600;   % [mm] x-Koordinate Center of Pressure in CAD
aero_ph = (l-(l-(x_COP-x_vA)))/l;   % [-] Anteil Aerokraft auf Hinterachse
aero_pv = 1-aero_ph;                % [-] Anteil Aerokraft auf Vorderachse

m_ph = 52;      % [%] Prozentualer Radlastanteil Hinterachse
lv = l*m_ph/100;        % [mm] Abstand Vorderachse zu Fahrzeugschwerpunkt
lh = l*(100-m_ph)/100;  % [mm] Abstand Hinterachse zu Fahrzeugschwerpunkt

FB = 5000;      % [N] Maximale Bremskraft
my = 1.4;       % [-] Reibkoeffizient (hier noch konstant)

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

% i = [2.579 1.720 1.429 1.238]; % [-] Übersetzungen der einzelnen Gänge
i = [1.857 1.565 1.35 1.238 1.136]; % [-] Übersetzungen der einzelnen Gänge
ig = i.*i_Primaer*i_Sekundaer;      % [-] Gesamtübersetzungen
Gang_max = length(i);               % [-] Höchster Gang

%% Initialisieren von benötigten Variablen
BPPsignal = zeros(1,length(Track)-1);
FB = FB*ones(1,length(Track));
Gangauswahl = zeros(1,length(Track)-1);
ni = zeros(1,length(Track)-1);
Faero = zeros(1,length(Track)-1);
Mi = zeros(1,length(Track)-1);
FVX = zeros(1,length(Track)-1);
FVXre = zeros(1,length(Track)-1);
FVXid = zeros(1,length(Track)-1);
FR = zeros(1,length(Track)-1);
FL = zeros(1,length(Track)-1);
Fdr = zeros(1,length(Track)-1);
FN = zeros(1,length(Track)-1);
FVY = zeros(1,length(Track)-1);
vVYmax = zeros(1,length(Track)-1);
aVX = zeros(1,length(Track)-1);
aVY = zeros(1,length(Track)-1);
vV = zeros(1,length(Track)-1);
t = zeros(1,length(Track)-1);
dFZhl_aero = zeros(1,length(Track)-1);
dFZhr_aero = zeros(1,length(Track)-1);
dFZvl_aero = zeros(1,length(Track)-1);
dFZvr_aero = zeros(1,length(Track)-1);
dFZhl_x = zeros(1,length(Track)-1);
dFZhr_x = zeros(1,length(Track)-1);
dFZvl_x = zeros(1,length(Track)-1);
dFZvr_x = zeros(1,length(Track)-1);
dFZhl_y = zeros(1,length(Track)-1);
dFZhr_y = zeros(1,length(Track)-1);
dFZvl_y = zeros(1,length(Track)-1);
dFZvr_y = zeros(1,length(Track)-1);
FZhl = zeros(1,length(Track)-1);
FZhr = zeros(1,length(Track)-1);
FZvl = zeros(1,length(Track)-1);
FZvr = zeros(1,length(Track)-1);
FZh = zeros(1,length(Track)-1);
FZv = zeros(1,length(Track)-1);
r_dyn = zeros(1,length(Track)-1);
alpha = zeros(1,length(Track)-1);
aReverse = zeros(1,length(Track)-1);
vReverse = zeros(1,length(Track)-1);

%% Reifenmodell - Magic Tire Formula Model 5.2

TireInertia = 0.2;    % [kg*m²] Massenträgheitsmoment eines Reifens

% Matrizen ALPHA, GAMMA, KAPPA müssen gleich groß sein wenn nicht =0!
ALPHA = 0; % [°] Schräglaufwinkel
GAMMA = 0; % [°] Sturz
KAPPA = 0:0.01:0.2; % [-] Längsschlupf

% Variable für Pfad des Tir-Files
FileNameLocation = ('C:\Users\Sven Weishaupt\Dropbox\MATLAB\_MATLAB_\C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

% Tir-File in Struktur laden
TIRparam = loadTIR(FileNameLocation);

% Vertikal
R0 = TIRparam.UNLOADED_RADIUS; % [m] Fertigungsradius des Reifens
TIRparam.LMUX = 0.65;          % [-] Skalierungsfaktor der Reifen (0.75 für optimale Reifentemperatur, 0.6 für niedrige Reifentemperatur)

cZ_tire = TIRparam.VERTICAL_STIFFNESS;  % [N/m] Vertikale Steifigkeit des Reifens
cY_tire = 450;  % [N/°] Schräglaufsteifigkeit des Reifens
%%%%% Conti-Dokument nutzen!

% (Statische) Radlasten
FZges = m_ges*g;          % [N] Statische Gesamtachslast 
FZh(1) = m_ph/100*FZges;  % [N] Statische Achslast hinten
FZv(1) = FZges-FZh(1);    % [N] Statische Achslast vorne
FZvr(1) = FZv(1)/2;       % [N] Statische Radlast vorne rechts  
FZvl(1) = FZv(1)/2;       % [N] Statische Radlast vorne links
FZhr(1) = FZh(1)/2;       % [N] Statische Radlast hinten rechts  
FZhl(1) = FZh(1)/2;       % [N] Statische Radlast hinten links  

R_dyn_hl(1) = 0.23;
r_dyn(1) = R0 - FZh(1) /cZ_tire/2;               % [m] Dynamischer Reifenradius im Stand(1) -> statisch
Fx_Achse = MF52_Fx_cs(ALPHA,FZh(1),GAMMA,KAPPA,TIRparam)*2; % [N] Übertragbare Achskraft, insbesondere abhängig vom Schlupf
F_Slip = Fx_Achse;                             % [N] Zwischenspeicherung zum Plotten des Reifenschlupfes

%% Berechnen der maximalen Kurvengeschwindigkeiten für alle Apexes

for i = 1:length(ApexIndizes)
    k = 0;
    v1 = 5;
    v2 = 100;

    while v1 < v2 
        v1 = v1 + k; 
        Faero = interp1(v,FA,v1*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert    
        FZ(i) = FZh(1) + Faero/2;   
        v2 = sqrt(max(MF52_Fy_cs(0:-0.1:-12,FZ(i),0,0,TIRparam))*R(ApexIndizes(i))/m_ges);  
        k = k+0.001; 
    end
    
    VYmax(i) = mean([v1 v2]);   % [m/s] Maximalgeschwindigkeit für jede Apex
    
end

%% Startwerte für Simulation
Gang = 1;      % [-] Fahrzeug startet im ersten Gang
vV(1) = 0;     % [m/s] Geschwindigkeit
t(1) = 0;      % [s] Zeit
Sector = 1;    % [-] Anfangssektor des Fahrzeugs

t_x = 0;    % [s] Hilfsgröße zur Berücksichtigung der Zugkraftunterbechung während des Schaltens
t_z = zeros(1,Gang_max);    % [s] Hilfsgröße zur Bestimmung der Zeit pro Gang   

%% Simulation
for i = 1:length(Track)-1
    
    % Bestimmen von Motordrehzahl und Gang
    Gangauswahl(i) = Gang; % [-] Speichern des Gangverlaufs
    ni(i) = vV(i)*30/pi*ig(Gang)/r_dyn(i); % [1/min] Aktuelle Motordrehzahl ermitteln
    
    if ni(i) < 5000 && Gang > 1 % Runterschalten
       Gang = Gang-1;    
       ni(i) = vV(i)*30/pi*ig(Gang)/r_dyn(i); % [1/min] Drehzahl nach Gangwechsel ermitteln
    end
    
    if ni(i) >= n_shift % Gang erhöhen, wenn Schaltdrehzahl erreicht
        Gang = Gang + 1;
        t_x = t(i)-t(i-1); % Schrittweite der Hilfsgröße festlegen
            if Gang >= Gang_max  % Begrenzen der höchsten Gangzahl
                Gang = Gang_max;
            end 
        ni(i)= vV(i)*30/pi*ig(Gang)/r_dyn(i); % [1/min] Drehzahl nach Gangwechsel ermitteln      
            if ni(i)>= n_Mmax % Drehzahlbegrenzer
                ni(i) = n_Mmax;
            end      
    end
    
    % Bestimmen von Aero-Kräften und Motormoment
    Faero(i) = interp1(v,FA,vV(i)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
    Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motormoment interpoliert
    
    % Berechnen der Zugkraft am Hinterrad
    if t_x == 0
        FVX(i)= Mi(i)*ig(Gang)/r_dyn(i); % [N] Zugkraft an der Hinterachse
    elseif t_x >=0 && t_shift > t_x
        FVX(i) = 0;
        t_x = t_x+t(i)-t(i-1);
    elseif t_x >= t_shift
        t_x = 0;                            
        FVX(i)= Mi(i)*ig(Gang)/r_dyn(i); % [N] Zugkraft an der Hinterachse
    end
    
    % Zugkrafthyperbeln
    FVXre(i) = P_Mmax/vV(i);       % [N] Berechnen der realen Zugkrafthyperbel
    FVXid(i) = FVXre(i)/eta_ges;  % [N] Berechnen der idealen Zugkrafthyperbel
    
    % Fahrwiderstände
    FR(i) = k_R*m_ges*g;              % [N] Rollwiderstand
    FL(i) = rho_L*vV(i)^2/2*c_w*A_S;  % [N] Luftwiderstand
    Fdr(i) = FR(i)+FL(i);             % [N] Gesamtwiderstand
    
    %%%%%%%%
    % Radlast
    FN(i) = FG+Faero(i); % [N] Radlast vereinfachend berechnet
     
    % Fahrzeug
    FVY(i) = m_ges*vV(i)^2/R(i);    % [N] Zentripetalkraft
    vVYmax(i) = ((my*FN(i))^2/((m_ges/R(i))^2+(rho_L/2*A_S*c_w)^2))^(1/4);  % [m/s] Max. Kurvengeschwindigkeit
    
    aVX(i) = (FVX(i)-Fdr(i))/(m_ges);   % [m/s²] (Longitudinal) Beschleunigung
    aVY(i) = vV(i)^2/R(i);              % [m/s²] Querbeschleunigung
    vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Gesamt-Fahrzeuggeschwindigkeit
    t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);    % [s] Zeit

    if vV(i+1) > vVYmax(i)   % Begrenzen auf maximale Kurvengeschwindigkeit
        vV(i+1) = vVYmax(i);
    end
        
%      % Radlastverlagerung in Folge von Aerokräften 
%      dFZhl_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf linkes Hinterrad
%      dFZhr_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf rechtes Hinterrad
%      dFZvl_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf linkes Vorderrad
%      dFZvr_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf rechtes Vorderrad
%      
%      % Dynamische Radlastverlagerungen in Längsrichtung
%      dFZhl_x(i) = m_ges*aXV(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung linkes Hinterrad
%      dFZhr_x(i) = m_ges*aXV(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung rechtes Hinterrad
%      dFZvl_x(i) = -m_ges*aXV(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung linkes Vorderrad
%      dFZvr_x(i) = -m_ges*aXV(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung rechtes Vorderrad
%                  
%      % Dynamische Radlastverlagerung in Querrichtung
%      dFZhl_y(i) = h_COG/B*lh/l*FVY(i);  % [N] Dynamische Radlastverlagerung linkes Hinterrad
%      dFZhr_y(i) = h_COG/B*lh/l*FVY(i);  % [N] Dynamische Radlastverlagerung rechtes Hinterrad
%      dFZvl_y(i) = h_COG/B*lv/l*FVY(i);  % [N] Dynamische Radlastverlagerung linkes Vorderrad
%      dFZvr_y(i) = h_COG/B*lv/l*FVY(i);  % [N] Dynamische Radlastverlagerung rechtes Vorderrad
%      
%      % Radlasten
%      FZhl(i+1) = FZhl(1) + dFZhl_aero(i) + dFZhl_x(i) + dFZhl_y(i); % [N] Radlast hinten links
%      FZhr(i+1) = FZhr(1) + dFZhr_aero(i) + dFZhr_x(i) + dFZhr_y(i); % [N] Radlast hinten rechts
%      FZvl(i+1) = FZvl(1) + dFZvl_aero(i) + dFZvl_x(i) + dFZvl_y(i); % [N] Radlast vorne links
%      FZvr(i+1) = FZvr(1) + dFZvr_aero(i) + dFZvr_x(i) + dFZvr_y(i); % [N] Radlast vorne rechts
%      
%      % Achslasten
%      FZh(i+1) = FZhl(i+1) + FZhr(i+1);  % [N] Achslast hinten
%      FZv(i+1) = FZvl(i+1) + FZvr(i+1);  % [N] Achslast vorne
%        
%     FX_Rad = MF52_Fx_cs(ALPHA,Fz_h(i+1),GAMMA,KAPPA,TIRparam); % [N] Übertragbare Radkraft
%     Fx_Achse(i+1) = max(FX_Rad)*2;        % [N] Übertragbare Achskraft  
%     
%     if FVX_x(k) >= Fx_Achse(k+1) % Begrenzung der Zugkraft auf die Traktionsgrenze
%         FVX_x(k) = Fx_Achse(k+1);
%     end
    
    %%%%%%%%
    r_dyn(i+1) = R0 - FN(1)/cZ_tire;    % [m] Veränderlicher dyn. Reifenradius  
    
    t_z(Gang) = t_z(Gang)+t(i+1)-t(i); % [s] Zeit pro Gang
end

vVYmax = [vVYmax vVYmax(end)]; % Erweitern um Wert des letzten Punkts

%% BREMSPUNKTBERECHNUNG - Berechnung erfolgt rückwärts
BrakeIndizes = [];
Bremsstrecke = [];

for p = length(Track)-1:-1:1
    Counter = 0;
    i_rev = p;
    vReverse(i_rev) = vVYmax(i_rev);
    while vReverse(i_rev) < vV(i_rev)
        FR(i_rev) = k_R*m_ges*g;              % [N] Rollwiderstand
        FL(i_rev) = rho_L*vReverse(i_rev)^2/2*c_w*A_S;  % [N] Luftwiderstand
        Fdr(i_rev) = FR(i_rev)+FL(i_rev);             % [N] Gesamtwiderstand
        aReverse(i_rev) = (-Fdr(i_rev)-FB(i_rev))/m_ges;   % NEGATIV
        vReverse(i_rev-1) = sqrt(vReverse(i_rev)^2-2*aReverse(i_rev)*(s(i_rev)-s(i_rev-1)));
        i_rev = i_rev-1;       
        Counter = Counter+1;
    end
    
    if Counter > 2
        BrakeIndizes = [BrakeIndizes i_rev];
        SectorExitIndizes(p) = p;
%         vV(i_rev:p) = vReverse(i_rev:p);
        Bremsstrecke = [Bremsstrecke i_rev:p];
    end
    
    
end

Bremsstrecke = unique(Bremsstrecke);  % Entfernen von doppelten Punkten

%% NEUE BREMSBERECHNUNG ÜBER APEXINDIZES
BrakeIndizes = [];
Bremsstrecke = [];

for k = 1:length(ApexIndizes)
    
    Counter = 0;
    j = ApexIndizes(k);
    vRev(j) = vVYmax(j);
    
    while vRev(j) < vV(j-1)
        FR(j) = k_R*m_ges*g;
        FL(j) = rho_L*vRev(j)^2/2*c_w*A_S;
        Fdr(j) = FR(j)+FL(j);
        aRev(j) = (-Fdr(j)-FB(j))/m_ges;
        vRev(j-1) = sqrt(vRev(j)^2-2*aRev(j)*(s(j)-s(j-1)));
        j = j - 1;
        Counter = Counter + 1;
    end
   
    if Counter > 1
        BrakeIndizes = [BrakeIndizes j];
        Bremsstrecke = [Bremsstrecke j:ApexIndizes(k)];
    end
    
end

%     j = ApexIndizes(4);
%     vRev(j) = 5;
%     
%     while vRev(j) < vV(j-1)
%         FR(j) = k_R*m_ges*g;
%         FL(j) = rho_L*vRev(j)^2/2*c_w*A_S;
%         Fdr(j) = FR(j)+FL(j);
%         aRev(j) = (-Fdr(j)-FB(j))/m_ges;
%         vRev(j-1) = sqrt(vRev(j)^2-2*aRev(j)*(s(j)-s(j-1)));
%         j = j - 1;
% 
%     end
% 
%         BrakeIndizes = [BrakeIndizes j];
%         Bremsstrecke = [Bremsstrecke j:ApexIndizes(4)];
    

% scatter(x_Track,y_Track,'b')
% hold on
% scatter(x_Track(Bremsstrecke),y_Track(Bremsstrecke),'r')
% scatter(x_Track(BrakeIndizes),y_Track(BrakeIndizes),40,'g')
% scatter(x_Track(ApexIndizes),y_Track(ApexIndizes),50,'g')

%% Startwerte für Simulation MIT BREMSEN
Gang = 1;      % [-] Fahrzeug startet im ersten Gang
vV(1) = 0;     % [m/s] Geschwindigkeit
t(1) = 0;      % [s] Zeit
omega = 0;
ff = 0;

t_x = 0;    % [s] Hilfsgröße zur Berücksichtigung der Zugkraftunterbechung während des Schaltens
t_z = zeros(1,Gang_max);    % [s] Hilfsgröße zur Bestimmung der Zeit pro Gang

%% Simulation MIT BREMSEN
for i = 1:length(Track)-1
    
    % Bestimmen von Motordrehzahl und Gang
    Gangauswahl(i) = Gang; % [-] Speichern des Gangverlaufs
    ni(i) = vV(i)*30/pi*ig(Gang)/R_dyn_hl(i); % [1/min] Aktuelle Motordrehzahl ermitteln
    
    if ni(i) < 5000 && Gang > 1 % Runterschalten
       Gang = Gang-1;    
       ni(i) = vV(i)*30/pi*ig(Gang)/R_dyn_hl(i); % [1/min] Drehzahl nach Gangwechsel ermitteln
    end
    
    if ni(i) >= n_shift % Gang erhöhen, wenn Schaltdrehzahl erreicht
        Gang = Gang + 1;
        t_x = t(i)-t(i-1); % Schrittweite der Hilfsgröße festlegen
            if Gang >= Gang_max  % Begrenzen der höchsten Gangzahl
                Gang = Gang_max;
            end 
        ni(i)= vV(i)*30/pi*ig(Gang)/R_dyn_hl(i); % [1/min] Drehzahl nach Gangwechsel ermitteln      
            if ni(i)>= n_Mmax % Drehzahlbegrenzer
                ni(i) = n_Mmax;
            end      
    end
    
    % Bestimmen von Aero-Kräften und Motormoment
    Faero(i) = interp1(v,FA,vV(i)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
    Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motormoment interpoliert
    
    % Berechnen der Zugkraft am Hinterrad
    if t_x == 0
        FVX(i)= Mi(i)*ig(Gang)/R_dyn_hl(i); % [N] Zugkraft an der Hinterachse
    elseif t_x >=0 && t_shift > t_x
        FVX(i) = 0;
        t_x = t_x+t(i)-t(i-1);
    elseif t_x >= t_shift
        t_x = 0;
        FVX(i)= Mi(i)*ig(Gang)/R_dyn_hl(i); % [N] Zugkraft an der Hinterachse
    end
    
    % Fahrwiderstände
    FR(i) = k_R*FZges(i);             % [N] Rollwiderstand
    FL(i) = rho_L*vV(i)^2/2*c_w*A_S;  % [N] Luftwiderstand
    Fdr(i) = FR(i)+FL(i);             % [N] Gesamtwiderstand
    
    %%%%%%%%
    % Radlast %%%%%%%%%%
    FN(i) = FG+Faero(i); % [N] Radlast vereinfachend berechnet
    
    % Prüfen, ob gebremst werden muss
    if ismember(i,Bremsstrecke) % Einleiten des Bremsvorgangs     
        BPPsignal(i) = 1;   %  Bremssignal
        FVX(i) = 0;         % [N] Zugkraft  
    else
        FB(i) = 0;          % [N] Bremskraft
    end
        
    % Gleichungen Fahrzeug
    FVY(i) = m_ges*vV(i)^2/R(i);    % [N] Zentripetalkraft
    vVYmax(i) = ((my*FN(i))^2/((m_ges/R(i))^2+(rho_L/2*A_S*c_w)^2))^(1/4);  % [m/s] Max. Kurvengeschwindigkeit
    
    aVX(i) = (FVX(i)-Fdr(i)-FB(i))/(m_ges);         % [m/s²] Longitudinal Beschleunigung
    aVY(i) = vV(i)^2/R(i);                          % [m/s²] Querbeschleunigung
    vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Gesamt-Fahrzeuggeschwindigkeit
    t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Zeit
 
%         % Fahrzeug
%         FVY(i) = 0;    % [N] Zentripetalkraft
%         vVYmax(i) = Inf;   % [m/s] Maximale Kurvengeschwindigkeit

    if vV(i+1) > vVYmax(i)   % Begrenzen auf maximale Kurvengeschwindigkeit
        vV(i+1) = vVYmax(i);
    end
    
     % Radlastverlagerung in Folge von Aerokräften 
     dFZhl_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf linkes Hinterrad
     dFZhr_aero(i) = Faero(i)/2*aero_ph;   % [N] Aerokraft auf rechtes Hinterrad
     dFZvl_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf linkes Vorderrad
     dFZvr_aero(i) = Faero(i)/2*aero_pv;   % [N] Aerokraft auf rechtes Vorderrad
     
     % Dynamische Radlastverlagerungen in Längsrichtung
     dFZhl_x(i) = m_ges*aVX(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung linkes Hinterrad
     dFZhr_x(i) = m_ges*aVX(i)*h_COG/l/2;   % [N] Dynamische Radlastverlagerung rechtes Hinterrad
     dFZvl_x(i) = -m_ges*aVX(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung linkes Vorderrad
     dFZvr_x(i) = -m_ges*aVX(i)*h_COG/l/2;  % [N] Dynamische Radlastverlagerung rechtes Vorderrad
     
     % Dynamische Radlastverlagerung in Querrichtung
     dFZhl_y(i) = h_COG/B*lh/l*FVY(i);   % [N] Dynamische Radlastverlagerung linkes Hinterrad
     dFZhr_y(i) = -h_COG/B*lh/l*FVY(i);  % [N] Dynamische Radlastverlagerung rechtes Hinterrad
     dFZvl_y(i) = h_COG/B*lv/l*FVY(i);   % [N] Dynamische Radlastverlagerung linkes Vorderrad
     dFZvr_y(i) = -h_COG/B*lv/l*FVY(i);  % [N] Dynamische Radlastverlagerung rechtes Vorderrad

     % Radlasten
     FZhl(i+1) = FZhl(1) + dFZhl_aero(i) + dFZhl_x(i) + dFZhl_y(i); % [N] Radlast hinten links
     FZhr(i+1) = FZhr(1) + dFZhr_aero(i) + dFZhr_x(i) + dFZhr_y(i); % [N] Radlast hinten rechts
     FZvl(i+1) = FZvl(1) + dFZvl_aero(i) + dFZvl_x(i) + dFZvl_y(i); % [N] Radlast vorne links
     FZvr(i+1) = FZvr(1) + dFZvr_aero(i) + dFZvr_x(i) + dFZvr_y(i); % [N] Radlast vorne rechts
     
     % Achslasten
     FZh(i+1) = FZhl(i+1) + FZhr(i+1);  % [N] Achslast hinten
     FZv(i+1) = FZvl(i+1) + FZvr(i+1);  % [N] Achslast vorne
     FZges(i+1) = FZh(i+1) + FZv(i+1);  % [N] Gesamtachs- und Radlast
     
     % Schräglaufwinkel aller Reifen muss hier berechnet werden
     alpha(i+1) = abs(FVY(i))/cY_tire;  % [°] Schräglaufwinkel
     FYmax(i) = MF52_Fy_cs(alpha(i),FZh(i),GAMMA,0,TIRparam);   % Test     
           
    %%%%%%%%
    r_dyn(i+1) = R0 - FN(1)/cZ_tire;    % [m] Veränderlicher dyn. Reifenradius  
    R_dyn_hl(i+1) = interp1(Fz,Rdyn,FZhl(i),'linear','extrap'); % [m] Dynamischer Reifenradius interpoliert
    
    FT(i) = FVX(i)/2-FR(i)/4-FB(i)/4;   % [N] Umfangskräfte an einem Hinterrad
    
    my0 = 2;    % [-] Haftreibungsbeiwert
    kappaCarter(i+1) = 0.05/(2*R_dyn_hl(i+1))*my0*(1-sqrt(1-FT(i)/(my0*FZhl(i+1))));    % [-] Längsschlupf nach Carter

    % Traktionsgrenzen aller Reifen in x- und y-Richtung
    FXhl(i+1) = MF52_Fx_cs(ALPHA,FZhl(i+1),GAMMA,kappaCarter(i+1),TIRparam);
    
    Fx_Rad = MF52_Fx_cs(ALPHA,FZh(i+1),GAMMA,KAPPA,TIRparam); % [N] Übertragbare Radkraft
    Fx_Achse(i+1) = max(Fx_Rad);        % [N] Übertragbare Achskraft  
    
%     if FVX_x(k) >= Fx_Achse(k+1) % Begrenzung der Zugkraft auf die Traktionsgrenze
%         FVX_x(k) = Fx_Achse(k+1);
%     end
    
    t_z(Gang) = t_z(Gang)+t(i+1)-t(i); % [s] Zeit pro Gang
end

%% Plotten der Ergebnisse
% Plotten der Strecke
for z = 1:16
   [Val,Index100(z)] = min(abs(z*50-s)); 
end
figure(1)
plot(x_Track,y_Track,'k')
hold on
title('AutoX-Strecke Formula Student Germany','FontSize',12) 
text(x_Track(1),y_Track(1),'Start')
text(x_Track(end),y_Track(end),'Ziel')
for z = 1:16
   scatter(x_Track(Index100(z)),y_Track(Index100(z)),20,'b','filled')
   hold on
   text(x_Track(Index100(z)),y_Track(Index100(z)),['  ' num2str(z*50) 'm'])
end
annotation('textbox',[0.6, 0.8, 0.25, 0.05],'String',...
    ['Track Length: ' num2str(max(s)) 'm'],'FitBoxToText','on')
hold off
box on

% Geschwindigkeit über Strecke
figure(2)
hold on
plot(s',vV*3.6)
scatter(s(ApexIndizes),VYmax*3.6,'r','filled')
title('Geschwindigkeitsverlauf','FontSize',12) 
xlabel('Strecke [m]','FontSize',10) 
ylabel('Geschwindigkeit [km/h]','FontSize',10)
grid on
box on

% Geschwindigkeitsverlauf auf Strecke
figure(3)
scatter(Track(:,1),Track(:,2),40,vV*3.6,'filled')
title('Geschwindigkeitsverlauf','FontSize',12) 
xlabel('X-Koordinate [m]','FontSize',10)
ylabel('Y-Koordinate [m]','FontSize',10)
colormap(jet)
c = colorbar;
ylabel(c,'Geschwindigkeit [km/h]')
box on
axis equal

% % Drehzahlverlauf über Strecke
% figure(4) 
% plot(s(1:end-1),ni,'b','LineWidth',1.5) 
% title('Motordrehzahl','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10) 
% ylabel('Drehzahl [min^-^1]','FontSize',10)
% grid on 
% box on
% 
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
% % Beschleunigungsverlauf auf Strecke
% figure(7)
% scatter(Track(1:end-1,1),Track(1:end-1,2),40,aVX/9.81,'filled')
% title('Beschleunigungsverlauf','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% colormap(jet)
% c = colorbar;
% ylabel(c,'Beschleunigung [g]')
% box on
% axis equal
% 
% % Gangauswahl 
% figure(8)
% plot(s(1:end-1),Gangauswahl,'b','LineWidth',1.5) 
% title('Gangauswahl','FontSize',12) 
% xlabel('Strecke [m]','FontSize',10)
% ylabel('Gang [-]','FontSize',10)
% grid on
% box on
% 
% % Gangverlauf auf Strecke
% figure(9)
% scatter(Track(1:end-1,1),Track(1:end-1,2),40,Gangauswahl,'filled')
% title('Gangauswahl','FontSize',12) 
% xlabel('X-Koordinate [m]','FontSize',10)
% ylabel('Y-Koordinate [m]','FontSize',10)
% colormap(jet(Gang_max))
% c = colorbar('Ticks',[1:1:Gang_max]);
% caxis([1 5])
% ylabel(c,'Gangauswahl')
% box on
% axis equal
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
% % % Radlasten
% % figure(11)
% % plot(s,FZhl,s,FZhr,s,FZvl,s,FZvr)
% % title('Radlasten')
% % xlabel('Strecke [m]','FontSize',10)
% % ylabel('Radlast [N]','FontSize',10)
% % legend('Radlast hinten links','Radlast hinten rechts',...
% %     'Radlast vorne links','Radlast vorne rechts')
% % grid on
% % box on
% 
% % % Achslasten
% % figure(12)
% % plot(s,FZh,s,FZv)
% % title('Achslasten')
% % xlabel('Strecke [m]','FontSize',10)
% % ylabel('Achslast [N]','FontSize',10)
% % legend('Achslast hinten','Achslast vorne')
% % grid on
% % box on
% 
% % % Schräglaufwinkel über Strecke %%%%%%%%%%%%
% % figure(13)
% % plot(s,alpha)
% % title('Schräglaufwinkel')
% % xlabel('Strecke [m]','FontSize',10)
% % ylabel('Schräglaufwinkel [°]','FontSize',10)
% % grid on
% % box on
% 
% % % Vergleich Reifenmodell und Zentripetalkraft
% % figure(14)
% % plot(s(1:end-1),FYmax,s(1:end-1),Fzp)
% % legend('MF Tire Model','Zentripetalkraft')
% 
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
% 
% % figure(16)
% % scatter(Track(1:end,1),Track(1:end,2),30,vVYmax,'filled')
% % hold on
% % scatter(Track(ApexIndizes,1),Track(ApexIndizes,2),50,'r','filled')
% 
% % figure(17)
% % plot(s,R_dyn_hl)
% % title('Dynamischer Reifenradius hinten links','FontSize',12) 
% % xlabel('Streckenlänge [m]','FontSize',10)
% % ylabel('Dynamischer Reifenradius [m]','FontSize',10)

% figure(18)
% plot(s(Bremsstrecke),kappaCarter(Bremsstrecke))
% 
% figure(19)
% plot(s,FXhl,'b',s(1:end-1),(FVX-Fdr),'r')

%%  Ausgabe der Werte
disp(['AutoX-Zeit: ' num2str(t(end)) ' s']);