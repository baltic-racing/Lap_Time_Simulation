clear all
close all
clc

RBG = imread("FSG 2019 AutoX.png"); % Einlesen des Fotos
BW = rgb2gray(RBG); % Konvertieren des Fotos in Graustufen
C = imcontour(BW,3);    % Erfassen der Konturen des Fotos
D = C(:,1502:end);      % Ermitteln der tatsächlich benötigten Punkte
D = D(:,6:3507);

close all

Scalefactor = 0.32; % Skalieren der Strecke auf Originallänge

X1 = D(2,:)*Scalefactor;  % Originalpunkte von Matlab
Y1 = D(1,:)*Scalefactor;  % Originalpunkte von Matlab

% 1. Spline durch Punkte legen
AnzahlXY_Punkte = 160;
t = [0 cumsum(sqrt(diff(X1).^2+diff(Y1).^2))];
t_int = linspace(0,t(end),AnzahlXY_Punkte);
X1_fit = spline(t,X1,t_int);    % X-Punkte des ersten Splines
Y1_fit = spline(t,Y1,t_int);    % Y-Punkte des ersten Splines

% 2. Geraden einbauen, gleiche Anzahl Punkte
IndexL = xlsread("Indizes_Straights_160.xlsx","B2:B15");
IndexH = xlsread("Indizes_Straights_160.xlsx","C2:C15");

X2_fit = X1_fit;
Y2_fit = Y1_fit;

for k = 1:length(IndexL)
    X2_fit(IndexL(k):1:IndexH(k)) = linspace(X2_fit(IndexL(k)),X2_fit(IndexH(k)),IndexH(k)-IndexL(k)+1);
    Y2_fit(IndexL(k):1:IndexH(k)) = linspace(Y2_fit(IndexL(k)),Y2_fit(IndexH(k)),IndexH(k)-IndexL(k)+1);
end

%%
close all
Y2_fit(15) = 15.6;
Y2_fit(17) = 14.9;
Y2_fit(23) = 5.4;
Y2_fit(24) = 7.2;
Y2_fit(25) = 12;
Y2_fit(26) = 14;
Y2_fit(27) = 13.5;
Y2_fit(42) = 47.5;
Y2_fit(48) = 72.5;
Y2_fit(59) = 40.2;
Y2_fit(64) = 58.8;
Y2_fit(65) = 61.5;
Y2_fit(80) = 75;
Y2_fit(81) = 73.1;
Y2_fit(82) = 70;
Y2_fit(90) = 37;
Y2_fit(133) = 87;
Y2_fit(134) = 90;
Y2_fit(136) = 94.5;
Y2_fit(148) = 85.4;

X2_fit(15) = 110;
X2_fit(17) = 102;
X2_fit(23) = 72;
X2_fit(24) = 67;
X2_fit(42) = 137.5;
X2_fit(43) = 136.5;
X2_fit(79) = 66;
X2_fit(86) = 49.3;
X2_fit(87) = 51.5;
X2_fit(90) = 47.1;
X2_fit(93) = 37;
X2_fit(101) = 21;
X2_fit(105) = 20.8;
X2_fit(106) = 22.5;
X2_fit(107) = 25;
X2_fit(124) = 5.6;
X2_fit(126) = 7.5;
X2_fit(133) = 43;
X2_fit(134) = 45.5;
X2_fit(160) = 141.7;

% 3. Finaler Spline durch Punkte
AnzahlXY_Punkte = 20000;
t = [0 cumsum(sqrt(diff(X2_fit).^2+diff(Y2_fit).^2))];
t_int = linspace(0,t(end),AnzahlXY_Punkte);
X = spline(t,X2_fit,t_int);
Y = spline(t,Y2_fit,t_int);

R = CurvatureRadius(X,Y);
ApexIndizes = Apexes(R);

% Krümmungsrichtung der Radien ermitteln
for k =2:length(X)-1
    
    m(k) = (Y(k)-Y(k-1))/(X(k)-X(k-1));
    
    if X(k) < X(k-1)
        Yfiktiv = m(k)*abs((X(k+1)-X(k)))+Y(k+1);
    else
       Yfiktiv = -m(k)*abs((X(k+1)-X(k)))+Y(k+1);
    end
    
    if X(k) < X(k-1) && Yfiktiv < Y(k)
        Vorzeichen(k) = 1;
        
    elseif X(k) < X(k-1) && Yfiktiv > Y(k)
        Vorzeichen(k) = -1;
        
    elseif X(k) > X(k-1) && Yfiktiv < Y(k)
        Vorzeichen(k) = -1;
        
    else X(k) > X(k-1) && Yfiktiv > Y(k)
        Vorzeichen(k) = 1;
    end
  
   
end

Vorzeichen(1) = Vorzeichen(2);
Vorzeichen = [Vorzeichen Vorzeichen(end)];

% Plotten
figure(1)
hold on
% scatter(X1_fit,Y1_fit,10,'r','filled')
% scatter(X,Y,5,'g','filled')
scatter(X(3050),Y(3050),100,'b','filled')
box on
axis equal
scatter(X,Y,40,Vorzeichen,'filled')
colormap([0 1 0;1 0 0])

% scatter(X(ApexIndizes),Y(ApexIndizes),60,'r')
% scatter(X2_fit,Y2_fit,30,'b','filled')
% for k = 1:1:length(X2_fit)
%     text(X2_fit(k),Y2_fit(k),[num2str(k)],'FontSize',7) 
% end
grid on
grid minor
axis equal

%% TRACK ERSTELLEN
Track = [X' Y']; % [x-y] Koordinaten der gesamten Strecke  
Track = unique(Track,'rows','stable');  % Entfernen von doppelten Punkten
Tracklength = [0; cumsum(sqrt(diff(Track(:,1)).^2+diff((Track(:,2))).^2))];
Track = [Track Tracklength];

% Manuelles Auffüllen der Kurvenradien bis zum letzen Element
R = [R R(end)];
% R = round(R,2);
R = R.*Vorzeichen;  % Hinzufügen des richtigen Vorzeichens vor die Kurvenradien
Track = [Track R'];

save('AutoXTrack.mat','Track')