clear all
close all
clc

addpath(genpath('D:\Users\Sven Weishaupt\Desktop\Vehicle Simulation'))

% Einlesen der Daten
load('Endurance Inner Coordinates.mat')
load('Endurance Outer Coordinates.mat')

% Anpassen der originalen Koordinaten
X_Coordinate(11) = [];
Y_Coordinate(11) = [];
Z_Coordinate(11) = [];
X_Coordinate = X_Coordinate(1:260);
Y_Coordinate = Y_Coordinate(1:260);
Z_Coordinate = Z_Coordinate(1:260);
X_Coordinate(end+1) = X_Coordinate(1);
Y_Coordinate(end+1) = Y_Coordinate(1);
Z_Coordinate(end+1) = Z_Coordinate(1);

X_outer(182) = [];
Y_outer(182) = [];
Z_outer(182) = [];
X_outer(181) = [];
Y_outer(181) = [];
Z_outer(181) = [];
X_outer(2) = [];
Y_outer(2) = [];
Z_outer(2) = [];

% Umsortieren der Punkte für Startlinie
X_Coordinate = [X_Coordinate(68:end); X_Coordinate(2:68)];
Y_Coordinate = [Y_Coordinate(68:end); Y_Coordinate(2:68)];
Z_Coordinate = [Z_Coordinate(68:end); Z_Coordinate(2:68)];

% Konvertieren der inneren Koordinaten
lla = [X_Coordinate Y_Coordinate Z_Coordinate];
p = lla2ecef(lla);
Xi = p(:,1);
Yi = p(:,2);
Zi = p(:,3);

% Skalieren der Rennstrecke
Scalefactor = 1.013; % Skalieren der Strecke auf Originallänge
Xi = Xi*Scalefactor;  % Originalpunkte von Matlab
Yi = Yi*Scalefactor;  % Originalpunkte von Matlab
Zi = Zi*Scalefactor;  % Originalpunkte von Matlab

% Anpassen der Rennstrecke zu umgänglichen Koordinaten
Xi = Xi-4.171*10^6;
Yi = Yi-6.28*10^5;
Zi = Zi-4.877*10^6;

% Verschieben einzelner Punkte + Einfügen des Slaloms
Yi(260) = Yi(260) + 0.1;
Yi(259) = Yi(259) + 0.3;
Yi(258) = Yi(258) + 0.4;
Yi(257) = Yi(257) + 0.3;

Xi(119) = 766;
Yi(119) = 387;
Xi(120) = 762.5;
Yi(120) = 390;
Yi(121) = 393;
Yi(122) = 395;
Xi(123) = 754;
Yi(123) = 395.5;
Xi(124) = 750;
Yi(124) = 395.5;
Xi(125) = 746;
Yi(125) = 398.5;
Xi(126) = 745;
Yi(126) = 401;
Xi(127) = 742;
Yi(127) = 402.5;

Xx(118:1:128) = linspace(Xi(118),Xi(128),11);
Yy(118:1:128) = linspace(Yi(118),Yi(128),11);

figure(1)
hold on
scatter(Xi,Yi,7,'filled')
scatter(Xx,Yy,7,'filled')
plot(Xi,Yi,'r')
for k = 1:2:length(Xi)
    text(Xi(k),Yi(k),num2str(k))
end
hold off
grid on
grid minor
box on
axis equal

%%
% TARGETLENGTH = 1.227e+03 m!!!!!!!!

% Value = [0;0];
% MaxRadius = 0;
% 
% for k = 1:500
% % 1. Spline durch Punkte legen
% AnzahlXY_Punkte = 120+k;
% t = [0; cumsum(sqrt(diff(Xi).^2+diff(Yi).^2+diff(Zi).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X1 = spline(t,Xi,t_int);    % X-Punkte des ersten Splines
% Y1 = spline(t,Yi,t_int);    % Y-Punkte des ersten Splines
% Z1 = spline(t,Zi,t_int);    % Y-Punkte des ersten Splines
% 
% % 2. Finaler Spline durch Punkte
% AnzahlXY_Punkte = 500;
% t = [0 cumsum(sqrt(diff(X1).^2+diff(Y1).^2+diff(Z1).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X2 = spline(t,X1,t_int);    % X-Punkte des ersten Splines
% Y2 = spline(t,Y1,t_int);    % Y-Punkte des ersten Splines
% Z2 = spline(t,Z1,t_int);    % Y-Punkte des ersten Splines
% 
% % 2. Finaler Spline durch Punkte
% AnzahlXY_Punkte = 1000;
% t = [0 cumsum(sqrt(diff(X2).^2+diff(Y2).^2+diff(Z2).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X3 = spline(t,X2,t_int);    % X-Punkte des ersten Splines
% Y3 = spline(t,Y2,t_int);    % Y-Punkte des ersten Splines
% Z3 = spline(t,Z2,t_int);    % Y-Punkte des ersten Splines
%  
% % 2. Finaler Spline durch Punkte
% AnzahlXY_Punkte = 4800;
% t = [0 cumsum(sqrt(diff(X3).^2+diff(Y3).^2+diff(Z3).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X4 = spline(t,X3,t_int);    % X-Punkte des ersten Splines
% Y4 = spline(t,Y3,t_int);    % Y-Punkte des ersten Splines
% Z4 = spline(t,Z3,t_int);    % Y-Punkte des ersten Splines
% 
% Track = [X4' Y4' Z4'];
% 
% % Radien und Apexes des Tracks ermitteln
% [L,R,K] = curvature(Track);
% 
% [Val,Index] = min(R);
% 
% if Val > 2.5
%     Value(:,end+1) = [Val;k];
% end
% 
% if Val > MaxRadius
%     MaxRadius = Val;
%     idx = k;
% end
% 
% end

%% Splines bei umgestellter Strecke
% 1. Spline durch Punkte legen
AnzahlXY_Punkte = 266;
t = [0; cumsum(sqrt(diff(Xi).^2+diff(Yi).^2+diff(Zi).^2))];
t_int = linspace(0,t(end),AnzahlXY_Punkte);
X1 = spline(t,Xi,t_int);    % X-Punkte des ersten Splines
Y1 = spline(t,Yi,t_int);    % Y-Punkte des ersten Splines
Z1 = spline(t,Zi,t_int);    % Y-Punkte des ersten Splines

% 2. Finaler Spline durch Punkte
AnzahlXY_Punkte = 500;
t = [0 cumsum(sqrt(diff(X1).^2+diff(Y1).^2+diff(Z1).^2))];
t_int = linspace(0,t(end),AnzahlXY_Punkte);
X2 = spline(t,X1,t_int);    % X-Punkte des ersten Splines
Y2 = spline(t,Y1,t_int);    % Y-Punkte des ersten Splines
Z2 = spline(t,Z1,t_int);    % Y-Punkte des ersten Splines

% 2. Spline durch Punkte legen
AnzahlXY_Punkte = 1000;
t = [0 cumsum(sqrt(diff(X2).^2+diff(Y2).^2+diff(Z2).^2))];
t_int = linspace(0,t(end),AnzahlXY_Punkte);
X3 = spline(t,X2,t_int);    % X-Punkte des ersten Splines
Y3 = spline(t,Y2,t_int);    % Y-Punkte des ersten Splines
Z3 = spline(t,Z2,t_int);    % Y-Punkte des ersten Splines

figure(4)
hold on
scatter(X3,Y3,7,'filled')
plot(X3,Y3,'r')
for k = 1:5:length(X3)
    text(X3(k),Y3(k),num2str(k))
end
hold off
grid on
box on
axis equal

%% Geraden einfügen

X3_fit = X3;
Y3_fit = Y3;

IndexL = [64 127 145 318 733 777 829 906 936];
IndexH = [90 138 176 344 748 814 850 922 956];

for k = 1:length(IndexL)
    X3_fit(IndexL(k):1:IndexH(k)) = linspace(X3_fit(IndexL(k)),X3_fit(IndexH(k)),IndexH(k)-IndexL(k)+1);
    Y3_fit(IndexL(k):1:IndexH(k)) = linspace(Y3_fit(IndexL(k)),Y3_fit(IndexH(k)),IndexH(k)-IndexL(k)+1);
end

figure(5)
% scatter(X3,Y3,5,'r','filled')
plot(X3,Y3,'r')
hold on
% scatter(X3_fit,Y3_fit,5,'b','filled')
plot(X3_fit,Y3_fit,'b')
hold off
grid on
box on
axis equal
 
%% 3. finalen Spline durch Punkte legen
AnzahlXY_Punkte = 4800;
t = [0 cumsum(sqrt(diff(X3_fit).^2+diff(Y3_fit).^2+diff(Z3).^2))];
t_int = linspace(0,t(end),AnzahlXY_Punkte);
X4 = spline(t,X3_fit,t_int);    % X-Punkte des ersten Splines
Y4 = spline(t,Y3_fit,t_int);    % Y-Punkte des ersten Splines
Z4 = spline(t,Z3,t_int);    % Y-Punkte des ersten Splines

Track = [X4' Y4' Z4'];
% Track = unique(Track,'rows','stable');  % Entfernen von doppelten Punkten

% Kurvenradien und Apexes bestimmen
[L,R,K] = curvature(Track);
min(R)
R(1) = R(2);
R(end) = R(1);
ApexIndizes = Apexes(R);

figure(3)
plot(1:length(R),R)
ylim([0 10])

%% Splines durch Original-Koordinaten legen
% % 1. Spline durch Punkte legen
% AnzahlXY_Punkte = 183;
% t = [0; cumsum(sqrt(diff(Xi).^2+diff(Yi).^2+diff(Zi).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X1 = spline(t,Xi,t_int);    % X-Punkte des ersten Splines
% Y1 = spline(t,Yi,t_int);    % Y-Punkte des ersten Splines
% Z1 = spline(t,Zi,t_int);    % Y-Punkte des ersten Splines
% 
% % 2. Finaler Spline durch Punkte
% AnzahlXY_Punkte = 500;
% t = [0 cumsum(sqrt(diff(X1).^2+diff(Y1).^2+diff(Z1).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X2 = spline(t,X1,t_int);    % X-Punkte des ersten Splines
% Y2 = spline(t,Y1,t_int);    % Y-Punkte des ersten Splines
% Z2 = spline(t,Z1,t_int);    % Y-Punkte des ersten Splines
% 
% % 2. Finaler Spline durch Punkte
% AnzahlXY_Punkte = 1000;
% t = [0 cumsum(sqrt(diff(X2).^2+diff(Y2).^2+diff(Z2).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X3 = spline(t,X2,t_int);    % X-Punkte des ersten Splines
% Y3 = spline(t,Y2,t_int);    % Y-Punkte des ersten Splines
% Z3 = spline(t,Z2,t_int);    % Y-Punkte des ersten Splines
%  
% % 2. Finaler Spline durch Punkte
% AnzahlXY_Punkte = 4800;
% t = [0 cumsum(sqrt(diff(X3).^2+diff(Y3).^2+diff(Z3).^2))];
% t_int = linspace(0,t(end),AnzahlXY_Punkte);
% X4 = spline(t,X3,t_int);    % X-Punkte des ersten Splines
% Y4 = spline(t,Y3,t_int);    % Y-Punkte des ersten Splines
% Z4 = spline(t,Z3,t_int);    % Y-Punkte des ersten Splines
% 
% Track = [X4' Y4' Z4'];
% % Track = unique(Track,'rows','stable');  % Entfernen von doppelten Punkten
% 
% % Kurvenradien und Apexes bestimmen
% [L,R,K] = curvature(Track);
% min(R)
% R(1) = R(2);
% R(end) = R(1);
% ApexIndizes = Apexes(R);
% 

%% Streckenlänge bestimmen
Tracklength = [0 cumsum(sqrt(diff(X4).^2+diff(Y4).^2+diff(Z4).^2))];
Tracklength(end)
Track = [Track Tracklength'];
Track = [Track R];

%% Plotten
figure(1)
plot(Xi,Yi)
hold on
plot(X1,Y1)
plot(X4,Y4)
scatter(X4(ApexIndizes),Y4(ApexIndizes),20,'g','filled')
scatter(X4(1),Y4(1),40,'b','filled')
legend('Original','Spline 1','Spline 4','Apexes','START')
grid on
box on

%% Steigungen bestimmen
% 
% % Vektor in Grundebene
% V0 = [X2(2)-X2(1);Y2(2)-Y2(1);Z2(2)-Z2(1)];
% 
% for k = 1:length(R)-1
%    
%     % Vektor erstellen
%     V = [X2(k+1)-X2(k);Y2(k+1)-Y2(k);Z2(k+1)-Z2(k)];
%     
%     Scalar(k) = dot(V,V0);
%     alpha(k) = acosd(Scalar(k)/(norm(V)*norm(V0)));
%     
% %     CosTheta = max(min(dot(V,V0)/(norm(V)*norm(V0)),1),-1);
% %     ThetaInDegrees = real(acosd(CosTheta));
%     
% end

% Track = [X' Y']; % [x-y] Koordinaten der gesamten Strecke  
% Track = unique(Track,'rows','stable');  % Entfernen von doppelten Punkten
% Tracklength = [0; cumsum(sqrt(diff(Track(:,1)).^2+diff((Track(:,2))).^2))];
% Tracklength(end)

%% Krümmungsrichtung der Radien ermitteln
% for k =2:length(X)-1
%     
%     m(k) = (Y(k)-Y(k-1))/(X(k)-X(k-1));
%     
%     if X(k) < X(k-1)
%         Yfiktiv = m(k)*abs((X(k+1)-X(k)))+Y(k+1);
%     else
%        Yfiktiv = -m(k)*abs((X(k+1)-X(k)))+Y(k+1);
%     end
%     
%     if X(k) < X(k-1) && Yfiktiv < Y(k)
%         Vorzeichen(k) = 1;
%         
%     elseif X(k) < X(k-1) && Yfiktiv > Y(k)
%         Vorzeichen(k) = -1;
        
%     elseif X(k) > X(k-1) && Yfiktiv < Y(k)
%         Vorzeichen(k) = -1;
%         
%     else X(k) > X(k-1) && Yfiktiv > Y(k)
%         Vorzeichen(k) = 1;
%     end
%   
%    
% end
% 
% Vorzeichen(1) = Vorzeichen(2);
% Vorzeichen = [Vorzeichen Vorzeichen(end)];

%% Speichern der .mat-Datei
save('EnduranceTrack.mat','Track')