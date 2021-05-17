clear all
close all
clc

RBG = imread("Hockenheimring.png"); % Einlesen des Fotos

BW = rgb2gray(RBG); % Konvertieren des Fotos in Graustufen

figure(1)
imshow(BW)  % Darstellen des Fotos in Graustufen

figure(2)
C = imcontour(BW,1);    % Erfassen der Konturen des Fotos
D = C(:,2:end);

[Max,Index] = max(D(2,:));

Kurve1_X = D(2,1:Index-1);
Kurve1_Y = D(1,1:Index-1);

Kurve2_X = D(2,Index+1:end);
Kurve2_Y = D(1,Index+1:end);

% Spline durch Punkte legen
t = [0 cumsum(sqrt(diff(Kurve1_X).^2+diff(Kurve1_Y).^2))];
t_int = linspace(0,t(end),500);
X1_fit = spline(t,Kurve1_X,t_int);
Y1_fit = spline(t,Kurve1_Y,t_int);

t = [0 cumsum(sqrt(diff(Kurve2_X).^2+diff(Kurve2_Y).^2))];
t_int = linspace(0,t(end),500);
X2_fit = spline(t,Kurve2_X,t_int);
Y2_fit = spline(t,Kurve2_Y,t_int);

% k=1;
% Idx = 0;
% Fehler = 1;

% for i=1:length(X1_fit)/100:length(X1_fit)
%     while abs(Fehler) > 0.5
%         Fehler = X1_fit(i) - X2_fit(k);
%         k=k+1;
%     end
%     Idx = [Idx k];
%     Fehler = 1;
% end

% Mx = mean([X1_fit;X2_fit]); 
% My = mean([Y1_fit;Y2_fit]);

% Plotten
plot(Kurve1_X,Kurve1_Y,'g',Kurve2_X,Kurve2_Y,'g')
hold on
plot(X1_fit,Y1_fit,'b',X2_fit,Y2_fit,'r')
% plot(Mx,My,'k')
title('Hockenheimring')
xlabel("X-Koordinate")
ylabel("Y-Koordinate")
axis equal

%% Methode mittels bwboundaries
% B = bwboundaries(BW);
% Curve1 = B{2};
% Curve2 = B{3};
% X1 = Curve1(:,1);
% Y1 = Curve1(:,2);
% X2 = Curve2(:,1);
% Y2 = Curve2(:,2);
% 
% figure(3)
% plot(X1,Y1,'b',X2,Y2,'r')
% title('Hockenheimring')
% xlabel("X-Koordinate")
% ylabel("Y-Koordinate")

%     % parameter - curve length
% t = [0; cumsum(sqrt(diff(X1).^2+diff(Y1).^2))];
% t_int = linspace(0,t(end),500);
% X1_fit = spline(t,X1,t_int);
% Y1_fit = spline(t,Y1,t_int);
% 
% t = [0; cumsum(sqrt(diff(X2).^2+diff(Y2).^2))];
% t_int = linspace(0,t(end),1000);
% X2_fit = spline(t,X2,t_int);
% Y2_fit = spline(t,Y2,t_int);
% 
% figure(4)
% plot(X1_fit,Y1_fit,'b',X1,Y1)
% legend('Fitted Curve','Data Points')
