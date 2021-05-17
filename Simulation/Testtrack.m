clear all
close all
clc

% Definieren der Kurven und Geraden

H = 50; % Anzahl Punkte je m für Geraden

Corner_1 = lefthander90(20,10,10);
Corner_2 = righthander180N(40,20,10);

Length_S2 = sqrt((Corner_2(1,1)-Corner_1(end,1))^2+(Corner_2(1,2)-Corner_1(end,2))^2);
Straight_2 = [linspace(Corner_1(end,1),Corner_2(1,1),Length_S2*H)'...
    linspace(Corner_1(end,2),Corner_2(1,2),Length_S2*H)'];

Corner_3 = righthander90S(30,-20,20);

Length_S3 = sqrt((Corner_3(1,1)-Corner_2(end,1))^2+(Corner_3(1,2)-Corner_2(end,2))^2);
Straight_3 = [linspace(Corner_2(end,1),Corner_3(1,1),Length_S3*H)'...
    linspace(Corner_2(end,2),Corner_3(1,2),Length_S3*H)'];

Corner_4 = righthander45W(-20,-30,10);

Length_S4 = sqrt((Corner_4(1,1)-Corner_3(end,1))^2+(Corner_4(1,2)-Corner_3(end,2))^2);
Straight_4 = [linspace(Corner_3(end,1),Corner_4(1,1),Length_S4*H)'...
    linspace(Corner_3(end,2),Corner_4(1,2),Length_S4*H)'];

Corner_5 = righthander135NW(-40,-10,10);

Length_S5 = sqrt((Corner_5(1,1)-Corner_4(end,1))^2+(Corner_5(1,2)-Corner_4(end,2))^2);
Straight_5 = [linspace(Corner_4(end,1),Corner_5(1,1),Length_S5*H)'...
    linspace(Corner_4(end,2),Corner_5(1,2),Length_S5*H)'];

Length_S1 = sqrt((Corner_1(1,1)-Corner_5(end,1))^2+(Corner_1(1,2)-Corner_5(end,2))^2);
Straight_1 = [linspace(Corner_5(end,1),Corner_1(1,1),Length_S1*H)'...
    linspace(Corner_5(end,2),Corner_1(1,2),Length_S1*H)'];

% Zusammenfügen der Strecke
Track = [Straight_1;Corner_1; % [x-y] Koordinaten der gesamten Strecke
    Straight_2;Corner_2;
    Straight_3;Corner_3;
    Straight_4;Corner_4;
    Straight_5;Corner_5];  
Track = unique(Track,'rows','stable');  % Entfernen von doppelten Punkten

% Längen der einzelnen Sektoren bestimmen
Sectorlength1 = [0; cumsum(sqrt(diff(Straight_1(:,1)).^2+diff((Straight_1(:,2))).^2))];
Sectorlength2 = [cumsum(sqrt(diff(Corner_1(:,1)).^2+diff((Corner_1(:,2))).^2))];
Sectorlength3 = [cumsum(sqrt(diff(Straight_2(:,1)).^2+diff((Straight_2(:,2))).^2))];
Sectorlength4 = [cumsum(sqrt(diff(Corner_2(:,1)).^2+diff((Corner_2(:,2))).^2))];
Sectorlength5 = [cumsum(sqrt(diff(Straight_3(:,1)).^2+diff((Straight_3(:,2))).^2))];
Sectorlength6 = [cumsum(sqrt(diff(Corner_3(:,1)).^2+diff((Corner_3(:,2))).^2))];
Sectorlength7 = [cumsum(sqrt(diff(Straight_4(:,1)).^2+diff((Straight_4(:,2))).^2))];
Sectorlength8 = [cumsum(sqrt(diff(Corner_4(:,1)).^2+diff((Corner_4(:,2))).^2))];
Sectorlength9 = [cumsum(sqrt(diff(Straight_5(:,1)).^2+diff((Straight_5(:,2))).^2))];
Sectorlength10 = [cumsum(sqrt(diff(Corner_5(:,1)).^2+diff((Corner_5(:,2))).^2))];

% Bestimmen der Länge der gesamten Strecke über einzelne Sektoren
Sectorlength = Sectorlength1;
Sectorlength = [Sectorlength; Sectorlength2+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength3+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength4+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength5+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength6+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength7+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength8+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength9+max(Sectorlength)];
Sectorlength = [Sectorlength; Sectorlength10+max(Sectorlength)];

Track = [Track Sectorlength(1:end-1)];

Tracklength = max(Sectorlength);    % [m] Gesamtlänge der Strecke

% %% Radien der gesamten Strecke berechnen
% for i = 1:length(Track)-2
%     
%     % Steigung zwischen zwei Punkten berechnen
%     m(i) = (Track(i+1,2)-Track(i,2))/(Track(i+1,1)-Track(i,1));
%     
%     % Winkel und Radien ausgeben lassen
%     a(i) = sqrt((Track(i+2,1)-Track(i,1))^2+(Track(i+2,2)-Track(i,2))^2);
%     b(i) = sqrt((Track(i+1,1)-Track(i,1))^2+(Track(i+1,2)-Track(i,2))^2);
%     c(i) = sqrt((Track(i+2,1)-Track(i+1,1))^2+(Track(i+2,2)-Track(i+1,2))^2);
%     A(i) = acos((b(i)^2+c(i)^2-a(i)^2)/(2*b(i)*c(i)));
%     A(i) = real(A(i));
%     R(i) = a(i)/(2*sin(pi-A(i)));
%         
% end
% 
% m(isinf(m)|isnan(m)) = 0;
% m = round(m,5);
% 
% for i = 1:length(Track)-2
%     if std([m(i) m(i+1)]) < 1e-13
%         R(i) = 0;
%     end
% end

% Manuelles Auffüllen der Kurvenradien bis zum letzen Element
% R = [R 10 10];
% R(2999) = 0;
% R(3499) = 10;
% R(3998) = 0;
% R(4998) = 10;
% R(6997) = 0;
% R(7497) = 20;
% R(9996) = 0;
% R(10246) = 10;
% R(11659) = 0;
% R = round(R,2);
% Track = [Track R'];

% SectorExitIndizes = [];
% for i = 1:length(Track)
%     if R(i) ~= R(i+1)
%         SectorExitIndizes = [SectorExitIndizes i];
%     end
% end
% 
% SectorExitIndizes = [SectorExitIndizes length(Track)];

% save('Testtrack_SectorExitIndizes.mat','SectorExitIndizes')
% save('Testtrack_Data.mat','Track')

%% TEST
for i = 2:length(Track)-1
    K(i) = 2*abs((Track(i,1)-Track(i-1,1)).*(Track(i+1,2)-Track(i-1,2))-(Track(i+1,1)-Track(i-1,1)).*(Track(i,2)-Track(i-1,2))) ./ ...
        sqrt(((Track(i,1)-Track(i-1,1)).^2+(Track(i,2)-Track(i-1,2)).^2)*((Track(i+1,1)-Track(i-1,1)).^2+(Track(i+1,2)-Track(i-1,2)).^2)*((Track(i+1,1)-Track(i,1)).^2+(Track(i+1,2)-Track(i,2)).^2));
end
R = 1./K;

% Apexes finden
i = 2;
% ApexIndizes = zeros(1,length(X4_fit));
ApexIndizes = [];

while i < length(Track)-1
    
    if R(i-1) > R(i) && R(i+1) > R(i) && R(i) < 100
        ApexIndizes = [ApexIndizes i];
    end
    
    i = i+1;
end

[Val,Index] = min(R);
figure(1)
plot(Track(:,1),Track(:,2))
hold on
scatter(Track(Index,1),Track(Index,2),40,'g')
scatter(Track(ApexIndizes,1),Track(ApexIndizes,2),20,'r','filled')
%%

% %% Plotten der Strecke
% StartFinish = [Track(1,1) Track(1,2)];  % [x-y] Start/Finish Punkt
% 
% % Plotten der gesamten Strecke
% figure(1)
% plot(Track(:,1),Track(:,2),'k',...
%     StartFinish(1,1),StartFinish(1,2),'o','MarkerSize',6,'MarkerFaceColor','r')
% grid on
% title('Testtrack')
% legend('Track','Start/Finish','Location','northwest')
% axis 'equal'
% xlabel('x [m]')
% ylabel('y [m]')

% % Plotten der Strecke abschnittweise
% figure(2)
% plot(Straight_1(:,1),Straight_1(:,2),'r',...
%     Corner_1(:,1),Corner_1(:,2),'y',...
%     Straight_2(:,1),Straight_2(:,2),'g',...
%     Corner_2(:,1),Corner_2(:,2),'m',...
%     Straight_3(:,1),Straight_3(:,2),'b',...
%     Corner_3(:,1),Corner_3(:,2),'r',...
%     Straight_4(:,1),Straight_4(:,2),'y',...
%     Corner_4(:,1),Corner_4(:,2),'b',...
%     Straight_5(:,1),Straight_5(:,2),'m',...
%     Corner_5(:,1),Corner_5(:,2),'g',...
%     StartFinish(1,1),StartFinish(1,2),'o','MarkerSize',6,'MarkerFaceColor','k')
% hold on
% grid on
% title('Testtrack')
% % legend('Mainstraight','Turn 1','Straight 2','Turn 2','Straight 3','Turn 3',...
% %     'Straight 4','Turn 4','Straight 5','Turn 5','Start/Finish','Location','northwest')
% axis equal
% xlabel('x [m]')
% ylabel('y [m]')