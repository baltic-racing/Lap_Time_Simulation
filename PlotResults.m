%% PlotResults.m
% Plots Result file data with given plotID.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function PlotResults(resultFile1,plotID,runID,saveID,resultFile2)
%PlotResults 
%   Plots the defined result graph with the result data in an extra plot
    switch nargin
        case 5            % Compare two result files
            
            % Loads first and second result file
            result1 = load(resultFile1, '-mat');
            result2 = load(resultFile2, '-mat');      

            % Checks if the track file is 4 or 5 columns and read it acordingly          
            [rows, columns] = size(result1.Track);
            
            if (columns == 4)
                x_Track = result1.Track(:,1);           % [m] X-Koordinate der Strecke
                y_Track = result1.Track(:,2);           % [m] Y-Koordinate der Strecke
                z_Track = [];                           % [m] Z-Koordinate der Strecke
                s = result1.Track(:,3);                 % [m] Verlauf der Streckenlänge
                R = result1.Track(:,4);                 % [m] Kurvenradien
            else
                x_Track = result1.Track(:,1);           % [m] X-Koordinate der Strecke
                y_Track = result1.Track(:,2);           % [m] Y-Koordinate der Strecke
                z_Track = result1.Track(:,3);           % [m] Z-Koordinate der Strecke
                s = result1.Track(:,4);                 % [m] Verlauf der Streckenlänge
                R = result1.Track(:,5);                 % [m] Kurvenradien
            end   
                 
            % Checks if the track file is 4 or 5 columns and read it acordingly          
            [rows, columns] = size(result2.Track);
            
            if (columns == 4)
                x_Track_2 = result2.Track(:,1);           % [m] X-Koordinate der Strecke
                y_Track_2 = result2.Track(:,2);           % [m] Y-Koordinate der Strecke
                z_Track_2 = [];                           % [m] Z-Koordinate der Strecke
                s_2 = result2.Track(:,3);                 % [m] Verlauf der Streckenlänge
                R_2 = result2.Track(:,4);                 % [m] Kurvenradien
            else
                x_Track_2 = result2.Track(:,1);           % [m] X-Koordinate der Strecke
                y_Track_2 = result2.Track(:,2);           % [m] Y-Koordinate der Strecke
                z_Track_2 = result2.Track(:,3);           % [m] Z-Koordinate der Strecke
                s_2 = result2.Track(:,4);                 % [m] Verlauf der Streckenlänge
                R_2 = result2.Track(:,5);                 % [m] Kurvenradien
            end   
            
            t_ges = result1.t(end) *(22000/s(end));
            t_ges_2 = result2.t(end) *(22000/s_2(end));
            
            switch plotID
                
                case 1
                    % Draw RPM Plot
                    figure(1)
                    ni = result1.ni(:,1);
                    ni(end+1) = result1.ni(end);
                    
                    ni_2 = result2.ni(:,1);
                    ni_2(end+1) = result2.ni(end);
                    
                    %figure(1)
                    plot(s,ni,'b',s_2,ni_2,'r')                
                    title('RPM','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Revolutions [1/min]','FontSize',10)
                    legend('RPM 1. Run','RPM 2. Run')
                    grid on
                    box on
                    
                case 2
                    % Draw Aero Plot
                    figure(2)
                    Faero = result1.Faero(:,1);
                    Faero(end+1,1) = result1.Faero(end,1);
                    
                    Faero_2 = result2.Faero(:,1);
                    Faero_2(end+1,1) = result2.Faero(end,1);
                    
                    plot(s,Faero,'b',s_2,Faero_2,'r')                
                    title('Downforce','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Downforce [N]','FontSize',10)
                    legend('Downforce 1. Run','Downforce 2. Run')        
                    grid on
                    box on
                    
                case 3
                    % Draw Torque Plot
                    figure(3)
                    Mi = result1.Mi(:,1);
                    Mi(end+1,1) = result1.Mi(end,1);
                    
                    Mi_2 = result2.Mi(:,1);
                    Mi_2(end+1,1) = result2.Mi(end,1);
                    
                    plot(s,Mi,'b',s_2,Mi_2,'r')                
                    title('Torque','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Torque [Nm]','FontSize',10)
                    legend('Torque 1. Run','Torque 2. Run')
                    grid on
                    box on
                    
                case 4
                    % Draw Power Plot
                    figure(4)
                    P_M = result1.P_M(:,1);
                    P_M(end+1,1) = result1.P_M(end,1);
                    
                    P_M_2 = result2.P_M(:,1);
                    P_M_2(end+1,1) = result2.P_M(end,1);
                    
                    %figure(1)
                    plot(s,P_M,'b',s_2,P_M_2,'r')                
                    title('Power','FontSize',12)
                    xlabel('Track [m]','FontSize',10) 
                    ylabel('Power [W]','FontSize',10)
                    legend('Power 1. Run','Power 2. Run')
                    grid on
                    box on
                    
                case 5
                    % Draw Drag Plot
                    figure(5)
                    FL = result1.FL(:,1);
                    FL(end+1,1) = result1.FL(end,1);
                    
                    FL_2 = result2.FL(:,1);
                    FL_2(end+1,1) = result2.FL(end,1);
                    
                    %figure(1)
                    plot(s,FL,'b',s_2,FL_2,'r')                
                    title('Drag','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Drag [N]','FontSize',10)
                    legend('Drag 1. Run','Drag 2. Run')
                    grid on
                    box on
                    
                case 6
                    % Draw Driving resistance Plot
                    figure(6)
                    Fdr = result1.Fdr(:,1);
                    Fdr(end+1,1) = result1.Fdr(end,1);
                    
                    Fdr_2 = result2.Fdr(:,1);
                    Fdr_2(end+1,1) = result2.Fdr(end,1);
                    
                    plot(s,Fdr,'b',s_2,Fdr_2,'r')                
                    title('Driving Resistance','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Driving Resistance [N]','FontSize',10)
                    legend('Driving resistance 1. Run','Driving resistance 2. Run')
                    grid on
                    box on
                    
                case 7                    
                    % Draw Velocity m/s
                    figure(7)
                    vV1 = result1.vV(:,1);
                    vV2 = result2.vV(:,1);
                    plot(s,vV1, s, vV2)              
                    title('Velocity m/s','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Velocity [m/s]','FontSize',10)
                    legend('Velocity 1. Run','Velocity 2. Run')
                    grid on
                    box on
                    
                case 8                    
                    % Draw Velocity km/h
                    figure(8)
                    vV1 = result1.vV(:,1);
                    vV2 = result2.vV(:,1);
                    plot(s, vV1*3.6, s, vV2*3.6)              
                    title('Velocity km/h','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Velocity [km/h]','FontSize',10)
                    legend('Velocity 1. Run','Velocity 2. Run')
                    grid on
                    box on
                    
                case 9
                    % Draw Cellcapacity Plot
                    figure(9)
                    Capacity_Cellpack = result1.Capacity_Cellpack(:,1);
                    % Capacity_Cellpack(end) = result1.Capacity_Cellpack(end);
                    
                    Capacity_Cellpack_2 = result2.Capacity_Cellpack(:,1);
                    % Capacity_Cellpack_2(end) = result2.Capacity_Cellpack(end);
                    
                    plot(s,Capacity_Cellpack,'b',s_2,Capacity_Cellpack_2,'r')                
                    title('Capacity','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Cellcapacity [As]','FontSize',10)
                    legend('Cellcapacity 1. Run', 'Cellcapacity 2. Run')
                    grid on
                    box on
                    
                case 10
                    % Draw Motor efficiency Plot
                    figure(10)
                    motor_eff = result1.motor_eff(:,1);
                    % motor_eff(end+1) = result1.motor_eff(end);
                    
                    motor_eff_2 = result2.motor_eff(:,1);
                    % motor_eff_2(end+1) = result2.motor_eff(end);
                    
                    %figure(10)
                    plot(s,motor_eff,'b',s_2,motor_eff_2,'r')                
                    title('Motor efficiency','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('efficiency [%]','FontSize',10)
                    legend('Efficiency 1. Run','Efficiency 2. Run')
                    grid on
                    box on
                    
                case 11
                    % Draw Acceleration m/s² Plot
                    figure(11)
                    aVX = result1.aVX(:,1);
                    aVX(end+1) = result1.aVX(end,1);
                    
                    aVY = result1.aVY(:,1);
                    aVY(end+1) = result1.aVY(end,1);
                    
                    aVX_2 = result2.aVX(:,1);
                    aVX_2(end+1) = result2.aVX(end,1);
                    
                    aVY_2 = result2.aVY(:,1);
                    aVY_2(end+1) = result2.aVY(end,1);
                    
                    plot(s,aVX,'b',s,aVY,'b--',s_2,aVX_2,'g',s_2,aVY_2,'g--')                
                    title('Accelerations','FontSize',12)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Acceleration [m/s²]','FontSize',10)
                    legend('Acceleration 1. Run','Acceleration 2. Run')
                    grid on
                    box on
            end
            
        case 4
            
            % Load result file
            result1 = load(resultFile1, '-mat');
            
            % Checks if the track file is 4 or 5 columns and read it acordingly          
            [rows, columns] = size(result1.Track);
            
            if (columns == 4)
                x_Track = result1.Track(:,1);           % [m] X-Koordinate der Strecke
                y_Track = result1.Track(:,2);           % [m] Y-Koordinate der Strecke
                z_Track = [];                           % [m] Z-Koordinate der Strecke
                s = result1.Track(:,3);                 % [m] Verlauf der Streckenlänge
                R = result1.Track(:,4);                 % [m] Kurvenradien
            else
                x_Track = result1.Track(:,1);           % [m] X-Koordinate der Strecke
                y_Track = result1.Track(:,2);           % [m] Y-Koordinate der Strecke
                z_Track = result1.Track(:,3);           % [m] Z-Koordinate der Strecke
                s = result1.Track(:,4);                 % [m] Verlauf der Streckenlänge
                R = result1.Track(:,5);                 % [m] Kurvenradien
            end   
            

            t_ges = result1.t(end) *(22000/s(end));
            
            switch plotID        
                
                case 1
                    % Cell current around track (Zellenstrom über Weg)
                    figure(saveID*10000 + runID*100 + 1)
                    Capacity_Cellpack = result1.Capacity_Cellpack(:,runID);
                    plot(s(end-1),Capacity_Cellpack)
                    title('Cellcapacity','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Capacity [As]','FontSize',10)
                    grid on
                    box on       
                    
                case 2
                    % Slip ----- needs update
                    figure(saveID*10000 + runID*100 + 2)
                    RutschenY_v = result1.RutschenY_v(:,runID);
                    RutschenY_h = result1.RutschenY_h(:,runID);
                    plot(s,RutschenY_v,s,RutschenY_h)
                    hold on
                    ApexIndizes = result1.ApexIndizes(:, runID);
                    scatter(s(ApexIndizes),ones(1,length(ApexIndizes)),20,'r','filled')
                    
                case 3
                    % Power Map on track (Leistungsverlauf auf Strecke)
                    figure(saveID*10000 + runID*100 + 3)
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,result1.P_M(:,runID)/1000,'filled')
                    title('Power Map for P_M','FontSize',12)            %% Leistungsverlauf
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Total consumed motor power [kW]') %%Gesamt-Motorleistung         
                    box on
                    axis equal
                    
                case 4
                    % RPM Plot
                    figure(saveID*10000 + runID*100 + 4)
                    ni = result1.ni(:,runID);
                    ni(end+1,1) = result1.ni(end,1);
                    plot(s,ni)
                    title('RPM','FontSize',12)             %%Drehzahl
                    xlabel('Track Length [m]','FontSize',10)    %%Strecke
                    ylabel('Revolutions [1/min]','FontSize',10) %%Umdrehungen
                    grid on
                    box on
                    
                case 5
                    % Draw track and brake signal
                    figure(saveID*10000 + runID*100 + 5)
                    BPPsignal = result1.BPPsignal(:,runID);
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,BPPsignal,'filled')
                    title('Brake Status','FontSize',12) 
                    xlabel('X-Koordinate [m]','FontSize',10)
                    ylabel('Y-Koordinate [m]','FontSize',10)
                    legend('No brake applied','FontSize',12) 
                    colormap([0 1 0;1 0 0])
                    box on
                    axis equal
                    
                case 6
                    % Tire load front
                    figure(saveID*10000 + runID*100 + 6)
                    FWZ_vl = result1.FWZ_vl(:,runID);
                    FWZ_vr = result1.FWZ_vr(:,runID);
                    plot(s,FWZ_vl,s,FWZ_vr)
                    title('Tireloads front')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Tireload [N]','FontSize',10)
                    legend('Tireload front left','Tirelod front right')
                    grid on
                    box on
                    
                case 7
                    % Tire load rear
                    figure(saveID*10000 + runID*100 + 7)
                    FWZ_hl = result1.FWZ_hl(:,runID);
                    FWZ_hr = result1.FWZ_hr(:,runID);
                    plot(s,FWZ_hl,s,FWZ_hr)
                    title('Tireload rear')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Tireload [N]','FontSize',10)
                    legend('Tireload rear left','Tireload rear right')
                    grid on
                    grid minor
                    box on
                    
                case 8
                    % Axle Load
                    figure(saveID*10000 + runID*100 + 8)
                    FWZh = result1.FWZh(:,runID);
                    FWZv = result1.FWZv(:,runID);
                    plot(s,FWZh,s,FWZv)
                    title('Axle Loads')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Axle load [N]','FontSize',10)
                    legend('axle load rear','axle load front')
                    grid on
                    box on
                    
                case 9
                    % Dynamic Tire Radius
                    figure(saveID*10000 + runID*100 + 9)
                    Rdyn_vl = result1.Rdyn_vl(:,runID);
                    Rdyn_vr = result1.Rdyn_vr(:,runID);
                    Rdyn_hl = result1.Rdyn_hl(:,runID);
                    Rdyn_hr = result1.Rdyn_hr(:,runID);
                    plot(s,Rdyn_vl,s,Rdyn_vr,s,Rdyn_hl,s,Rdyn_hr)
                    title('Dynamic Tire Radius','FontSize',12) 
                    legend('Front left','Front right','Rear left','Rear right')
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Dynamic Tire Radius [m]','FontSize',10)
                    grid on
                    box on
                    grid on
                    
                case 10
                    % Maximal übertragbare und tatsächliche Querkräfte Vorderachse
                    figure(saveID*10000 + runID*100 + 10)
                    FWYmax_v = result1.FWYmax_v(:,runID);
                    FWYv = result1.FWYv(:,runID);
                    plot(s,FWYmax_v,s(1:end-1),abs(FWYv))
                    title('Maximal übertragbare und tatsächliche Querkraft Vorderachse','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Querkraft [N]','FontSize',10)
                    legend('Max. übertragbare Querkraft','Tatsächliche Querkraft')
                    grid on
                    grid minor
                    box on
                    
                case 11
                    % Maximal übertragbare und tatsächliche Querkräfte Hinterachse
                    figure(saveID*10000 + runID*100 + 11)
                    FWYmax_h = result1.FWYmax_h(:,runID);
                    FWYh = result1.FWYh(:,runID);
                    plot(s,FWYmax_h,s(1:end-1),abs(FWYh))
                    title('Maximal übertragbare und tatsächliche Querkraft Hinterachse','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Querkraft [N]','FontSize',10)
                    legend('Max. übertragbare Querkraft','Tatsächliche Querkraft')
                    grid on
                    grid minor
                    box on
                    
                case 12
                    % Verlauf der übertragbaren Querkräfte (Conti-Plot)   
                    
                    TIRparam = loadTIR('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');
                    
                    figure(saveID*10000 + runID*100 + 12)
                    
                    for f = 2000:-200:0
                    plot([-12:0.01:12],MF52_Fy_ps([-12:0.01:12],f,result1.GAMMA,0,TIRparam))
                    hold on
                    end
                    
                    legend('2000 N','1800 N','1600 N','1400 N','1200 N','1000 N','800 N','600 N','400 N','200 N','0 N','FontSize',15)
                    grid on
                    box on
                    
                case 13
                    % Tractive Power Rear wheels
                    figure(saveID*10000 + runID*100 + 13)
                    FVX_hl = result1.FVX_hl(:,runID);
                    FVX_hr = result1.FVX_hr(:,runID);
                    plot(s(1:end-1),FVX_hl,s(1:end-1),FVX_hr)
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Tractive Force rear wheels [N]','FontSize',10)
                    
                    vV = result1.vV(:,runID);
                    yyaxis right
                    plot(s,vV*3.6);
                    ylabel('Velocity [km/h]','FontSize',10)
                    legend('RL','RR','vV [km/h]')
                    grid on
                    box on
                    
                case 14
                    % Traction limit
                    figure(saveID*10000 + runID*100 + 14)
                    % plot(s(1:end-1),FVX,s,FWXmax_h)
                    % legend('Zugkraft','Traktionsgrenze')
                    TC = result1.TC(:,runID);
                    scatter(result1.Track(:,1),result1.Track(:,2),40,TC,'filled')
                    title('Tractioncontrol usage','FontSize',12) 
                    xlabel('X-coordinate [m]','FontSize',10)
                    ylabel('Y-coordinate [m]','FontSize',10)
                    legend('TC on','FontSize',12) 
                    colormap([0 1 0;1 0 0])
                    box on
                    axis equal
                    
                case 15 
                    % Plot Track
                    figure(saveID*10000 + runID*100 + 15)
                    
                    for z = 1:24
                       [Val,Index50(z)] = min(abs(z*50-s)); 
                    end

                    plot(x_Track,y_Track,'k')
                    hold on
                    scatter(x_Track(1),y_Track(1),10,'r','filled')
                    text(x_Track(1),y_Track(1),' Start/Ziel')
                    title('Endurance-track Formula Student Germany','FontSize',12) 
                    for z = 1:24
                       scatter(x_Track(Index50(z)),y_Track(Index50(z)),20,'b','filled')
                       hold on
                       text(x_Track(Index50(z)),y_Track(Index50(z)),['  ' num2str(z*50) 'm'])
                    end
                    annotation('textbox',[0.6, 0.8, 0.25, 0.05],'String',...
                        ['Track Length: ' num2str(max(s)) 'm'],'FitBoxToText','on')
                    hold off
                    axis equal
                    box on
                    
                case 16
                    % Consumed energy by battery over the course of the track (Verbrauchte Akku-Energie über Strecke)
                    figure(saveID*10000 + runID*100 + 16)
                    E_Akku = result1.E_Akku(:,runID);
                    E_Akku_Reku = result1.E_Akku_Reku(:,runID);
                    E_Waerme = result1.E_Waerme(:,runID);
                    E_ges = result1.E_ges(:,runID);
                    plot(s,E_Akku,s,E_Akku_Reku,s,E_Akku-E_Akku_Reku,s,E_Waerme)
                    title('Battery Energy','FontSize',12) 
                    legend('Consumed','Recuperated','Resultant','Thermal') 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Energy [kWh]','FontSize',10)
                    annotation('textbox', [0.2, 0.8, 0.1, 0.1],...
                        'String', {['Consumed Energy = ' num2str(round(E_Akku(end),2)) ' kWh'],...
                        ['Recuperated Energy = ' num2str(round(E_Akku_Reku(end),2)) ' kWh'],...
                        ['Thermal Energy = ' num2str(round(E_Waerme(end),4)) ' kWh'],...
                        ['Resultant Energy = ' num2str(round(E_ges(end),2)) ' kWh']})
                    grid on
                    box on
                    
                case 17
                    % g-g-V Plot  (Not working as planned)
                    figure(saveID*10000 + runID*100 + 17)
                    aVX = result1.aVX(:,runID);
                    aVX(end+1,runID) = result1.aVX(end,runID);
                    
                    aVY = result1.aVY(:,runID);
                    aVY(end+1,runID) = result1.aVY(end,runID);
                    
                    vV = result1.vV(:,runID);
                    figure(17)
                    plot3(aVX,aVY,vV)
                    
                case 18
                    % Velcoity Plot m/s
                    figure(saveID*10000 + runID*100 + 18)
                    vV = result1.vV(:,runID);
                    %%vV(end+1,1) = result1.vV(end,1);
                    plot(s,vV)
                    title('Velocity m/s','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Velocity [m/s]','FontSize',10)
                    grid on
                    grid minor
                    box on
                    
                case 19
                    % Velcoity Plot km/h
                    figure(saveID*10000 + runID*100 + 19)
                    vV = result1.vV(:,runID);
                    plot(s,vV*3.6)
                    title('Velocity km/h','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Velocity [km/h]','FontSize',10)
                    grid on
                    grid minor
                    box on
                    
                case 20
                    % Aero Tire load
                    figure(saveID*10000 + runID*100 + 20)
                    
                    dFWZhl_aero = result1.dFWZhl_aero(:,runID);
                    dFWZhl_aero(end+1,runID) = result1.dFWZhl_aero(end,runID);
                    
                    dFWZhr_aero = result1.dFWZhr_aero(:,runID);
                    dFWZhr_aero(end+1,runID) = result1.dFWZhr_aero(end,runID);
                    
                    dFWZvl_aero = result1.dFWZvl_aero(:,runID);
                    dFWZvl_aero(end+1,runID) = result1.dFWZvl_aero(end,runID);
                    
                    dFWZvr_aero = result1.dFWZvr_aero(:,runID);
                    dFWZvr_aero(end+1,runID) = result1.dFWZvr_aero(end,runID);
                    
%                     % Aero Load all 4 Tires
%                     plot(s, dFWZhl_aero, s, dFWZhr_aero, s, dFWZvl_aero, s, dFWZvr_aero)
%                     legend('HL','HR','VL','VR')

                    % Aero Load both axes + all
                    plot(s, 2*dFWZhl_aero,'b', s, 2*dFWZvl_aero,'r', s, 2*(dFWZhl_aero+dFWZvl_aero),'m')
                    
                    aero_pv = result1.aero_pv(:,runID);
                    aero_ph = result1.aero_ph(:,runID);
                    annotation('textbox', [0.2, 0.8, 0.1, 0.1],...
                        'String', {['Aero Front = ' num2str(aero_pv)],...
                        ['Aero Rear = ' num2str(aero_ph)]})
                    
                    title('Aero Tire Load','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Tire Load [N]','FontSize',10)
                    legend('Rear','Front','Overall')
                    grid on
                    grid minor
                    box on
                    
                case 21
                    % Motor Efficency
                    figure(saveID*10000 + runID*100 + 21)
                    motor_eff = result1.motor_eff(:,runID);
                    motor_eff(end,runID) = motor_eff(end-1,runID);
                    
                    plot(s,motor_eff*100)
                    title('Motor Efficency','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Efficency [%]','FontSize',10)
                    
                case 22
                    % Accelearions m/s
                    figure(saveID*10000 + runID*100 + 22)
                    aVX = result1.aVX(:,runID);
                    aVX(end+1,runID) = result1.aVX(end,runID);
                    
                    aVY = result1.aVY(:,runID);
                    aVY(end+1,runID) = result1.aVY(end,runID);
                    
                    plot(s, aVX, s, aVY);
                    
                    title('longitudinal and Lateral Acceleration','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Acceleration [m/s²]','FontSize',10)
                    legend('AVX','AVY')
                    
                case 23
                    % P_el : Energy consumed including tractive power
                    figure(saveID*10000 + runID*100 + 23)
                    P_el = result1.P_el(:,runID);
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,P_el/1000,'filled')
                    title('Power Map for P_el','FontSize',12)            %%Leistungsverlauf
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Total consumed motor power [kW]') %%Gesamt-Motorleistung         
                    box on
                    axis equal
                    
                 case 24
                    % Consumed energy by battery over the course of the track (Verbrauchte Akku-Energie über Strecke)
                    figure(saveID*10000 + runID*100 + 24)
                    P_el = result1.P_el(:,runID); 
                    P_el(end+1,runID) = result1.P_el(end,runID);
                    plot(s,P_el)
                    
                    t = result1.t(:,runID); 
                    t(1) = [];
                    
                    t_2 = result1.t(:,runID);
                    t_2(end) = [];
                    
                    t_3 = t - t_2;

                    P_el_mean = sum(result1.P_el(:,runID) .* t_3)/t(end); 
                    
                     annotation('textbox', [0.2, 0.8, 0.1, 0.1],...
                        'String', ['Mean Power = ' num2str(P_el_mean) ' W'])
                    
                    title('Power consumed','FontSize',12)   %% P_el including tractive power
                    legend('P_el') 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Energy [kWh]','FontSize',10)
                    grid on
                    box on
                    
                case 25
                    % Lateral Acceleration on track G
                    figure(saveID*10000 + runID*100 + 25)
                    aVY = result1.aVY(:,runID);
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,aVY/9.81,'filled')
                    title('Lateral Acceleration [g]','FontSize',12) 
                    xlabel('X-Koordinate [m]','FontSize',10)
                    ylabel('Y-Koordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Lateral Acceleration [g]')     
                    box on
                    axis equal
                    
                case 26
                    % Longitudinal Acceleration on track G
                    figure(saveID*10000 + runID*100 + 26)
                    aVX = result1.aVX(:,runID);
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,aVX/9.81,'filled')
                    title('Longitudinal Acceleration [g]','FontSize',12) 
                    xlabel('X-Koordinate [m]','FontSize',10)
                    ylabel('Y-Koordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Longitudinal Acceleration [g]')     
                    box on
                    axis equal
                    
                case 27
                    % Accelearions g
                    figure(saveID*10000 + runID*100 + 27)
                    aVX = result1.aVX(:,runID);
                    aVX(end+1,runID) = result1.aVX(end,runID);
                    
                    aVY = result1.aVY(:,runID);
                    aVY(end+1,runID) = result1.aVY(end,runID);
                    
                    plot(s, aVX/9.81, s, aVY/9.81);
                    
                    title('longitudinal and Lateral Acceleration','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Acceleration [g]','FontSize',10)
                    legend('AVX','AVY')
                    
                case 28
                    % Accelearions x/y plot
                    figure(saveID*10000 + runID*100 + 28)
                    aVX = result1.aVX(:,runID);
                    aVX(end+1,runID) = result1.aVX(end,runID);
                    
                    aVY = result1.aVY(:,runID);
                    aVY(end+1,runID) = result1.aVY(end,runID);
                    
                    plot(aVX/9.81, aVY/9.81);
                    
                    title('Longitudinal and Lateral Acceleration','FontSize',12) 
                    xlabel('Longitudinal Acceleration [g]','FontSize',10) 
                    ylabel('Lateral Acceleration [g]','FontSize',10)
                    legend('AVX','AVY')
                    
                case 29
                    % Virtual Current cellpack
                    figure(saveID*10000 + runID*100 + 29)
                    Energy_Cellpack = result1.Energy_Cellpack(:,runID);
                    VirtualCurrent_Cellpack = result1.VirtualCurrent_Cellpack(:,runID);
                    plot(s,Energy_Cellpack,s,VirtualCurrent_Cellpack)
                    title('Current flow','FontSize',12) 
                    legend('Energy Cellpack','VirtualCurrent Cellpack')
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Current [A]','FontSize',10)
                    grid on
                    box on    
                    
                case 30
                    % Accumulator Voltage
                    figure(saveID*10000 + runID*100 + 30)
                    V_i = result1.V_i(:,runID);
                    plot(s,V_i)
                    title('Accumulator voltage','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Capacity [As]','FontSize',10)
                    grid on
                    box on  
                    
                case 31
                    % P_tractive
                    figure(saveID*10000 + runID*100 + 31)
                    P_tractive = result1.P_tractive(:,runID); %%
                    P_tractive(end+1,runID) = result1.P_tractive(end,runID);
                    plot(s,P_tractive)
                    title('Power consumed','FontSize',12)   %% P_el including tractive power
                    legend('P_tractive') 
                    xlabel('Track Length [m]','FontSize',10) 
                    ylabel('Energy [kWh]','FontSize',10)
                    grid on
                    box on
                    
                case 32
                    % Velocity Map
                    figure(saveID*10000 + runID*100 + 32)
                    vV = result1.vV(:,runID);
                    scatter(result1.Track(1:end,1),result1.Track(1:end,2),40,vV,'filled')
                    title('Velocity Map','FontSize',12)            %%Beschleunigungsverlauf
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Velocity [m/s]') %%Velocity        
                    box on
                    axis equal
                    
                case 33
                    % DRS Status Map
                    figure(saveID*10000 + runID*100 + 33)
                    DRS_status = result1.DRS_status(:,runID);
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,DRS_status,'filled')
                    title('DRS Map','FontSize',12)            %%Beschleunigungsverlauf
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap([1 0 0;0 1 0])
                    caxis([0 1])
                    legend('DRS closed','FontSize',12) 
                    box on
                    axis equal
                    
                case 34
                    % Track Analyzer
                    figure(saveID*10000 + runID*100 + 34)
                    
                    apex = (zeros(1,length(result1.ApexIndizes)));
                    
                    % Creates an array of every apex aon the track
                    k = 1;
                    for i = 1:length(R)
                        if ismember(i,result1.ApexIndizes)
                            apex(k,1) = result1.Track(i,1);
                            apex(k,2) = result1.Track(i,2);
                            k = k + 1;
                        end          
                    end
                    
                    % Scatters the track
                    scatter(result1.Track(1:end,1),result1.Track(1:end,2),40,R,'filled')
                    
                    annotation('textbox', [0.6, 0.8, 0.1, 0.1],...
                        'String', {['minimum R = ' num2str(min(R)) ' m'],...
                        ['maximum R = ' num2str(max(R)) ' m'],...             
                        ['num of apexes = ' num2str(length(apex))]})  
                    
                    title('Track Analyzer','FontSize',12)           
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet)
                    caxis([0 4.5]) % Sets the minimum and maximum value for the colorbar 
                    c = colorbar;
                    ylabel(c,'Corner Radius') %%Velocity        
                    
                    % Hold plot to scatter the apexes of the track
                    hold
                    
                    try % Try catch is needed because of tracks without apexes (Accel)
                        scatter(apex(1:end,1),apex(1:end,2),80) % Draws points for every apex on the track
                    catch
                    end
                    
                    box on
                    axis equal
                    
                  case 35
                    % RPM map on track (Drehzahlverlauf auf Strecke)
                    figure(saveID*10000 + runID*100 + 35)
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,result1.ni(:,runID),'filled')
                    title('RPM Map','FontSize',12)            % Drehzahlverlauf
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Total RPM [1/min]') % Overall-RPM      
                    box on
                    axis equal
                    
                 case 36
                    % Torque map on track (Momentenverlauf auf Strecke)
                    figure(saveID*10000 + runID*100 + 36)
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,result1.Mi(:,runID),'filled')
                    title('Torque Map','FontSize',12)            % Momentenverlauf
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet)
                    c = colorbar;
                    ylabel(c,'Total torque [Nm]') % Overll-Torque        
                    box on
                    axis equal
                    
                  case 37
                    % torque Plot
                    figure(saveID*10000 + runID*100 + 37)
                    Mi = result1.Mi(:,runID);
                    Mi(end+1) = result1.Mi(end,runID);
                    plot(s,Mi)
                    title('torque','FontSize',12)               % torque
                    xlabel('Track Length [m]','FontSize',10)    % track
                    ylabel('Torque [Nm]','FontSize',10) % Newtonmeters
                    grid on
                    box on
                    
                case 38
                    % Gear selection on track
                    figure(saveID*10000 + runID*100 + 38)
                    scatter(result1.Track(1:end-1,1),result1.Track(1:end-1,2),40,result1.gear(:,runID),'filled')
                    title('Gearselection','FontSize',12) 
                    xlabel('X-Coordinate [m]','FontSize',10)
                    ylabel('Y-Coordinate [m]','FontSize',10)
                    colormap(jet(length(result1.i_param(:,runID))))
                    c = colorbar('Ticks',[1:1:length(result1.i_param(:,runID))]);
                    caxis([1 length(result1.i_param(:,runID))])
                    ylabel(c,'Gearselection')
                    box on
                    axis equal
                    
                case 39
                    % Gear selection Plot
                    figure(saveID*10000 + runID*100 + 39)
                    plot(s(1:end-1),result1.gear(:,runID),'b','LineWidth',1.5) 
                    title('Gearselection','FontSize',12) 
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Gang [-]','FontSize',10)
                    grid on
                    box on
                    
                case 40
                    % Steering Angle front
                    figure(saveID*10000 + runID*100 + 40)
                    plot(s(1:end),result1.delta(:,runID)*180/pi)
                    title('Steering Angle front')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Steering Angle','FontSize',10)
                    grid on 
                    grid minor
                    
                case 41
                    % Sideslip angle Plot
                    figure(saveID*10000 + runID*100 + 41)
                    plot(s(1:end),result1.beta(:,runID)*180/pi)
                    title('Sideslip angle')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Sideslip Angle','FontSize',10)
                    grid on 
                    grid minor
                    
                case 42
                    % Slip angle
                    figure(saveID*10000 + runID*100 + 42)
                    plot(s(1:end),result1.alpha_f(:,runID),s(1:end),result1.alpha_r(:,runID))
                    legend('Vorne','Hinten')
                    title('Slip Angle')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Slip Angle','FontSize',10)
                    grid on 
                    grid minor
                    
                case 43
                    % vV vs vWoBrake
                    figure(saveID*10000 + runID*100 + 43)
                    
                    vVlength = length(result1.vV);
                    Counter = 0;
                    
                    for i=1:vVlength
                        Diff = result1.vV(i)-result1.vWoBrake(i)
 
                        if Diff > 0
                            vDiff(i) = Diff;
                            Counter = Counter + 1;
                        else
                            vDiff(i) = 0;
                        end
                    end
                    
                    BrakeError = sum(vDiff) / vVlength;
                    ErrorPoints = Counter / vVlength * 100;
                    
                    plot(s(1:end),result1.vV(:,runID),s(1:end),result1.vWoBrake(:,runID),s(1:end),vDiff)
                    legend('vV','vWoBrake','vError')
                    title('Braking Velocity')
                    xlabel('Track Length [m]','FontSize',10)
                    ylabel('Velocity [m/s]','FontSize',10)
                    grid on 
                    grid minor
                                   
                    annotation('textbox', [0.2, 0.8, 0.1, 0.1],...
                        'String', {['Brake Error = ' num2str(BrakeError) ' m/s (' num2str(BrakeError * 3.6) ' km/h)'],...
                        ['Error Points = ' num2str(ErrorPoints) ' %']})
                    
                otherwise
                    return
            end
            
        otherwise
            return
    end

end


