function Vehiclesim_Endurance_GUI_Version2(file, path, trackID)

    global lr lf
    
    % Loads setup file containg all car parameters needed for the
    % simulation
    setup = load(file, '-mat');
    
    % Starts timer
    tic
    
    switch trackID
        case 1
            load('EnduranceTrack.mat','Track');     % Track Data of Endurance Track (Streckendaten Endurance Track)
            
            track.x_Track = Track(:,1);           % [m] X-Coordinate of the Track
            track.y_Track = Track(:,2);           % [m] Y-Coordinate of the Track
            track.z_Track = Track(:,3);           % [m] Z-Coordinate of the Track
            
            track.s = Track(:,4);                 % [m] Track Pathway (Verlauf der Streckenlänge)
            track.R = Track(:,5);                 % [m] Radius of curves (Kurvenradien)
        case 2
            load('AutoXTrack.mat','Track');         % Track Data of AutoX Track (Streckendaten AutoX Track)
            
            track.x_Track = Track(:,1);           % [m] X-Coordinate of the Track
            track.y_Track = Track(:,2);           % [m] Y-Coordinate of the Track
            track.s = Track(:,3);                 % [m] Track Pathway (Verlauf der Streckenlänge)
            track.R = Track(:,4);                 % [m] Radius of curves (Kurvenradien)
        case 3
            load('Testtrack_Data.mat','Track');     % Track Data of Test Track (Streckendaten Test Track)
            
            track.x_Track = Track(:,1);           % [m] X-Coordinate der Strecke
            track.y_Track = Track(:,2);           % [m] Y-Coordinate der Strecke
            track.s = Track(:,3);                 % [m] Track Pathway (Verlauf der Streckenlänge)
            track.R = Track(:,4);                 % [m] Radius of curves (Kurvenradien)
        case 4
            track.s(1) = 0;
            track.x_Track(1) = 0;
            track.y_Track(1) = 0;
            track.R(1) = 0;
            track(1) = 0;
            
            for i = 2:300
                track.y_Track(i) = 0;
                track.z_Track(i) = 0;
                track.s(i) = s(i-1) + 0.25;
                track.x_Track (i) = s(i);
                track.R(i) = 0;
                track(i) = 0; 
            end     
    end
   
    track.ApexIndizes = Apexes(abs(track.R));   % [-] Indices of the Apexes  
    
    %% Initialisation of all variables
    load('Emraxefficiencydata_interpolated.mat');
    load('CorrectedDischargeInterpolated.mat');
    load('RandomizedCellData.mat');
    load('CellparametersVoltageInterpolation.mat');
    
    aRev = zeros(1,length(track)-1);
    aVX = zeros(1,length(track)-1);
    aVY = zeros(1,length(track)-1);
    BPPsignal = zeros(1,length(track)-1);
    cZ_fl = zeros(1,length(track));
    cZ_fr = zeros(1,length(track));
    cZ_rl = zeros(1,length(track));
    cZ_rr = zeros(1,length(track));
    dFWZrl_aero = zeros(1,length(track)-1);
    dFWZrr_aero = zeros(1,length(track)-1);
    dFWZfl_aero = zeros(1,length(track)-1);
    dFWZfr_aero = zeros(1,length(track)-1);
    dFWZrl_x = zeros(1,length(track)-1);
    dFWZrr_x = zeros(1,length(track)-1);
    dFWZfl_x = zeros(1,length(track)-1);
    dFWZfr_x = zeros(1,length(track)-1);
    dFWZrl_y = zeros(1,length(track)-1);
    dFWZrr_y = zeros(1,length(track)-1);
    dFWZvl_y = zeros(1,length(track)-1);
    dFWZvr_y = zeros(1,length(track)-1);
    E_Akku = zeros(1,length(track));
    E_Waerme = zeros(1,length(track));
    E_Akku_Reku = zeros(1,length(track));
    Faero = zeros(1,length(track)-1);
    FB = setup.FB*ones(1,length(track));
    FVX = zeros(1,length(track)-1);
    FVX_rl = zeros(1,length(track)-1);
    FVX_rr = zeros(1,length(track)-1);
    FVXre = zeros(1,length(track)-1);
    FVXid = zeros(1,length(track)-1);
    FR = zeros(1,length(track)-1);
    FL = zeros(1,length(track)-1);
    Fdr = zeros(1,length(track)-1);
    FVY = zeros(1,length(track)-1);
    FWXmax_f = zeros(1,length(track));
    FWXmax_r = zeros(1,length(track));
    FWXmax_fl = zeros(1,length(track));
    FWXmax_fr = zeros(1,length(track));
    FWXmax_rl = zeros(1,length(track));
    FWXmax_rr = zeros(1,length(track));
    FWYf = zeros(1,length(track)-1);
    FWYr = zeros(1,length(track)-1);
    FWYmax_f = zeros(1,length(track));
    FWYmax_r = zeros(1,length(track));
    FWYmax_fl = zeros(1,length(track));
    FWYmax_fr = zeros(1,length(track));
    FWYmax_rl = zeros(1,length(track));
    FWYmax_rr = zeros(1,length(track));
    FWZges = zeros(1,length(track));
    FWZ_fl = zeros(1,length(track));
    FWZ_fr = zeros(1,length(track));
    FWZ_rl = zeros(1,length(track));
    FWZ_rr = zeros(1,length(track));
    FWZr = zeros(1,length(track)-1);
    FWZf = zeros(1,length(track)-1);
    Mi = zeros(1,length(track)-1);
    M_tractive = zeros(1,length(track)-1);
    ni = zeros(1,length(track)-1);
    P_M = zeros(1,length(track)-1);
    P_el = zeros(1,length(track)-1);
    P_tractive = zeros(1,length(track)-1);
    Rdyn_vl = zeros(1,length(track));
    Rdyn_vr = zeros(1,length(track));
    Rdyn_rl = zeros(1,length(track));
    Rdyn_rr = zeros(1,length(track));
    RutschenY_v = zeros(1,length(track));
    RutschenY_h = zeros(1,length(track));
    t = zeros(1,length(track)-1);
    TC = zeros(1,length(track));
    Tirelimit = zeros(1,length(track));
    vAPEXmax = zeros(1,length(track.ApexIndizes));
    vV = zeros(1,length(track)-1);
    vRev = zeros(1,length(track)-1);
    vVYmax = zeros(1,length(track)-1);
    A_Akkuzelle = zeros(1,length(track));
    P_Mloss = zeros(1, length(track));
    P_Bh = zeros(1, length(track));
    motor_eff = zeros(1, length(track));
    Capacity_Cellpack = Parrallelcellgroups(1:131,1);
    SOC_Cellpack(1:131,1) = 1;
    Voltage_Cellpack(1:131,1) = 4.2;
    V_i = zeros(1,length(track));
    VirtualCurrent_Cellpack = zeros(1,length(track));
    Current_Cellpack_Pointer = zeros(1,length(track));
    Energy_Cellpack = zeros(1,length(track));
    Energy_Cellpack_Total = zeros(1,length(track));
    
    %% Tire Model - Magic Tire Formula Model 5.2 (Reifenmodell - Magic Tire Formula Model 5.2)
    
    load('VerticalStiffness_65kPA_IA0.mat');    % Vertical tire stiffness Lookup-Table (Lookup-Table für vertikale Reifensteifigkeit)

    % Variable for path of the tir-file (Variable für Pfad des Tir-Files)
    FileNameLocation = ('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

    % Load tir-file in structure (Tir-File in Struktur laden)
    setup.TIRparam = loadTIR(FileNameLocation);
    R0 = setup.TIRparam.UNLOADED_RADIUS;      % [m] Tire radius - Manufacturing (Fertigungsradius des Reifens)
    bW = setup.TIRparam.WIDTH;                % [m] Tire width (contact area) (Breite des Reifens (Aufstandsbreite))
    
    % Additional factors for longitudinal/lateral slip interaction (Zusätzliche Faktoren für Wechselwirkung Längs/Querschlupf)
    setup.TIRparam.LXAL = 1.2;    % [-] Influence of slip angle on transmissible longitudinal force (Einfluss Schräglaufwinkel auf übertragbare Längskraft)
    setup.TIRparam.LYKA = 1.2;    % [-] Influence of longitudinal slip on transmissible lateral force (Einfluss Längsschlupf auf übertragbare Querkraft)
    setup.TIRparam.RBX3 = 0;      % [-] Additional factor Fx for combined slip (Zusätzlicher Faktor für combined slip Fx)

    % Axle and wheel loads (static) ((Statische) Achs- und Radlasten)
    FG = setup.m_ges*setup.g;                   % [N] Force due to weight of the vehicle
    FWZges(1) = FG;                 % [N] Static total axle load (Statische Gesamtachslast)
    FWZr(1) = setup.m_ph/100*FWZges(1);   % [N] Static rear axle load (Statische Achslast hinten)
    FWZf(1) = FWZges(1)-FWZr(1);    % [N] Static front axle load (Statische Achslast vorne)
    FWZ_fr(1) = FWZf(1)/2;          % [N] Static front right wheel load (Statische Radlast vorne rechts)  
    FWZ_fl(1) = FWZf(1)/2;          % [N] Static front left wheel load (Statische Radlast vorne links)
    FWZ_rr(1) = FWZr(1)/2;          % [N] Static rear right wheel load (Statische Radlast hinten rechts)
    FWZ_rl(1) = FWZr(1)/2;          % [N] Static rear left wheel load (Statische Radlast hinten links)
    
    % Interpolated vertical tire stiffness
    cZ_fl = interp1(Fz,cZ_tire,FWZ_fl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
    cZ_fr = interp1(Fz,cZ_tire,FWZ_fr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
    cZ_rl = interp1(Fz,cZ_tire,FWZ_rl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
    cZ_rr = interp1(Fz,cZ_tire,FWZ_rr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)  
    
    % Dynamic tire radii (stationary = static)
    Rdyn_fl = R0 - FWZ_fl/cZ_fl;    % [m] Dynamic front left tire radius (Dynamischer Reifenradius)
    Rdyn_fr = R0 - FWZ_fr/cZ_fr;    % [m] Dynamic front right tire radius (Dynamischer Reifenradius)
    Rdyn_rl = R0 - FWZ_rl/cZ_rl;    % [m] Dynamic rear left tire radius (Dynamischer Reifenradius)
    Rdyn_rr = R0 - FWZ_rr/cZ_rr;    % [m] Dynamic rear right tire radius (Dynamischer Reifenradius)     
    
    FWXmax_r(1) = Inf;
    
    % Static tire loads for calculation (Statische Reifenlasten für Berechnung)
    setup.FWZ_fl_stat = FWZ_fl(1); 
    setup.FWZ_fr_stat = FWZ_fr(1);
    setup.FWZ_rl_stat = FWZ_rl(1);
    setup.FWZ_rr_stat = FWZ_rr(1);
    
    lf = setup.wheelbase*setup.m_ph/100;     % [mm] Distance from front axle to CoG (Abstand Vorderachse zu Fahrzeugschwerpunkt)
    lr = setup.wheelbase-lf;           % [mm] Distance from rear axle to CoG (Abstand Hinterachse zu Fahrzeugschwerpunkt)

    calculateApexSpeeds(setup);
    
    simulationWithoutBrakes(track, setup);
    
    simulationWithBrakes(setup);
    
%     outputData();
    
    writeSaveFile();
end   

function writeSaveFile()
    %% Writing the results to Mat File (Schreiben der Ergebnisse in Mat File)
    result.tEnd = tEnd;
    result.t_ges = t_ges;    
    result.Track = Track;
    result.aRev = aRev;
    result.aVX = aVX;
    result.aVY = aVY;
    result.BPPsignal = BPPsignal;
    result.cZ_vl = cZ_vl;
    result.cZ_vr = cZ_vr;
    result.cZ_hl = cZ_hl;
    result.cZ_hr = cZ_hr;
    result.dFWZhl_aero = dFWZhl_aero;
    result.dFWZhr_aero = dFWZhr_aero;
    result.dFWZvl_aero = dFWZvl_aero;
    result.dFWZvr_aero = dFWZvr_aero;
    result.dFWZhl_x = dFWZhl_x;
    result.dFWZhr_x = dFWZhr_x;
    result.dFWZvl_x = dFWZvl_x;
    result.dFWZvr_x = dFWZvr_x;
    result.dFWZhl_y = dFWZhl_y;
    result.dFWZhr_y = dFWZhr_y;
    result.dFWZvl_y = dFWZvl_y;
    result.dFWZvr_y = dFWZvr_y;
    result.E_Akku = E_Akku;
    result.E_Waerme = E_Waerme;
    result.E_Akku_Reku = E_Akku_Reku;
    result.E_ges = E_ges;
    result.Faero = Faero;
    %result.FB = FB*ones(1,length(Track));
    result.FVX = FVX;
    result.FVX_hl = FVX_hl;
    result.FVX_hr = FVX_hr;
    result.FVXre = FVXre;
    result.FVXid = FVXid;
    result.FR = FR;
    result.FL = FL;
    result.Fdr = Fdr;
    result.FVY = FVY;
    result.FWXmax_v = FWXmax_v;
    result.FWXmax_h = FWXmax_h;
    result.FWXmax_vl = FWXmax_vl;
    result.FWXmax_vr = FWXmax_vr;
    result.FWXmax_hl = FWXmax_hl;
    result.FWXmax_hr = FWXmax_hr;
    result.FWYv = FWYv;
    result.FWYh = FWYh;
    result.FWYmax_v = FWYmax_v;
    result.FWYmax_h = FWYmax_h;
    result.FWYmax_vl = FWYmax_vl;
    result.FWYmax_vr = FWYmax_vr;
    result.FWYmax_hl = FWYmax_hl;
    result.FWYmax_hr = FWYmax_hr;
    result.FWZges = FWZges;
    result.FWZ_hl = FWZ_hl;
    result.FWZ_hr = FWZ_hr;
    result.FWZ_vl = FWZ_vl;
    result.FWZ_vr = FWZ_vr;
    result.FWZh = FWZh;
    result.FWZv = FWZv;
    result.Mi = Mi;
    result.ni = ni;
    result.P_M = P_M;
    result.P_el = P_el;
    result.Rdyn_vl = Rdyn_vl;
    result.Rdyn_vr = Rdyn_vr;
    result.Rdyn_hl = Rdyn_hl;
    result.Rdyn_hr = Rdyn_hr;
    result.RutschenY_v = RutschenY_v;
    result.RutschenY_h = RutschenY_h;
    result.t = t;
    result.TC = TC;
    result.Tirelimit = Tirelimit;
    result.vAPEXmax = vAPEXmax;
    result.vV = vV;
    result.vRev = vRev;
    result.vVYmax = vVYmax;
    result.A_Akkuzelle = A_Akkuzelle;
    result.ApexIndizes = ApexIndizes;
%     result.M_eff_inter = M_eff_inter;
    result.motor_eff = motor_eff;
    
    result.Energy_Cellpack = Energy_Cellpack;
    result.VirtualCurrent_Cellpack = VirtualCurrent_Cellpack;
    result.V_i = V_i;
    result.Capacity_Cellpack = Capacity_Cellpack;
    
    % Aero
    result.aero_ph = aero_ph;
    result.aero_pv = aero_pv;
    
    % Car Parameters needed to draw result plots
    result.P_max = max_power;
    result.GAMMA = GAMMA;
    
    
    [~, name, ~] = fileparts(file);
    
    savefilename = name + "_result.mat";
    
    save(savefilename, '-struct','result');
    
    disp('File succesfully written');
end

function outputData()
    %%  Output of the values
    disp(['Endurance Time for ONE Lap: ' num2str(t(end)) ' s = ' num2str(t(end)/60) ' min']);
    t_ges = t(end) *(22000/s(end));      
    disp(['Endurance Total Time: ' num2str(t_ges(end)) ' s = ' num2str(t_ges(end)/60) ' min']);
    disp(['Total Travel Distance: ' num2str(s(end)) ' m']);
    E_Akku_gesamt_ohne_reku = E_Akku(end) * (22000/s(end));
    disp(['Total energy consumption W/O recuperator: ' num2str(E_Akku_gesamt_ohne_reku) ' kWh']);
    E_Akku_gesamt = E_ges(end) * (22000/s(end));
    disp(['Total energy consumption W recuperator: ' num2str(E_Akku_gesamt) ' kWh']);
    tEnd = toc;
end

function calculateApexSpeeds(setup)
    global ApexIndizes R lr lf
    
    %% Calculation of the maximum apex speed for all apexes (numerically) (Berechnen der maximalen Kurvengeschwindigkeiten für alle Apexes (numerisch))
    for i = 1:length(ApexIndizes)

        FWYf(i) = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
        FWYr(i) = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
        FWYmax_f(i) = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal übertragbare Querkraft Vorderachse)
        FWYmax_r(i) = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal übertragbare Querkraft Hinterachse)
        vV(i) = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)

        while  FWYf(i) < FWYmax_f(i) && FWYr(i) < FWYmax_r(i) && vV(i) < 30

            vV(i) = vV(i) + 0.01;   % [m/s] Increaing vehicle speed (Erhöhen der Fahrzeuggeschwindigkeit)
            
            FVY(i) = setup.m_ges*vV(i)^2/R(ApexIndizes(i));    % [N] Centrifugal force (Zentrifugalkraft)

            aVY(i) = vV(i)^2/R(i);  % [m/s²] Lateral acceleration (Querbeschleunigung)

            % Lateral forces to be applied on front and rear axle (Aufzubringende Querkräfte an Vorder- und Hinterachse)
            FWYf(i) = lr/setup.wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
            FWYr(i) = lf/setup.wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

            % Wheel load transfer due to aerodynamic forces         
            [Faero(i), dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = aeroForceWheels(setup, vV(i));

            % calculates Dynamic wheel load displacement and wheel loads
            [~, ~, FWYmax_f(i), FWYmax_r(i)] = wheelload_transfer(setup, 0, FVY(i), dFWZfl_aero(i), dFWZfr_aero(i), dFWZrl_aero(i), dFWZrr_aero(i));                     
        end

        vAPEXmax(i) = vV(i);   % [m/s] Maximum speed for any apex (Maximalgeschwindigkeit für jede Apex)
    end
    
    disp('caclculated Apex Speeds!')
end

function simulationWithoutBrakes(setup, track)
    %% Start/Initial values for first simulation run WITHOUT BRAKES 
    vV(1) = 15;    % [m/s] Speed 
    t(1) = 0;      % [s] Time

    % Supporting variables
    z = 1;         % [-] Determination of the upcoming apex

    %% Simulation WITHOUT BRAKES (Simulation OHNE BREMSEN)
    NonBrakeApexes = [];
    for i = 1:length(track)-1

        % Determination of motor speed and gear (Bestimmen von Motordrehzahl und Gang)
        ni(i) = vV(i)*30/pi*i_G/Rdyn_hl(i); % [1/min] Determine current motor rpm
       
        Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Interpolated motor torque 
          
        wheelPower();

        FR = k_R*FWZges;                % [N] Rolling resistance (Rollwiderstand)
        FL = rho_L*vV^2/2*c_w*A_S;      % [N] Air Resistance/Drag (Luftwiderstand)
        Fdr = FR + FL;                  % [N] Total resistance (Gesamtwiderstand)

        FVY = m_ges*vV^2/R;             % [N] Centrifugal force (Zentrifugalkraft)
        aVX = (FVX-Fdr-FB)/m_ges;       % [m/s²] Longitudinal acceleration (Längsbeschleunigung)
        aVY = vV^2/R;                   % [m/s²] Lateral acceleration (Querbeschleunigung)
        
        if ismember(i,ApexIndizes)
            if vV(i) > vAPEXmax(z)   % Limiting curve maximum speeds at apexes (Begrenzen auf maximale Kurvengeschwindigkeit in Apexes)
                vV(i) = vAPEXmax(z);
            end
            z = z + 1;
        end
        
        vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle seed (Gesamt-Fahrzeuggeschwindigkeit)
        t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Time (Zeit)

        % Lateral forces to be applied on front and rear axles (Aufzubringende Querkräfte an Vorder- und Hinterachse)
        FWYv(i) = lh/l*FVY(i);   % [N] Lateral force to be applied on front axle (Aufzubringende Querkraft der Vorderachse)
        FWYh(i) = lv/l*FVY(i);   % [N] Lateral force to be applied on rear axle (Aufzubringende Querkraft der Hinterachse)

        % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokräften)
        
        aeroForceWheels(setup, vV(i));

        % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerungen in Längsrichtung)
        
        [dFWZvl_x(i), dFWZvr_x(i), dFWZhl_x(i), dFWZhr_x(i)] = wheelload_longdisp(h_COG, m_ges, aVX(i), l);
        
        % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
        
        [dFWZvl_y(i), dFWZvr_y(i), dFWZhl_y(i), dFWZhr_y(i)] = wheelload_latdisp(h_COG, B, lh, lv, l, FVY(i));
        
        % Wheel loads (Radlasten)
        FWZ_vl(i+1) = FWZ_vl(1) + dFWZvl_aero(i) + dFWZvl_x(i) + dFWZvl_y(i); % [N] Front left wheel load (Radlast vorne links)
        FWZ_vr(i+1) = FWZ_vr(1) + dFWZvr_aero(i) + dFWZvr_x(i) + dFWZvr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
        FWZ_hl(i+1) = FWZ_hl(1) + dFWZhl_aero(i) + dFWZhl_x(i) + dFWZhl_y(i); % [N] Rear left wheel load (Radlast hinten links)
        FWZ_hr(i+1) = FWZ_hr(1) + dFWZhr_aero(i) + dFWZhr_x(i) + dFWZhr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)
        
        % Limiting the wheel loads to (almost) zero (Begrenzen der Radlasten auf (quasi) Null)
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

        % Axle loads - for dynamic radii (Achslasten)
        
        [FWZh(i+1), FWZv(i+1), FWZges(i+1)] = axleloads(FWZ_hl(i+1), FWZ_hr(i+1), FWZ_vl(i+1), FWZ_vr(i+1));
        
        % Vertical tire stiffnesses - for dynamic radii (Vertikale Reifensteifigkeiten)
       
        [cZ_vl(i+1), cZ_vr(i+1), cZ_hl(i+1), cZ_hr(i+1)] = vtirestiff(Fz, cZ_tire, FWZ_vl(i+1), FWZ_vr(i+1), FWZ_hl(i+1), FWZ_hr(i+1));
        
        % Dynamic tire radii (Dynamische Reifenradien)
        
        [Rdyn_vl(i+1), Rdyn_vr(i+1), Rdyn_hl(i+1), Rdyn_hr(i+1)] = dyn_radii(R0, FWZ_vl(i+1), FWZ_vr(i+1), FWZ_hl(i+1), FWZ_hr(i+1), cZ_vl(i+1), cZ_vr(i+1), cZ_hr(i+1), cZ_hl(i+1));
        
        % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)
      
        [FWXmax_vl(i+1), FWXmax_vr(i+1), FWXmax_hl(i+1), FWXmax_hr(i+1), FWXmax_v(i+1), FWXmax_h(i+1)] = longi_tireforces(FWZ_vl(i+1), FWZ_vr(i+1),FWZ_hl(i+1), FWZ_hr(i+1), GAMMA, TIRparam);
        
        % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)
       
        [FWYmax_vl(i+1), FWYmax_vr(i+1), FWYmax_hl(i+1), FWYmax_hr(i+1), FWYmax_v(i+1), FWYmax_h(i+1)] = lat_tireforces(FWZ_vl(i), FWZ_vr(i),FWZ_hl(i), FWZ_hr(i), GAMMA, TIRparam);
        
    end
    
    disp('Simulated without brakes!')

    % Temporary saving of various variables for plots (Zwischenspeichern von diversen Variablen für Plots)
    vV_NoBrake = vV;    

    %% BRAKING POINT CALCULATION (BREMSPUNKTBERECHNUNG)
    BrakeIndizes = [];

    k = length(track.ApexIndizes);

    while k >= 1

        Counter = 0;
        j = track.ApexIndizes(k);
        vRev(j-1) = vAPEXmax(k);

        while vRev(j-1) < vV(j-1)
%             Faero(j-1) = interp1(v,FA,vRev(j-1)*3.6,'linear','extrap'); % [N] Interpolated downforce (Abtriebskraft interpoliert)
          
          [Faero(i)] = aeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i)); % [N] Aerodynamic force
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

        if k > 1 && j < ApexIndizes(k-1)
            NonBrakeApexes = [NonBrakeApexes k-1];
            k = k - 1;
        end

        k = k - 1;

    end
end

function simulationWithBrakes()
 %% Start values for simulation WITH BRAKES (Startwerte für Simulation MIT BREMSEN)
    vV(1) = 15;     % [m/s] Speed (Geschwindigkeit)
    t(1) = 0;       % [s] Time (Zeit)

    E_Akku(1) = 0;       % [J] Energy consumed by battery (Verbrauchte Energie Akku)
    E_Waerme(1) = 0; 
    E_Akku_Reku(i) = 0;  % [J] Energy recuperated by battery (Rekuperierte Energie Akku)

    % Supporting variables (Hilfsgrößen)
    z = 1;          % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

    %% SIMULATION WITH BRAKES (SIMULATION MIT BREMSEN)
    for i = 1:length(Track)-1

        % Checking at which apex the vehicle is (Überprüfen, vor welcher Apex das Auto ist)
        if ismember(i,ApexIndizes)  
            z = z + 1;
        end
        
        % Motor RPM (Motordrehzahl)
        ni(i) = vV(i)*30/pi*i_G/Rdyn_hl(i); % [1/min]     

        % Checking if braking is required (Prüfen, ob gebremst werden muss)
        if ismember(i,BrakeIndizes)                       % Initiaion of braking process (Einleiten des Bremsvorgangs)     
            Mi(i) = 0;                                    % [Nm] Motor torque (Motormoment)
            BPPsignal(i) = 1;                             % [-] Brake signal (Bremssignal)
            P_Bh(i) = FWZh(i)/FWZges(i)*FB(i)*vV(i);      % [W] Rear braking power for recuperation (Bremsleistung hinten für Rekuperation)
        else
            Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motor torque (Motormoment)
            FB(i) = 0;                                    % [N] Braking force (Bremskraft)
            P_Bh(i) = 0;                                  % [W] Rear braking power (Bremsleistung hinten)
        end

        % Motor power & limitation to 80 kW from FS-Rules (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
        drehzahlpointer = round(ni(i));                   % Pointer for efficiency table (Pointer für effizienztabelle)
        
        if Mi(i) <= 0
            Momentenpointer = 1;
        else
            Momentenpointer = round(Mi(i));               % Pointer for efficiency table (Pointer für effizienztabelle)
        end

        P_M(i) = num_motors * Mi(i) * ni(i) / 60 * 2 * pi;% [W] Total motor power (Gesamt-Motorleistung)
        if P_M(i) > max_power
            P_M(i) = max_power;                           % [W] Limited power (Begrenzte Leistung)
            Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Begrenzen des Moments)
        end
        
        if(drehzahlpointer > setup.n_max)
            drehzahlpointer = setup.n_max;
        end

        % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
        motor_eff(i) = M_eff_inter(drehzahlpointer,Momentenpointer);
        
        P_Mloss(i) = P_M(i)*(1-(motor_eff(i)*eta_ges*eta_inv));% Calculation of power loss (berechnung der Verlustleistung)
        
        %P_el(i) = P_M(i);
        
        P_M(i) = P_M(i) * eta_ges * M_eff_inter(i) * eta_inv;  % Calculation of motor power after deduction of 0.95 efficiency for inverter (berechnung der Motorleistung nach abzug der Effizienz  0.95 für inverter)
        Mi(i) = P_M(i)*(60/ni(i)/2/pi);                       % some issue dont understand
        
        % Calculation of tractive forces (Berechnen der Zugkraft)
        FVX_hl(i) = Mi(i)*i_G/Rdyn_hl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
        FVX_hr(i) = Mi(i)*i_G/Rdyn_hr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
        FVX(i) = FVX_hl(i) + FVX_hr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

        if FVX(i) > FWXmax_h(i)     % Limiting the tractive force to the traction limit (Begrenzen der Zugkraft auf Traktionsgrenze)
            FVX(i) = FWXmax_h(i);
            TC(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
        end
        
        M_tractive(i) = FVX(i)/(i_G/Rdyn_hr(i));            % [Nm] Torque including tractive force
        P_tractive(i) = M_tractive(i)/(60/ni(i)/2/pi);      % [kW] Motor power required for traction 
        P_el(i) = (P_tractive(i)/(eta_ges * M_eff_inter(i) * eta_inv))/2;     % [kW] Motor power including efficiencies
        
        % Driving resistances (Fahrwiderstände) & Vehicle (Fahrzeug)
     
        [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = vehicle_resistances_forces(k_R, FWZges(i), rho_L, vV(i), c_w, A_S, m_ges, R(i), FVX(i), FB(i));
        vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i)));     % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
       
        if ismember(i,BrakeIndizes) && vV(i+1) < vAPEXmax(z) && not(ismember(z,NonBrakeApexes))  % Begrenzen der Geschwindigkeit auf ApexGeschwindigkeit (Bremst solange bis Geschwindigkeiten gleich)
             vV(i+1) = vAPEXmax(z);                         % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
        end

        t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);                % [s] Time (Zeit)
        
        % Battery energy capacity (Energiemenge Akku)
        E_Akku_Reku(i+1) = E_Akku_Reku(i) + (t(i+1)-t(i)) * P_Bh(i);                    % [J]
        E_Akku(i+1) = E_Akku(i) + (t(i+1)-t(i)) * P_el(i);                              % [J] 
        E_Waerme(i+1) = P_el(i)*(1-M_eff_inter(i))+(E_Akku(i) + (t(i+1)-t(i)))*0.05;    % [J] Motor losses + 5% for invertor losses; Drivetrain losses with heat (Motorverluste + 5% flat für inverterverlsute % Drivetrain Verluste auch Wärme)


        % Lateral forces on front and rear axle (Querkräfte an Vorder- und Hinterachse)
        FWYv(i) = lh/l*FVY(i);   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
        FWYh(i) = lv/l*FVY(i);   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

        if FWYv(i) > FWYmax_v(i)
           RutschenY_v(i) = 1;  
        end

         if FWYh(i) > FWYmax_h(i)
           RutschenY_h(i) = 1;  
        end

         aeroForceWheels();

         wheelload_transfer();

         % Limiting the wheel loads to (almost) zero (Begrenzen der Radlasten auf (quasi) Null)
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

        % Axle loads - for dynamic radii (Achslasten)       
        [FWZh(i+1), FWZv(i+1), FWZges(i+1)] = axleloads(FWZ_hl(i+1), FWZ_hr(i+1), FWZ_vl(i+1), FWZ_vr(i+1)); 

         
        % Vertical tire stiffness - for dynamic radii (Vertikale Reifensteifigkeiten)        
        [cZ_vl(i+1), cZ_vr(i+1), cZ_hl(i+1), cZ_hr(i+1)] = vtirestiff(Fz, cZ_tire, FWZ_vl(i+1), FWZ_vr(i+1), FWZ_hl(i+1), FWZ_hr(i+1));
        
        % Dynamic tire radii (Dynamische Reifenradien)
       
        [Rdyn_vl(i+1), Rdyn_vr(i+1), Rdyn_hl(i+1), Rdyn_hr(i+1)] = dyn_radii(R0, FWZ_vl(i+1), FWZ_vr(i+1), FWZ_hl(i+1), FWZ_hr(i+1), cZ_vl(i+1), cZ_vr(i+1), cZ_hr(i+1), cZ_hl(i+1));
        
        % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)
        
        [FWXmax_vl(i+1), FWXmax_vr(i+1), FWXmax_hl(i+1), FWXmax_hr(i+1),FWXmax_v(i+1), FWXmax_h(i+1)] = longi_tireforces(FWZ_vl(i+1), FWZ_vr(i+1),FWZ_hl(i+1), FWZ_hr(i+1), GAMMA, TIRparam);
        
        % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)
        [FWYmax_vl(i+1), FWYmax_vr(i+1), FWYmax_hl(i+1), FWYmax_hr(i+1), FWYmax_v(i+1), FWYmax_h(i+1)] = lat_tireforces(FWZ_vl(i), FWZ_vr(i),FWZ_hl(i), FWZ_hr(i), GAMMA, TIRparam);
        
        
      
        % Battery Currents (Akkuströme)
        A_Akkuzelle(i) = P_el(i) / V_i(i) / nZellen_parralel;
        
        %Akkuströme
        V_i(i) = sum(Voltage_Cellpack(:,i));
        
        Current_Cellpack_Pointer(i) = P_M(i) / V_i(i) *10 ; %Strombelastung eines %er Parrallel Paketes in 0,1A parameter für die berechnung der korrigierten belastung mit höhren verlusten durch höhere zellströme
        if Current_Cellpack_Pointer(i) <= 1
            Current_Cellpack_Pointer(i)=1;
        end

        if Current_Cellpack_Pointer(i) >= 1500 %begrenzen des max Zellstromes auf 30A pro Zelle im 5er parralelverbund also 150A
            Current_Cellpack_Pointer(i)=1500;
        end
        
        VirtualCurrent_Cellpack(i) = CorrectedDischargeInterpolated(1,round(Current_Cellpack_Pointer(i))); %Berechnung der Virtuell höheren zellströme basierend auf den höheren verlsuten durch höhere Ströme
    
        Energy_Cellpack(i) = (VirtualCurrent_Cellpack(i)*(t(i+1)-t(i))) - ((P_Bh(i)/V_i(i))*(t(i+1)-t(i))) ; %Energieverbrauch in As für ein 5erpacket an akkuzellen -> Akkustrom zum zeitpunkt i mal Zeitdifferenz zwischen i und i+1
        Energy_Cellpack_Total(i+1) = Energy_Cellpack_Total(i) + Energy_Cellpack(i); % Über Endurance Run Integrierte Energieverbrauch in As für ein 5erpacket an akkuzellen

        Capacity_Cellpack(1:131,i+1) =  Capacity_Cellpack(1:131,i)- Energy_Cellpack(i); 

        SOC_Cellpack(1:131,i+1) = Capacity_Cellpack(1:131,i)./Capacity_Cellpack(1:131,1); %Berechnung des SOC für den nächsten tick basierend auf der aktuellen cellcapacity und der im nächsten tick

        SOC_Pointer(1:131,i+1) = round(SOC_Cellpack(1:131,i+1)*1000);
        Current_Cellpack_Pointer_Voltage(1,i+1) = round(Current_Cellpack_Pointer(i)/5);

        if Current_Cellpack_Pointer_Voltage(i) <= 3
            Current_Cellpack_Pointer_Voltage(i)=3;
        end

        Voltage_Cellpack(1:131,i+1) = Voltage_inter(Current_Cellpack_Pointer_Voltage(1,i),SOC_Pointer(1:131,i+1));
    end
end

 function wheelPower(i)
        % Motor power & limitation to 80 kW from FS-Rules (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
        drehzahlpointer = round(ni(i));                   % Pointer for efficiency table (Pointer für effizienztabelle)
        
        if Mi(i) <= 0
            Momentenpointer = 1;
        else
            Momentenpointer = round(Mi(i));               % Pointer for efficiency table (Pointer für effizienztabelle)
        end

        P_M(i) = num_motors * Mi(i) * ni(i) / 60 * 2 * pi;% [W] Total motor power (Gesamt-Motorleistung)
        if P_M(i) > max_power
            P_M(i) = max_power;                           % [W] Limited power (Begrenzte Leistung)
            Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Begrenzen des Moments)
        end
        
        if(drehzahlpointer > setup.n_max)
            drehzahlpointer = setup.n_max;
        end

        % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
        motor_eff(i) = M_eff_inter(drehzahlpointer,Momentenpointer);
        
        P_Mloss(i) = P_M(i)*(1-(motor_eff(i)*eta_ges*eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)
        
        %P_el(i) = P_M(i);
        
        P_M(i) = P_M(i) * eta_ges * M_eff_inter(i) * eta_inv;  % Calculation of motor power after deduction of efficiency of the inverter
        Mi(i) = P_M(i)*(60/ni(i)/2/pi);                       
        
        % Calculation of tractive forces (Berechnen der Zugkraft)
        if num_motors == 4
            FVX_fl(i) = Mi(i)*i_G/Rdyn_hl(i);                   % [N] Tractive Force on front left wheel (AWD) 
            FVX_fr(i) = Mi(i)*i_G/Rdyn_hr(i);                   % [N] Tractive Force on front right wheel (AWD) 
            FVX(i) = FVX_hl(i) + FVX_hr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)
            
            FVX_hl(i) = Mi(i)*i_G/Rdyn_hl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
            FVX_hr(i) = Mi(i)*i_G/Rdyn_hr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
            FVX(i) = FVX_hl(i) + FVX_hr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)
            
            if FVX(i) > FWXmax_h(i)     % Limiting the tractive force to the traction limit front axle
                FVX(i) = FWXmax_v(i);
                TC_front(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
            end
        elseif num_motors == 2
            FVX_hl(i) = Mi(i)*i_G/Rdyn_hl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
            FVX_hr(i) = Mi(i)*i_G/Rdyn_hr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
            FVX(i) = FVX_hl(i) + FVX_hr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)
        else
            FVX_hl(i) = Mi(i)/2*i_G/Rdyn_hl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
            FVX_hr(i) = Mi(i)/2*i_G/Rdyn_hr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
            FVX(i) = FVX_hl(i) + FVX_hr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)
        end

        if FVX(i) > FWXmax_h(i)     % Limiting the tractive force to the traction limit rear axle
            FVX(i) = FWXmax_h(i);
            TC(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
        end
end


function [FWXmax_f, FWXmax_r, FWYmax_f, FWYmax_r] = wheelload_transfer(setup, aVX, FVY, dFWZfl_aero, dFWZfr_aero, dFWZrl_aero, dFWZrr_aero)
    % Dynamic wheel load displacement in longitudinal direction
    dFWZfl_x = -setup.m_ges*aVX*setup.h_cog/setup.wheelbase/2;  % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
    dFWZfr_x = -setup.m_ges*aVX*setup.h_cog/setup.wheelbase/2;  % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
    dFWZrl_x = setup.m_ges*aVX*setup.h_cog/setup.wheelbase/2;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
    dFWZrr_x = setup.m_ges*aVX*setup.h_cog/setup.wheelbase/2;   % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
    
    lv = setup.wheelbase*setup.m_ph/100;        % [mm] Distance from front axle to CoG (Abstand Vorderachse zu Fahrzeugschwerpunkt)
    lh = setup.wheelbase-lv;                    % [mm] Distance from rear axle to CoG (Abstand Hinterachse zu Fahrzeugschwerpunkt)
    
    % Dynamic wheel load displacement in lateral direction
    dFWZfl_y = -setup.h_cog/setup.track*lh/setup.wheelbase*FVY;   % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
    dFWZfr_y = setup.h_cog/setup.track*lh/setup.wheelbase*FVY;    % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
    dFWZrl_y = -setup.h_cog/setup.track*lv/setup.wheelbase*FVY;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
    dFWZrr_y = setup.h_cog/setup.track*lv/setup.wheelbase*FVY;    % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
    
    % Wheel loads (Radlasten)
    FWZ_fl = setup.FWZ_fl_stat + dFWZfl_aero + dFWZfl_x + dFWZfl_y; % [N] Front left wheel load (Radlast vorne links)
    FWZ_fr = setup.FWZ_fr_stat + dFWZfr_aero + dFWZfr_x + dFWZfr_y; % [N] Front right wheel load (Radlast vorne rechts)
    FWZ_rl = setup.FWZ_rl_stat + dFWZrl_aero + dFWZrl_x + dFWZrl_y; % [N] Rear left wheel load (Radlast hinten links)
    FWZ_rr = setup.FWZ_rr_stat + dFWZrr_aero + dFWZrr_x + dFWZrr_y; % [N] Rear right wheel load (Radlast hinten rechts)   
    
    % Maximum transmissible tire forces in lateral direction
    FWYmax_fl = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_fl,setup.camber,0,setup.TIRparam))); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
    FWYmax_fr = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_fr,setup.camber,0,setup.TIRparam))); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
    FWYmax_rl = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_rl,setup.camber,0,setup.TIRparam))); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
    FWYmax_rr = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_rr,setup.camber,0,setup.TIRparam))); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)

    FWYmax_f = FWYmax_fl + FWYmax_fr;    % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
    FWYmax_r = FWYmax_rl + FWYmax_rr;    % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
    
    % Maximum transmissible tire forces in longitudinal direction   
    FWXmax_fl = max(abs(MF52_Fx_cs(0,FWZ_fl,setup.camber,0:0.01:0.2,setup.TIRparam))); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
    FWXmax_fr = max(abs(MF52_Fx_cs(0,FWZ_fr,setup.camber,0:0.01:0.2,setup.TIRparam))); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
    FWXmax_rl = max(abs(MF52_Fx_cs(0,FWZ_rl,setup.camber,0:0.01:0.2,setup.TIRparam))); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
    FWXmax_rr = max(abs(MF52_Fx_cs(0,FWZ_rr,setup.camber,0:0.01:0.2,setup.TIRparam))); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)        

    FWXmax_f = FWXmax_fl + FWXmax_fr;                                % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
    FWXmax_r = FWXmax_rl + FWXmax_rr;                                % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
end

function [Faero, dFWZrl_aero, dFWZrr_aero, dFWZfl_aero, dFWZfr_aero] = aeroForceWheels(setup, vV)
    
    rho_L = setup.p_L*10^5/(setup.R_L*(setup.t_L+273.15));

    % Aerodynamic force
    Faero = setup.downforce_multiplier * setup.c_l * setup.A * rho_L * (vV)^2 / 2; 

    dFWZrl_aero = Faero/2*(1-setup.aero_pv);   % [N] Aerodynamic force on rear left wheel (Aerokraft auf linkes Hinterrad)
    dFWZrr_aero = Faero/2*(1-setup.aero_pv);   % [N] Aerodynamic force on rear right wheel (Aerokraft auf rechtes Hinterrad)
    dFWZfl_aero = Faero/2*setup.aero_pv;   % [N] Aerodynamic force on front left wheel (Aerokraft auf linkes Vorderrad)
    dFWZfr_aero = Faero/2*setup.aero_pv;   % [N] Aerodynamic force on front right wheel (Aerokraft auf rechtes Vorderrad)
end