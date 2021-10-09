%% Vehiclesim_Endurance_GUI_Version.m
% Main function of the simulation. 
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function Vehiclesim_Endurance_GUI_Version(setupFile, path, TrackFileName, disciplineID, sensitivityID, minValue, stepSize, numSteps, processDataButtonHandle, textAreaHandle, sensitivityID2, minValue2, stepSize2, logCellData, Debug, StartingSpeed, numOfLaps, brakeFunction)

    % DEBUG
    %brakeFunction = 2;

    % Adds search path of the setup file and then loads the setup.
    setup = load(path+"/"+setupFile, '-mat');

    tic   
    
    %% Initalize GUI loading bar
    % Store original button text
    originalButtonText = processDataButtonHandle.Text;
    % When the function ends, return the original button state
    cleanup = onCleanup(@()set(processDataButtonHandle,'Text',originalButtonText,'Icon',''));
    % Change button name to "Processing"
    processDataButtonHandle.Text = 'Simulating...';
    % Put text on top of icon
    processDataButtonHandle.IconAlignment = 'bottom';
    % Create waitbar with same color as button
    wbar = permute(repmat(processDataButtonHandle.BackgroundColor,15,1,200),[1,3,2]);
    % Black frame around waitbar
    wbar([1,end],:,:) = 0;
    wbar(:,[1,end],:) = 0;
    % Load the empty waitbar to the button
    processDataButtonHandle.Icon = wbar;
    
    %% loads the selected track
    [~, ~, ~, s, R, Track, ApexIndexes, lapLength] = loadTrack(TrackFileName, disciplineID, numOfLaps);
    
    %% Vehicle Data (Fahrzeugdaten)   
    m_tot = setup.m_ges;    % [kg] Total mass of the vehicle including driver (Fahrzeuggesamtmasse inkl. Fahrer)
    g = setup.g;            % [m/s²] Acceleration due to Gravity (Erdbeschleunigung)
    FG = m_tot*g;           % [N] Force due to weight of the vehicle (Gewichtskraft des Fahrzeugs)
    m_balast = setup.m_ballast;
    m_driver = setup.m_driver;
    h_COG_balast = setup.h_cog_ballast;
    h_COG_driver = setup.h_cog_driver;
    
    wheelbase = setup.wheelbase;    % [mm] Wheelbase (Radstand)
    track = setup.track;            % [mm] Trackwidth - Front & Rear (Spurweite für vorne und hinten)
    h_COG = setup.h_cog;            % [mm] Height of vehicle's Center of Gravity(COG) (Höhe Fahrzeugschwerpunkt)
    x_COG = setup.x_cog;            % [mm] x-Coordinate of vehicle's CoG in CAD
    x_vA = setup.x_va;              % [mm] x-Coordinate of the front axle in CAD
%     x_COP = 1600;           % [mm] x-Coordinate of Center of Pressure in CAD
    
%     aero_ph = (l-(l-(x_COP-x_vA)))/l;   % [-] Aerodynamic Force on rear axle (Anteil Aerokraft auf Hinterachse)
%     aero_pv = 1-aero_ph;                % [-] Aerodynamic Force on front axle (Anteil Aerokraft auf Vorderachse)

    % Replaced Formula with direct downforce percentage front
    aero_pv = setup.aero_pv;
    aero_ph = 1-aero_pv;

    thetaV_X = setup.thetaV_X;                          % [kg*m²] Moment of Inertia of vehicle  about X-Axis (Trägheitsmoment Fahrzeug um X-Achse)
    thetaV_Y = setup.thetaV_Y;                          % [kg*m²] Moment of Inertia of vehicle about Y-Axis (Trägheitsmoment Fahrzeug um Y-Achse)
    thetaV_Z = setup.thetaV_Z;                          % [kg*m²] Moment of Inertia of vehicle about Z-Axis (Trägheitsmoment Fahrzeug um Z-Achse)  

    m_ph = setup.m_ph;                                  % [%] Percentage of rear axle wheel load (Prozentualer Radlastanteil Hinterachse)
    lf = wheelbase*m_ph/100;                            % [mm] Distance from front axle to CoG (Abstand Vorderachse zu Fahrzeugschwerpunkt)
    lr = wheelbase-lf;                                  % [mm] Distance from rear axle to CoG (Abstand Hinterachse zu Fahrzeugschwerpunkt)

    FB = setup.FB;                                      % [N] Maximum Braking Force (Maximale Bremskraft)

    k_R = setup.k_R;                                    % [-] Co-efficient of rolling resistance (Rollwiderstandsbeiwert)
    c_w = setup.c_w;                                    % [-] cw value (cw-Wert)
    c_l = setup.c_l;                                    % [-] cl value 
    A_S = setup.A;                                      % [m²] Front Surface Area (Stirnfläche)
    downforce_multiplier = setup.downforce_multiplier;  % [-] multiplier for car downforce
    c_d_DRS = setup.c_d_DRS;                            % [-] Drag coefficient with DRS
                                
    DRS = setup.DRS;                                    % [-]   DRS on/off
    
    try 
        ConstantDownforce = setup.ConstantDownforce;    % [N] Constant Downforce for example from ground effect fans
    catch
        ConstantDownforce = 0;
    end
    
    try
        brakeBias_setup = setup.brakeBias_setup;
    catch
        brakeBias_setup = -1;
    end
    
    %% Initalise DRS
    try
        DRS_Radius = setup.DRS_Radius;                  % [m] DRS 
    catch
        DRS_Radius = 300;
    end

    try
        c_l_DRS = setup.c_l_DRS;                        % [-] Lift coefficient with DRS
    catch
        c_l_DRS = c_l;
    end
    
    for i = 1:length(R)
        if R(i) > DRS_Radius && DRS
            DRS_status(i) = 1;
        else
            DRS_status(i) = 0;
        end

    end      
    
    %% Environmental Conditions (Umgebungsbedingungen)
    t_L = setup.t_L;                        % [°C] Ambient air temperature (Umgebungslufttemperatur)
    p_L = setup.p_L/100000;                 % [bar] Ambient air pressure (Umgebungsluftdruck)
    R_L = setup.R_L;                        % [J/(kg*K)] Air gas constant (Gaskonstante Luft)
    rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m³] Air density (Luftdichte)

    %% Motor Data (Motordaten)
    
    n = setup.engine_param(:,1);
    M = setup.engine_param(:,2);   
    
    num_motors = setup.num_motors;
%     P = num_motors * M.*n*2*pi/60;        % [W] Read power matrix 
%     P_Mmax = max(P);                      % [W] Determination of the Maximum Power 
    n_Mmax = setup.n_max;                   % [1/min] Maximum RPM of the motor 
    drivetrain_eff = setup.drivetrain_eff;  % [-] Overall efficiency of powertrain
    max_power = setup.p_max * 1000;         % [W] Power Limit in endurance 
    trq_multiplier = setup.trq_multiplier;
    
    ptype = setup.ptype;                    % variable for powertrain type : 1 = Electric and 0 = Combustion
    
    if ptype
        eta_inv = setup.invertor_eff;       % [-] Inverter efficiency
    else
        eta_inv = 1;
    end
    
    n_shift = setup.n_shift;                % engine rpm when shifting gears [1/min]
    n_downshift = setup.n_downshift;
    t_shift = setup.t_shift;                % time needed to shift gears [s]
    gr = setup.i_param(:);                  % Gear ratio for each gear  
    gearbox = setup.gearbox;
    t_x = 0;                                % Initialise the shift dead time
    idleRPM = setup.idleRPM;
    
    % Sets gearratio to a constant 1 if no gearbox is present
    if ~gearbox
        gr(1) = 1;
    end
    
    % Gearbox Data (Getriebedaten)
    z_chaindrive = setup.z_chaindrive;      % [-] Number of teeth on sprocket
    z_pinion = setup.z_sprocket;            % [-] Number of teeth on pinion 
    i_P = setup.i_P;
    i_G = z_chaindrive / z_pinion * i_P;          % [-] Gear ratio (Motor to wheel)
    
    % HV-Accumulator Data
%     V_i = 550;                                    % [V] Voltage 
    Energy_i = setup.Energy_i;                      % [kWh] Energy 
    ncells_parallel = setup.nZellen_Parallel;       % [-] Number of parallel cell rows 
    capacity_singlecell = setup.capacity_cell;      % [Wh] Capacity of single cell
    number_cells_in_a_row = setup.nZellen_Reihe;    % [-] Number of cells in one row
    capacity_accumulator = number_cells_in_a_row*ncells_parallel*capacity_singlecell; % [Wh] Capacity of total battery
    
    writeToLogfile('loaded Car Data!', Debug, textAreaHandle);

    %% Initialisation of all variables
    
    load('Emraxefficiencydata_interpolated.mat');
    load('CorrectedDischargeInterpolated.mat');
    load('RandomizedCellData.mat');
    load('CellparametersVoltageInterpolation.mat'); 

    %% Tire Model - Magic Tire Formula Model 5.2 (Reifenmodell - Magic Tire Formula Model 5.2)
    
    load('VerticalStiffness_65kPA_IA0.mat');    % Vertical tire stiffness Lookup-Table (Lookup-Table für vertikale Reifensteifigkeit)

    % Load tir-file in structure (Tir-File in Struktur laden)
    TIRparam = loadTIR('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

    % Tire Data (Reifendaten)
    J_Tire = setup.J_Tire;              % [kg*m²] Moment of inertia of tire about axis of rotaion (Trägheitsmoment des Reifens um Drehachse)
    R0 = TIRparam.UNLOADED_RADIUS;      % [m] Tire radius - Manufacturing (Fertigungsradius des Reifens)
    bW = TIRparam.WIDTH;                % [m] Tire width (contact area) (Breite des Reifens (Aufstandsbreite))
    p_infl = setup.p_Tire;              % [Pa] Tire pressure (Luftdruck des Reifens)
    GAMMA = setup.camber;               % [°] Camber (Sturz)

    % Scaling factors for grip level (Skalierungsfaktoren für Grip-Niveau)
    % 0.75 for optimum tire temperature; 0.6 for low tire temperature (0.75 für optimale Reifentemperatur, 0.6 für niedrige Reifentemperatur)
    TIRparam.LMUX = setup.LMUX;          % [-] Longitudinal scaling factor (Skalierungsfaktor Längsrichtung)
    TIRparam.LMUY = setup.LMUY;          % [-] Lateral scaling factor (Skalierungsfaktor Querrichtung)

    % Additional factors for longitudinal/lateral slip interaction (Zusätzliche Faktoren für Wechselwirkung Längs/Querschlupf)
    TIRparam.LXAL = 1.2;    % [-] Influence of slip angle on transmissible longitudinal force (Einfluss Schräglaufwinkel auf übertragbare Längskraft)
    TIRparam.LYKA = 1.2;    % [-] Influence of longitudinal slip on transmissible lateral force (Einfluss Längsschlupf auf übertragbare Querkraft)
    TIRparam.RBX3 = 0;      % [-] Additional factor Fx for combined slip (Zusätzlicher Faktor für combined slip Fx)

    %% Set progress
    % sets button progress (progressbar)
    currentProg = min(round((size(wbar,2)-2)*(0/numSteps)),size(wbar,2)-2); 
    processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 1) = 0.25391; % (royalblue)
    processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 2) = 0.41016;
    processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 3) = 0.87891;
    drawnow; % updates button progress
    
    %% Sensitivity Analysis
    steps = 0;
    
    % Checks if sensitvity analysis with two variables is started or
    % not, if so the the second loop is activated with the same length
    % as the outer loop.
    if sensitivityID2 ~= 0 
        numSteps2 = numSteps;
    else
        numSteps2 = 1;
    end
    

    %% For used for sensitivity analysis
    for steps1 = 1:numSteps 
        
        switch sensitivityID
            case 1
                A_S = minValue + stepSize*(steps1-1);
            case 2
                FB = minValue + stepSize*(steps1-1);
            case 3
                TIRparam.LMUX = minValue + stepSize*(steps1-1);
            case 4
                TIRparam.LMUY = minValue + stepSize*(steps1-1);
            case 5
                aero_pv = minValue + stepSize*(steps1-1);
            case 6
                c_l = minValue + stepSize*(steps1-1);
            case 7
                c_w = minValue + stepSize*(steps1-1);
            case 8
                GAMMA = minValue + stepSize*(steps1-1);
            case 9
                downforce_multiplier = minValue + stepSize*(steps1-1);
            case 10
                drivetrain_eff = minValue + stepSize*(steps1-1);
            case 11
                m_balast = minValue + stepSize*(steps1-1);
            case 12
                m_driver = minValue + stepSize*(steps1-1);
            case 13
                m_tot = minValue + stepSize*(steps1-1);
            case 14
                m_ph = minValue + stepSize*(steps1-1);
            case 15
                n_Mmax = minValue + stepSize*(steps1-1);
            case 16
                num_motors = minValue + stepSize*(steps1-1);
            case 17
                p_infl = minValue + stepSize*(steps1-1);
            case 18
                max_power = minValue + stepSize*(steps1-1);
            case 19
                t_L = minValue + stepSize*(steps1-1);
                rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m³] Air density (Luftdichte)
            case 20
                thetaV_X = minValue + stepSize*(steps1-1);
            case 21
                thetaV_Y = minValue + stepSize*(steps1-1);
            case 22
                thetaV_Z = minValue + stepSize*(steps1-1);
            case 23
                track = minValue + stepSize*(steps1-1);
            case 24
                trq_multiplier = minValue + stepSize*(steps1-1);
            case 25
                wheelbase = minValue + stepSize*(steps1-1);
            case 26
                h_COG = minValue + stepSize*(steps1-1);
            case 27
                x_COG = minValue + stepSize*(steps1-1);
            case 28
                h_COG_balast = minValue + stepSize*(steps1-1);
            case 29
                h_COG_driver = minValue + stepSize*(steps1-1);
            case 30
                x_vA = minValue + stepSize*(steps1-1);
            case 31
                z_chaindrive = minValue + stepSize*(steps1-1);
                z_pinion = setup.z_sprocket;        % [-] Number of teeth on pinion (Zähnezahl des Ritzels)
                i_G = z_chaindrive / z_pinion * i_P;     % [-] Gear ratio (Motor to wheel) (Übersetzung Motor zu Rad)
            case 32
                z_pinion = minValue + stepSize*(steps1-1);
                z_chaindrive = setup.z_chaindrive; % [-] Number of teeth on sprocket (Zähnezahl des Kettenblatts)
                i_G = z_chaindrive / z_pinion * i_P;     % [-] Gear ratio (Motor to wheel) (Übersetzung Motor zu Rad)
            case 33       
                i_G = minValue + stepSize*(steps1-1);
            case 34
                ConstantDownforce = minValue + stepSize*(steps1-1);
            case 35
                DRS_Radius = minValue + stepSize*(steps1-1);
        end
        
        
        % Second for loop for 2 variable sensitivity analysis
        for steps2 = 1:numSteps2
            steps = steps + 1;
        
           
            if sensitivityID2 ~= 0 
                switch sensitivityID2
                    case 1
                        A_S = minValue2 + stepSize2*(steps2-1);
                    case 2
                        FB = minValue2 + stepSize2*(steps2-1);
                    case 3
                        TIRparam.LMUX = minValue2 + stepSize2*(steps2-1);
                    case 4
                        TIRparam.LMUY = minValue2 + stepSize2*(steps2-1);
                    case 5
                        aero_pv = minValue2 + stepSize2*(steps2-1);
                    case 6
                        c_l = minValue2 + stepSize2*(steps2-1);
                    case 7
                        c_w = minValue2 + stepSize2*(steps2-1);
                    case 8
                        GAMMA = minValue2 + stepSize2*(steps2-1);
                    case 9
                        downforce_multiplier = minValue2 + stepSize2*(steps2-1);
                    case 10
                        drivetrain_eff = minValue2 + stepSize2*(steps2-1);
                    case 11
                        m_balast = minValue2 + stepSize2*(steps2-1);
                    case 12
                        m_driver = minValue2 + stepSize2*(steps2-1);
                    case 13
                        m_tot = minValue2 + stepSize2*(steps2-1);
                    case 14
                        m_ph = minValue2 + stepSize2*(steps2-1);
                    case 15
                        n_Mmax = minValue2 + stepSize2*(steps2-1);
                    case 16
                        num_motors = minValue2 + stepSize2*(steps2-1);
                    case 17
                        p_infl = minValue2 + stepSize2*(steps2-1);
                    case 18
                        max_power = minValue2 + stepSize2*(steps2-1);
                    case 19
                        t_L = minValue2 + stepSize2*(steps2-1);
                        rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m³] Air density (Luftdichte)
                    case 20
                        thetaV_X = minValue2 + stepSize2*(steps2-1);
                    case 21
                        thetaV_Y = minValue2 + stepSize2*(steps2-1);
                    case 22
                        thetaV_Z = minValue2 + stepSize2*(steps2-1);
                    case 23
                        track = minValue2 + stepSize2*(steps2-1);
                    case 24
                        trq_multiplier = minValue2 + stepSize2*(steps2-1);
                    case 25
                        wheelbase = minValue2 + stepSize2*(steps2-1);
                    case 26
                        h_COG = minValue2 + stepSize2*(steps2-1);
                    case 27
                        x_COG = minValue2 + stepSize2*(steps2-1);
                    case 28
                        h_COG_balast = minValue2 + stepSize2*(steps2-1);
                    case 29
                        h_COG_driver = minValue2 + stepSize2*(steps2-1);
                    case 30
                        x_vA = minValue2 + stepSize*(steps2-1);
                    case 31
                        z_chaindrive = minValue2 + stepSize2*(steps2-1);
                        z_pinion = setup.z_sprocket;        % [-] Number of teeth on pinion (Zähnezahl des Ritzels)
                        i_G = z_chaindrive / z_pinion;     % [-] Gear ratio (Motor to wheel) (Übersetzung Motor zu Rad)
                    case 32
                        z_pinion = minValue2 + stepSize2*(steps2-1);
                        z_chaindrive = setup.z_chaindrive; % [-] Number of teeth on sprocket (Zähnezahl des Kettenblatts)
                        i_G = z_chaindrive / z_pinion;     % [-] Gear ratio (Motor to wheel) (Übersetzung Motor zu Rad)
                    case 33       
                        i_G = minValue2 + stepSize2*(steps2-1);
                    case 34
                        ConstantDownforce = minValue2 + stepSize2*(steps2-1);
                    case 35
                        DRS_Radius = minValue2 + stepSize2*(steps2-1);
                end
            end    
            
            % Initalize simulation start values
            [aRev,aVX,aVY,BPPsignal,cZ_fl,cZ_fr,cZ_rl,cZ_rr,dFWZrl_aero,dFWZrr_aero,dFWZfl_aero,dFWZfr_aero,dFWZrl_x,dFWZrr_x,dFWZfl_x,dFWZfr_x...
            ,dFWZrl_y,dFWZrr_y,dFWZfl_y,dFWZfr_y,E_Accu,E_Accu_Recu,E_heat,Faero,FB,FVX,FVX_f,FVX_fr,FVX_fl,FVX_rl,FVX_rr,FVXre,FVXid,FR,FL,Fdr,FVY,...
            FWXmax_f,FWXmax_r,FWXmax_fl,FWXmax_fr,FWXmax_rl,FWXmax_rr,FWYf,FWYr,FWYmax_f,FWYmax_r,FWYmax_fl,FWYmax_fr,FWYmax_rl,FWYmax_rr,FWZtot,FWZ_rl,...
            FWZ_rr,FWZ_fl,FWZ_fr,FWZr,FWZf,Mi,M_tractive,ni,P_M,P_el,P_tractive,Rdyn_fl,Rdyn_fr,Rdyn_rl,Rdyn_rr,slipY_f,slipY_r,t,l_contact_patch_fl,...
            l_contact_patch_fr,l_contact_patch_rl,l_contact_patch_rr,kappa_rl,kappa_rr,kappa_fl,kappa_fr,delta,beta,psi1,alpha_f,alpha_r,TC,TC_front,ABS,...
            gearSelection,Tirelimit,vAPEXmax,vV,vVYmax,A_accu_cell,P_Mloss,P_Bh,motor_eff,Capacity_Cellpack,SOC_Cellpack,Voltage_Cellpack,V_i,VirtualCurrent_Cellpack,...
            Current_Cellpack_Pointer,Energy_Cellpack,Energy_Cellpack_Total,SOC_Pointer,Current_Cellpack_Pointer_Voltage] = initializeStartValues(FB, Track, ApexIndexes);
            
            % Axle and wheel loads (static) ((Statische) Achs- und Radlasten)
            FWZtot(1) = FG;                 % [N] Static total axle load (Statische Gesamtachslast)
            FWZr(1) = m_ph/100*FWZtot(1);   % [N] Static rear axle load (Statische Achslast hinten)
            FWZf(1) = FWZtot(1)-FWZr(1);    % [N] Static front axle load (Statische Achslast vorne)
            FWZ_fr(1) = FWZf(1)/2;          % [N] Static front right wheel load (Statische Radlast vorne rechts)  
            FWZ_fl(1) = FWZf(1)/2;          % [N] Static front left wheel load (Statische Radlast vorne links)
            FWZ_rr(1) = FWZr(1)/2;          % [N] Static rear right wheel load (Statische Radlast hinten rechts)
            FWZ_rl(1) = FWZr(1)/2;          % [N] Static rear left wheel load (Statische Radlast hinten links)

            % Interpolated vertical tire stiffness (Vertikale Reifensteifigkeiten interpoliert) 
            [cZ_fl(1), cZ_fr(1), cZ_rl(1), cZ_rr(1)] = calculateVtirestiff(Fz, cZ_tire, FWZ_fl(1), FWZ_fr(1), FWZ_rl(1), FWZ_rr(1));

            % Dynamic tire radii (stationary = static) (Dynamische Reifenradien (im Stand = statisch))
            [Rdyn_fl(1), Rdyn_fr(1), Rdyn_rl(1), Rdyn_rr(1)] = calculateDynRadii(R0, FWZ_fl(1), FWZ_fr(1), FWZ_rl(1), FWZ_rr(1), cZ_fl(1), cZ_fr(1), cZ_rr(1), cZ_rl(1));

            FWXmax_r(1) = Inf;

            % Static tire loads for calculation (Statische Reifenlasten für Berechnung)
            FWZ_fl_stat = FWZ_fl(1); 
            FWZ_fr_stat = FWZ_fr(1);
            FWZ_rl_stat = FWZ_rl(1);
            FWZ_rr_stat = FWZ_rr(1);
            
            % Calculate Skidpad Time and Speed
            [t_skidpad, vV_skidpad] = calculateSkidPad(downforce_multiplier, c_l, A_S, rho_L, ConstantDownforce, c_l_DRS, DRS_status, m_tot, lr, lf, wheelbase, track, aero_ph, aero_pv, h_COG, GAMMA, TIRparam, FWZ_fl_stat, FWZ_fr_stat, FWZ_rl_stat, FWZ_rr_stat);

            
            % DEBUG
            if c_w == 1 && c_l == 1.6
                AFWawf = 3;
            end
            
            %% Calculation of the maximum apex speed for all apexes (numerically) (Berechnen der maximalen Kurvengeschwindigkeiten für alle Apexes (numerisch))
            for i = 1:length(ApexIndexes)

                FWYf(i) = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
                FWYr(i) = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
                FWYmax_f(i) = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal übertragbare Querkraft Vorderachse)
                FWYmax_r(i) = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal übertragbare Querkraft Hinterachse)
                vV(i) = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)
                
                if R(ApexIndexes(i)) > 0
                    f = -1;
                else
                    f = 1;
                end

                while  FWYf(i) < FWYmax_f(i) && FWYr(i) < FWYmax_r(i) && vV(i) < 40

                    vV(i) = vV(i) + 0.01;   % [m/s] Increaing vehicle speed (Erhöhen der Fahrzeuggeschwindigkeit)
              
                    Faero(i) = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i), ConstantDownforce, c_l_DRS, DRS_status(ApexIndexes(i))); % [N] Aerodynamic force

                    FVY(i) = m_tot*vV(i)^2/R(ApexIndexes(i));    % [N] Centrifugal force (Zentrifugalkraft)

                    %aVY(i) = vV(i)^2/R(ApexIndexes(i));  % [m/s²] Lateral acceleration (Querbeschleunigung)

                    % Lateral forces to be applied on front and rear axle (Aufzubringende Querkräfte an Vorder- und Hinterachse)
                    FWYf(i) = lr/wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                    FWYr(i) = lf/wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)
                    
                    % Lenkwinkel, Schwimmwinkel, Gierrate, Schräglaufwinkel  
                    delta(i) = f*atan((wheelbase/1000)/sqrt(R(ApexIndexes(i))^2-(lr/1000)^2));   % [rad] Lenkwinkel
                    beta(i) = f*atan((lr/1000)/sqrt(R(ApexIndexes(i))^2-(lr/1000)^2));  % [rad] Schwimmwinkel
                    psi1(i) = vV(i)/R(ApexIndexes(i));                               % [rad/s] Gierrate
                    alpha_f(i) = 180/pi*(delta(i)-(lf/1000)/vV(i)*psi1(i)-beta(i));  % [°] Schräglaufwinkel vorne
                    alpha_r(i) = 180/pi*((lr/1000)/vV(i)*psi1(i)-beta(i));           % [°] Schräglaufwinkel hinten

                    % Wheel load transfer due to aero forces (Radlastverlagerung in Folge von Aerokräften) 
                    [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = calculateAeroforceOnWheels(Faero(i), aero_ph, aero_pv);

                    % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in Längsrichtung = 0 angenommen)
                    [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = calculateWheelloadLongDisp(h_COG, 0, aVX(i), wheelbase); % Loads = 0 assumed

                    % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                    [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = calculateWheelloadLatDisp(h_COG, track, lr, lf, wheelbase, FVY(i));

                    % Wheel loads (Radlasten)
                    FWZ_fl(i) = FWZ_fl_stat + dFWZfl_aero(i) + dFWZfl_x(i) + dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                    FWZ_fr(i) = FWZ_fr_stat + dFWZfr_aero(i) + dFWZfr_x(i) + dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                    FWZ_rl(i) = FWZ_rl_stat + dFWZrl_aero(i) + dFWZrl_x(i) + dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                    FWZ_rr(i) = FWZ_rr_stat + dFWZrr_aero(i) + dFWZrr_x(i) + dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)   

                    % Maximum transmissible tire forces in longitudinal direction = 0 assumed (because longitudinal wheel loads = 0 assumed) 

                    % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)    
                    [FWYmax_f(i), FWYmax_r(i)] = calculateLatTireforces(FWZ_fl(i), FWZ_fr(i),FWZ_rl(i), FWZ_rr(i), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

                end

                vAPEXmax(i) = vV(i);   % [m/s] Maximum speed for any apex (Maximalgeschwindigkeit für jede Apex)
            end

            writeToLogfile('caclculated Apex Speeds!', Debug, textAreaHandle);
            
            %% Start/Initial values for first simulation run WITHOUT BRAKES (Startwerte für ersten Simulationslauf OHNE BREMSEN)
            vV(1) = StartingSpeed + 0.005;
            
            gear = 1;

            t(1) = 0;      % [s] Time (Zeit)

            % Supporting variables (Hilfsgrößen)
            z = 1;         % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

            %% Simulation WITHOUT BRAKES (Simulation OHNE BREMSEN)
            for i = 1:length(Track)-1
                
                if i > 1    % if i > 1 use real rpm instead of idle rpm
                    [ni(i), gear, t_x] = calculateGearbox(gearbox, idleRPM, n_shift, n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, n_Mmax, t_x, ni(i-1), t(i), t(i-1));      % Calculates Gearbox data and rpm
                else
                    [ni(i), gear, t_x] = calculateGearbox(gearbox, idleRPM, n_shift, n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, n_Mmax, t_x);      % Calculates Gearbox data and rpm
                end 
                    
                % Determination of aero forces and motor torque (Bestimmen von Aero-Kräften und Motormoment)      
                Faero(i) = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i), ConstantDownforce, c_l_DRS, DRS_status(i)); % [N] Aerodynamic force

                % [Nm] Interpolated motor torque (Motormoment interpoliert)
                Mi(i) = interp1(n,M,ni(i),'linear','extrap'); 

                % Pointer for efficiency table (Pointer für effizienztabelle)
                rpmpointer = round(ni(i));                          

                if Mi(i) <= 0
                    torquepointer = 1;
                else
                    torquepointer = round(Mi(i));               % Pointer for efficiency table (Pointer für effizienztabelle)
                end

                % Motor power & limitation to 80 kW from FS-Rules (for electric cars) (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
                P_M(i) = num_motors * Mi(i) * ni(i) / 60 * 2 * pi;% [W] Total motor power (Gesamt-Motorleistung)
                if ptype && P_M(i) > max_power
                    P_M(i) = max_power;                           % [W] Limited power (Begrenzte Leistung)
                    %Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Begrenzen des Moments)
                end

                if(rpmpointer > n_Mmax)
                    rpmpointer = n_Mmax;
                elseif(rpmpointer < 1)
                    rpmpointer = 1;
                end

                % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
                if ptype
                    motor_eff(i) = M_eff_inter(rpmpointer,torquepointer);
                else
                    motor_eff(i) = 1;
                end

                P_Mloss(i) = P_M(i)*(1-(motor_eff(i)*drivetrain_eff*eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)

                P_M(i) = P_M(i) - P_Mloss(i);  % Calculation of motor power after deduction of efficiency of the inverter
                
                % Calculation Overall Torque with real power
                Mi(i) = P_M(i)*(60/ni(i)/2/pi);  
                
                % Calculate the tractive forces on the wheels
                [FVX_fl(i), FVX_fr(i), FVX_rl(i), FVX_rr(i), FVX(i), FVX_f(i), TC_front(i), TC(i)] = calculateTractiveForces(Mi(i), num_motors, i_G, gr, Rdyn_fl(i), Rdyn_fr(i), Rdyn_rl(i), Rdyn_rr(i), t_x, gear, FWXmax_f(i), FWXmax_r(i), t_shift);

                % Driving resistances (Fahrwiderstände) & Vehicle (Fahrzeug)        
                [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = calculateVehicleResistancesForces(k_R, FWZtot(i), rho_L, vV(i), c_w, A_S, m_tot, R(i), FVX(i), FVX_f(i), c_d_DRS, DRS_status(i), rpmpointer, n_Mmax, 0, 0, 0, 0);

                if ismember(i,ApexIndexes)
                    if vV(i) > vAPEXmax(z)   % Limiting maximum speeds at apexes (Begrenzen auf maximale Kurvengeschwindigkeit in Apexes)
                        vV(i) = vAPEXmax(z);
                    end
                    z = z + 1;
                end

                vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Time (Zeit)

                % Lateral forces to be applied on front and rear axles (Aufzubringende Querkräfte an Vorder- und Hinterachse)
                FWYf(i) = lr/wheelbase*FVY(i);   % [N] Lateral force to be applied on front axle (Aufzubringende Querkraft der Vorderachse)
                FWYr(i) = lf/wheelbase*FVY(i);   % [N] Lateral force to be applied on rear axle (Aufzubringende Querkraft der Hinterachse)
                
                delta(i) = f*atan((wheelbase/1000)/sqrt(R(i)^2-(lr/1000)^2));       % [rad] Lenkwinkel
                beta(i) = f*atan((lr/1000)/sqrt(R(i)^2-(lr/1000)^2));               % [rad] Schwimmwinkel
                psi1(i) = vV(i)/R(i);                                               % [rad/s] Gierrate
                alpha_f(i) = 180/pi*(delta(i)-(lf/1000)/vV(i)*psi1(i)-beta(i));     % [°] Schräglaufwinkel vorne
                alpha_r(i) = 180/pi*((lr/1000)/vV(i)*psi1(i)-beta(i));              % [°] Schräglaufwinkel hinten

                % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokräften)  
                [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = calculateAeroforceOnWheels(Faero(i), aero_ph, aero_pv);

                % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerungen in Längsrichtung)
                [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = calculateWheelloadLongDisp(h_COG, m_tot, aVX(i), wheelbase);

                % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = calculateWheelloadLatDisp(h_COG, track, lr, lf, wheelbase, FVY(i));

                % Wheel loads (Radlasten)
                FWZ_fl(i+1) = FWZ_fl(1) + dFWZfl_aero(i) + dFWZfl_x(i) + dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                FWZ_fr(i+1) = FWZ_fr(1) + dFWZfr_aero(i) + dFWZfr_x(i) + dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                FWZ_rl(i+1) = FWZ_rl(1) + dFWZrl_aero(i) + dFWZrl_x(i) + dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                FWZ_rr(i+1) = FWZ_rr(1) + dFWZrr_aero(i) + dFWZrr_x(i) + dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)

                % Limiting the wheel loads to (almost) zero (Begrenzen der Radlasten auf (quasi) Null)
                if FWZ_fl(i+1) < 0
                    FWZ_fl(i+1) = 0.001;
                end
                if FWZ_fr(i+1) < 0
                    FWZ_fr(i+1) = 0.001;
                end
                if FWZ_rl(i+1) < 0
                    FWZ_rl(i+1) = 0.001;
                end
                if FWZ_rr(i+1) < 0
                    FWZ_rr(i+1) = 0.001;
                end

                % Axle loads - for dynamic radii (Achslasten)
                [FWZr(i+1), FWZf(i+1), FWZtot(i+1)] = calculateAxleloads(FWZ_rl(i+1), FWZ_rr(i+1), FWZ_fl(i+1), FWZ_fr(i+1));

                % Vertical tire stiffnesses - for dynamic radii (Vertikale Reifensteifigkeiten)  
                [cZ_fl(i+1), cZ_fr(i+1), cZ_rl(i+1), cZ_rr(i+1)] = calculateVtirestiff(Fz, cZ_tire, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1));

                % Dynamic tire radii (Dynamische Reifenradien)
                [Rdyn_fl(i+1), Rdyn_fr(i+1), Rdyn_rl(i+1), Rdyn_rr(i+1)] = calculateDynRadii(R0, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1), cZ_fl(i+1), cZ_fr(i+1), cZ_rr(i+1), cZ_rl(i+1));

                % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)
                [FWXmax_fl(i+1), FWXmax_fr(i+1), FWXmax_rl(i+1), FWXmax_rr(i+1), FWXmax_f(i+1), FWXmax_r(i+1)] = calculateLongiTireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

                % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)
                [FWYmax_fl(i+1), FWYmax_fr(i+1), FWYmax_rl(i+1), FWYmax_rr(i+1), FWYmax_f(i+1), FWYmax_r(i+1)] = calculateLatTireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

            end
            
            vWoBrake = vV;                                      % Save velocity without braking for log file.

            writeToLogfile('Simulated without brakes!', Debug, textAreaHandle);

            %% BRAKING POINT CALCULATION (BREMSPUNKTBERECHNUNG)
            [BrakeIndexes, NonBrakeApexes, vRev] = calculateBrakepoints(FB, Track, ApexIndexes, vAPEXmax, m_tot, downforce_multiplier, c_l, c_w, A_S, rho_L, ConstantDownforce, c_l_DRS, DRS_status, aero_ph, aero_pv, vV, k_R, FG, h_COG, wheelbase, track, lr, lf, GAMMA, TIRparam, FWZ_fl_stat, FWZ_fr_stat, FWZ_rl_stat, FWZ_rr_stat, R, s, brakeBias_setup, brakeFunction);

            %% Start values for simulation WITH BRAKES
            vV(1) = StartingSpeed + 0.005;
            
            t_x = 0;        % Reset Shift time
            
            gear = 1;
            
            t(1) = 0;       % [s] Time (Zeit)

            E_Accu(1) = 0;       % [J] Energy consumed by battery (Verbrauchte Energie Akku)
            E_heat(1) = 0; 
            E_Accu_Recu(i) = 0;  % [J] Energy recuperated by battery (Rekuperierte Energie Akku)

            % Supporting variables (Hilfsgrößen)
            z = 1;          % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

            %% SIMULATION WITH BRAKES (SIMULATION MIT BREMSEN)
            for i = 1:length(Track)-1

                % Checking at which apex the vehicle is (Überprüfen, vor welcher Apex das Auto ist)
                if ismember(i,ApexIndexes)  
                    z = z + 1;
                end
                
                % Saving the current gear for the result file
                gearSelection(i) = gear;
                
                if i > 1    % if i > 1 use real rpm instead of idle rpm
                    [ni(i), gear, t_x] = calculateGearbox(gearbox, idleRPM, n_shift, n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, n_Mmax, t_x, ni(i-1), t(i), t(i-1));      % Calculates Gearbox data and rpm
                else
                    [ni(i), gear, t_x] = calculateGearbox(gearbox, idleRPM, n_shift, n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, n_Mmax, t_x);      % Calculates Gearbox data and rpm
                end 

                [Faero(i)] = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i), ConstantDownforce, c_l_DRS, DRS_status(i)); % [N] Aerodynamic force
                
               %% Braking
                % Checking if braking is required (Prüfen, ob gebremst werden muss)
                if ismember(i,BrakeIndexes) && not(ismember(z,NonBrakeApexes))                      % Initiaion of braking process (Einleiten des Bremsvorgangs)     
                    % Braking 
                    Mi(i) = 0;                                    % [Nm] Motor torque (Motormoment)
                    BPPsignal(i) = 1;                             % [-] Brake signal (Bremssignal)
                    P_Bh(i) = FWZr(i)/FWZtot(i)*FB(i)*vV(i);      % [W] Rear braking power for recuperation (Bremsleistung hinten für Rekuperation)

                    P_M(i) = 0;
                    Mi(i) = 0;
                    motor_eff(i) = 0;
                    P_Mloss(i) = 0;

                    M_tractive(i) = 0;
                    P_tractive(i) = 0;
                    P_el(i) = 0;

                    FVX_fl(i) = 0;                   % [N] Tractive Force on front left wheel (AWD) 
                    FVX_fr(i) = 0;                   % [N] Tractive Force on front right wheel (AWD) 
                    FVX_f(i) = 0;              % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                    FVX_rl(i) = 0;                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                    FVX_rr(i) = 0;                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                    FVX(i) = 0;                % [N] Traction on rear axle (Zugkraft an der Hinterachse)
                     
                    
                    [FB_fl(i), FB_fr(i), FB_rl(i), FB_rr(i), ~, BrakeBias(i), ABS(i)] = calculateDeceleration(FB(i), m_tot, Fdr(i), FWXmax_fl(i), FWXmax_fr(i), FWXmax_rl(i), FWXmax_rr(i), brakeBias_setup);
                else
                %% Accelerating 
                
                    Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motor torque (single motor!)
                    FB(i) = 0;                                    % [N] Braking force
                    P_Bh(i) = 0;                                  % [W] Rear braking power (Bremsleistung hinten)
                    FB_fl(i) = 0;    
                    FB_fr(i) = 0;    
                    FB_rl(i) = 0;    
                    FB_rr(i) = 0;    

                    % Motor power & limitation to 80 kW from FS-Rules (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
                    rpmpointer = round(ni(i));                   % Pointer for efficiency table (Pointer für effizienztabelle)

                    if Mi(i) <= 0
                        torquepointer = 1;
                    else
                        torquepointer = round(Mi(i));               % Pointer for efficiency table (Pointer für effizienztabelle)
                    end 

                    P_M(i) = num_motors * Mi(i) * ni(i) / 60 * 2 * pi; % [W] Total motor power (P_el!)

                    % Limiting the maximal power when using an electric
                    % drivetrain
                    if ptype && P_M(i) > max_power
                        P_M(i) = max_power;                           % [W] Limited power (P_el!)
                        %Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Total Motor Torque!)
                    end

                    % Checks if the rpmpointer is higher than the maximum rpm and
                    % adjusts it if needed
                    if(rpmpointer > n_Mmax)
                        rpmpointer = n_Mmax;
                    elseif(rpmpointer < 1)
                        rpmpointer = 1;
                    end

                    % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
                    if ptype
                        motor_eff(i) = M_eff_inter(rpmpointer,torquepointer);
                    else
                        motor_eff(i) = 1;
                    end

                    P_Mloss(i) = P_M(i)*(1-(motor_eff(i)*drivetrain_eff*eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)

                    %P_el(i) = P_M(i);

                    % Calculation of motor power after deduction of efficiency of the inverter (ALL MOTORS!)
                    P_M(i) = P_M(i) - P_Mloss(i);  

                    % Calculation Overall Torque with real power
                    Mi(i) = P_M(i)*(60/ni(i)/2/pi);    

                    % Calculate the tractive forces on the wheels
                    [FVX_fl(i), FVX_fr(i), FVX_rl(i), FVX_rr(i), FVX(i), FVX_f(i), TC_front(i), TC(i)] = calculateTractiveForces(Mi(i), num_motors, i_G, gr, Rdyn_fl(i), Rdyn_fr(i), Rdyn_rl(i), Rdyn_rr(i), t_x, gear, FWXmax_f(i), FWXmax_r(i), t_shift);

                    M_tractive(i) = (FVX(i)+FVX_f(i))/(i_G*gr(gear)/Rdyn_rr(i));            % [Nm] Torque including tractive force
                    P_tractive(i) = M_tractive(i)/(60/ni(i)/2/pi);      % [kW] Motor power required for traction 
                    P_el(i) = (P_tractive(i)/(drivetrain_eff * motor_eff(i) * eta_inv));     % [kW] Motor power including efficiencies
                end

                % Driving resistances (Fahrwiderstände) & Vehicle (Fahrzeug)
                [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = calculateVehicleResistancesForces(k_R, FWZtot(i), rho_L, vV(i), c_w, A_S, m_tot, R(i), FVX(i), FVX_f(i), c_d_DRS, DRS_status(i), rpmpointer, n_Mmax, FB_fl(i), FB_fr(i), FB_rl(i), FB_rr(i));

                % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i)));     

                % ToDo Check vehicle speed before applying brakes 
                % Limit Braking before Apex if car is allready slower than
                % needed
                if ismember(i,BrakeIndexes) && vV(i+1) < vAPEXmax(z) && not(ismember(z,NonBrakeApexes))  % Begrenzen der Geschwindigkeit auf ApexGeschwindigkeit (Bremst solange bis Geschwindigkeiten gleich)
                     vV(i+1) = vAPEXmax(z);                         % [m/s] Total vehicle speed 
                end

                t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);                % [s] Time (Zeit)

                % Battery energy capacity (Energiemenge Akku)
                if (t(i+1)-t(i)) * P_Bh(i) > 16.5*10^3 * (t(i+1)-t(i)) 
                    E_Accu_Recu(i+1) = E_Accu_Recu(i) + 16.5*10^3 * (t(i+1)-t(i)); % [J] 
                else
                    E_Accu_Recu(i+1) = E_Accu_Recu(i) + (t(i+1)-t(i)) * P_Bh(i); % [J]
                end
                
                % Latschlängen für Längsschlupfberechnung nach Carter
                FU_fl(i) = FVX(i)/2-k_R*FWZ_fl(i+1)-FB_fl(i);   % [N] Umfangskräfte an einem Hinterrad
                FU_fr(i) = FVX(i)/2-k_R*FWZ_fr(i+1)-FB_fr(i);   % [N] Umfangskräfte an einem Hinterrad
                FU_rl(i) = FVX(i)/2-k_R*FWZ_rl(i+1)-FB_rl(i);   % [N] Umfangskräfte an einem Hinterrad
                FU_rr(i) = FVX(i)/2-k_R*FWZ_rr(i+1)-FB_rr(i);   % [N] Umfangskräfte an einem Hinterrad
                
                l_contact_patch_fl(i) = FWZ_fl(i)/(p_infl*bW);   % [m] Latschlänge vorne links
                l_contact_patch_fr(i) = FWZ_fr(i)/(p_infl*bW);   % [m] Latschlänge vorne rechts
                l_contact_patch_rl(i) = FWZ_rl(i)/(p_infl*bW);   % [m] Latschlänge hinten links
                l_contact_patch_rr(i) = FWZ_rr(i)/(p_infl*bW);   % [m] Latschlänge hinten rechts

                % Längsschlupf nach Carter
                my0 = 2;    % [-] Haftreibungsbeiwert   
                kappa_fl(i) = l_contact_patch_fl(i)/(2*Rdyn_fl(i))*my0*(wheelbase-sqrt(1-FU_fl(i)/(my0*FWZ_fl(i))));
                kappa_fr(i) = l_contact_patch_fr(i)/(2*Rdyn_fr(i))*my0*(wheelbase-sqrt(1-FU_fr(i)/(my0*FWZ_fr(i))));
                kappa_rl(i) = l_contact_patch_rl(i)/(2*Rdyn_rl(i))*my0*(wheelbase-sqrt(1-FU_rl(i)/(my0*FWZ_rl(i))));
                kappa_rr(i) = l_contact_patch_rr(i)/(2*Rdyn_rr(i))*my0*(wheelbase-sqrt(1-FU_rr(i)/(my0*FWZ_rr(i))));

                delta(i) = f*atan((wheelbase/1000)/sqrt(R(i)^2-(lr/1000)^2));       % [rad] Lenkwinkel
                beta(i) = f*atan((lr/1000)/sqrt(R(i)^2-(lr/1000)^2));               % [rad] Schwimmwinkel
                psi1(i) = vV(i)/R(i);                                               % [rad/s] Gierrate
                alpha_f(i) = 180/pi*(delta(i)-(lf/1000)/vV(i)*psi1(i)-beta(i));     % [°] Schräglaufwinkel vorne
                alpha_r(i) = 180/pi*((lr/1000)/vV(i)*psi1(i)-beta(i));              % [°] Schräglaufwinkel hinten


                E_Accu(i+1) = E_Accu(i) + (t(i+1)-t(i)) * P_el(i); % [J] 
                E_heat(i+1) = E_heat(i) + (P_el(i)-P_tractive(i)) * (t(i+1)-t(i));
                %E_Waerme(i+1) = P_el(i)*(1-M_eff_inter(i))+(E_Akku(i) + (t(i+1)-t(i)))*(1-eta_inv);    % [J] Motor losses + 5% for inverter losses; Drivetrain losses with heat (Motorverluste + 5% flat für inverterverlsute % Drivetrain Verluste auch Wärme)


                % Lateral forces on front and rear axle (Querkräfte an Vorder- und Hinterachse)
                FWYf(i) = lr/wheelbase*FVY(i);   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                FWYr(i) = lf/wheelbase*FVY(i);   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

                if FWYf(i) > FWYmax_f(i)
                   slipY_f(i) = 1;  
                end

                 if FWYr(i) > FWYmax_r(i)
                   slipY_r(i) = 1;  
                end

                 % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokräften)        
                 [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = calculateAeroforceOnWheels(Faero(i), aero_ph, aero_pv);

                 % Dynamic wheel load displacements in longitudinal direction (Dynamische Radlastverlagerungen in Längsrichtung)       
                 [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = calculateWheelloadLongDisp(h_COG, m_tot, aVX(i), wheelbase);

                 % Dynamic wheel load displacements in lateral direction (Dynamische Radlastverlagerung in Querrichtung)  
                 [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = calculateWheelloadLatDisp(h_COG, track, lr, lf, wheelbase, FVY(i));

                 % Wheel Loads (Radlasten)
                 FWZ_fl(i+1) = FWZ_fl(1) + dFWZfl_aero(i) + dFWZfl_x(i) + dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                 FWZ_fr(i+1) = FWZ_fr(1) + dFWZfr_aero(i) + dFWZfr_x(i) + dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                 FWZ_rl(i+1) = FWZ_rl(1) + dFWZrl_aero(i) + dFWZrl_x(i) + dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                 FWZ_rr(i+1) = FWZ_rr(1) + dFWZrr_aero(i) + dFWZrr_x(i) + dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)

                 % Limiting the wheel loads to (almost) zero (Begrenzen der Radlasten auf (quasi) Null)
                 if FWZ_fl(i+1) < 0
                     FWZ_fl(i+1) = 0.001;
                 end
                 if FWZ_fr(i+1) < 0
                     FWZ_fr(i+1) = 0.001;
                 end
                 if FWZ_rl(i+1) < 0
                     FWZ_rl(i+1) = 0.001;
                 end
                 if FWZ_rr(i+1) < 0
                     FWZ_rr(i+1) = 0.001;
                 end

                % Axle loads - for dynamic radii (Achslasten)      
                [FWZr(i+1), FWZf(i+1), FWZtot(i+1)] = calculateAxleloads(FWZ_rl(i+1), FWZ_rr(i+1), FWZ_fl(i+1), FWZ_fr(i+1)); 

                % Vertical tire stiffness - for dynamic radii (Vertikale Reifensteifigkeiten)  
                [cZ_fl(i+1), cZ_fr(i+1), cZ_rl(i+1), cZ_rr(i+1)] = calculateVtirestiff(Fz, cZ_tire, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1));

                % Dynamic tire radii (Dynamische Reifenradien)   
                [Rdyn_fl(i+1), Rdyn_fr(i+1), Rdyn_rl(i+1), Rdyn_rr(i+1)] = calculateDynRadii(R0, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1), cZ_fl(i+1), cZ_fr(i+1), cZ_rr(i+1), cZ_rl(i+1));

                % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)       
                [FWXmax_fl(i+1), FWXmax_fr(i+1), FWXmax_rl(i+1), FWXmax_rr(i+1),FWXmax_f(i+1), FWXmax_r(i+1)] = calculateLongiTireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

                % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)      
                [FWYmax_fl(i+1), FWYmax_fr(i+1), FWYmax_rl(i+1), FWYmax_rr(i+1), FWYmax_f(i+1), FWYmax_r(i+1)] = calculateLatTireforces(FWZ_fl(i), FWZ_fr(i),FWZ_rl(i), FWZ_rr(i), GAMMA, TIRparam, alpha_f(i), alpha_r(i));
                
                % Maximum cornering velocity 
                vVYmax(i+1) = sqrt((FWYmax_f(i+1)+FWYmax_r(i+1))*R(i+1)/m_tot); % Calculating maximum possible lateral velocity with given Tire forces [m/s] (inaccuaracy because tire force is based on aero force)

%                 %Akkuströme
%                 V_i(i) = sum(Voltage_Cellpack(:,i));
% 
%                 % Battery Currents (Akkuströme)
%                 A_accu_cell(i) = P_el(i) / V_i(i) / ncells_parallel;  
% 
%                 Current_Cellpack_Pointer(i) = P_M(i) / V_i(i) *10 ; %Strombelastung eines %er Parrallel Paketes in 0,1A parameter für die berechnung der korrigierten belastung mit höhren verlusten durch höhere zellströme
%                 if Current_Cellpack_Pointer(i) <= 1
%                     Current_Cellpack_Pointer(i)=1;
%                 end
% 
%                 if Current_Cellpack_Pointer(i) >= 1500 %begrenzen des max Zellstromes auf 30A pro Zelle im 5er parralelverbund also 150A
%                     Current_Cellpack_Pointer(i)=1500;
%                 end
% 
%                 VirtualCurrent_Cellpack(i) = CorrectedDischargeInterpolated(1,round(Current_Cellpack_Pointer(i))); %Berechnung der Virtuell höheren zellströme basierend auf den höheren verlsuten durch höhere Ströme
% 
%                 Energy_Cellpack(i) = (VirtualCurrent_Cellpack(i)*(t(i+1)-t(i))) - ((P_Bh(i)/V_i(i))*(t(i+1)-t(i))) ; %Energieverbrauch in As für ein 5erpacket an akkuzellen -> Akkustrom zum zeitpunkt i mal Zeitdifferenz zwischen i und i+1
%                 Energy_Cellpack_Total(i+1) = Energy_Cellpack_Total(i) + Energy_Cellpack(i); % Über Endurance Run Integrierte Energieverbrauch in As für ein 5erpacket an akkuzellen
% 
%                 Capacity_Cellpack(1:131,i+1) =  Capacity_Cellpack(1:131,i)- Energy_Cellpack(i); 
% 
%                 SOC_Cellpack(1:131,i+1) = Capacity_Cellpack(1:131,i)./Capacity_Cellpack(1:131,1); %Berechnung des SOC für den nächsten tick basierend auf der aktuellen cellcapacity und der im nächsten tick
% 
%                 SOC_Pointer(1:131,i+1) = round(SOC_Cellpack(1:131,i+1)*1000);
%                 Current_Cellpack_Pointer_Voltage(1,i+1) = round(Current_Cellpack_Pointer(i)/5);
% 
%                 if Current_Cellpack_Pointer_Voltage(i) <= 3
%                     Current_Cellpack_Pointer_Voltage(i)=3;
%                 end
%                 
%                 if size(Track,1) < SOC_Pointer(1:131,i+1)
%                     SOC_Pointer(1:131,i+1) = size(Track,1); 
%                 end
%                 
%                 if size(Track,1) < Current_Cellpack_Pointer_Voltage(1,i)
%                     Current_Cellpack_Pointer_Voltage(1,i) = size(Track,1);
%                 end
%                 
%                 %% @Lukas bitte prüfen if, wegen sonst auftretender Fehler bei langen Strecken!
%                 if (SOC_Pointer(1:131,i+1)>1001) 
%                     Voltage_Cellpack(1:131,i+1) = Voltage_inter(Current_Cellpack_Pointer_Voltage(1,i),1001);    
%                 else
%                     Voltage_Cellpack(1:131,i+1) = Voltage_inter(Current_Cellpack_Pointer_Voltage(1,i),SOC_Pointer(1:131,i+1));
%                 end
                
             
            end
            
            writeToLogfile('Simulated with brakes!', Debug, textAreaHandle);
            
            % Conversion of battery energy capacity (Umrechnen der Energiemengen des Akkus)
            E_Accu = E_Accu/(3.6e6);                % [kWh] Energy consumed by battery per lap (Verbrauchte Akku-Energie je Runde)
            E_heat = E_heat/(3.6e6);                % [kWh] 3.6e6 Joule conversion (3.6e6 Umrechnung Joule)
            E_Accu_Recu = E_Accu_Recu/(3.6e6);      % [kWh] Energy recuperated by batter per lap (Rekuperierte Akku-Energie je Runde)
            E_res = E_Accu - E_Accu_Recu;           % [kWh] Resulting energy consumption per lap (Resultierender Verbrauch je Runde)

            %%  Output of the values (Ausgabe der Werte)
            tEnd = toc;
            
            if Debug
                if (disciplineID == 2)
                    textAreaHandle.Value{end+1} = ['Endurance Time for ONE Lap: ' num2str(t(end)) ' s = ' num2str(t(end)/60) ' min'];
                    t_tot = t(end) *(22000/s(end));       % t(end) = Time for one lap ; 22000 [m] = length of Endurance ; s(end) length of the track.
                    textAreaHandle.Value{end+1} =  ['Endurance Total Time: ' num2str(t_tot(end)) ' s = ' num2str(t_tot(end)/60) ' min'];
                    textAreaHandle.Value{end+1} =  ['Total Travel Distance: ' num2str(s(end)) ' m'];
                    E_Accu_total_without_recu = E_Accu(end) * (22000/s(end));
                    textAreaHandle.Value{end+1} =  ['Total energy consumption W/O recuperation: ' num2str(E_Accu_total_without_recu) ' kWh'];
                    E_Accu_total = E_res(end) * (22000/s(end));
                    textAreaHandle.Value{end+1} =  ['Total energy consumption W recuperation: ' num2str(E_Accu_total) ' kWh'];
                    E_AccuCellPackage_energy_totalEndurance = (Energy_Cellpack_Total(end) * (22000/s(end)))/60^2;
                    textAreaHandle.Value{end+1} =  ['Total energy consumption in Ah: ' num2str(E_AccuCellPackage_energy_totalEndurance) ' Ah'];
                    textAreaHandle.Value{end+1} =  ['Average heat output: ' num2str(E_heat(end)/(t(end)/60^2)) ' kW'];
                else
                    %%  Output of the values (Ausgabe der Werte)
                    textAreaHandle.Value{end+1} = ['Endurance Time for ONE Lap: ' num2str(t(end)) ' s = ' num2str(t(end)/60) ' min'];
                    textAreaHandle.Value{end+1} =  ['Laptime: ' num2str(t(end)) ' s = ' num2str(t(end)/60) ' min'];
                    textAreaHandle.Value{end+1} =  ['Total Travel Distance: ' num2str(s(end)) ' m'];
                    textAreaHandle.Value{end+1} =  ['Total energy consumption W/O recuperation: ' num2str(E_Accu(end)) ' kWh'];
                    textAreaHandle.Value{end+1} =  ['Total energy consumption W recuperation: ' num2str(E_res(end)) ' kWh'];
                    E_AccuCellPackage_energy_totalEndurance = (Energy_Cellpack_Total(end))/60^2;
                    textAreaHandle.Value{end+1} =  ['Total energy consumption in Ah: ' num2str(E_AccuCellPackage_energy_totalEndurance) ' Ah'];
                    textAreaHandle.Value{end+1} =  ['Average heat output: ' num2str(E_heat(end)/(t(end)/60^2)) ' kW'];   
                end
            end

            t_tot = t(end);

            %% Writing the results to Mat File (Schreiben der Ergebnisse in Mat File)
            result.t_skidpad(steps) = t_skidpad;
            result.vV_skidpad(steps) = vV_skidpad;
            result.tEnd(steps) = tEnd;
            result.t_ges(steps) = t_tot;    
            result.Track = Track;
            result.aRev(:,steps) = aRev(:);
            result.aVX(:,steps) = aVX(:);
            result.aVY(:,steps) = aVY(:);
            result.BPPsignal(:,steps) = BPPsignal(:);
            result.cZ_vl(:,steps) = cZ_fl(:);
            result.cZ_vr(:,steps) = cZ_fr(:);
            result.cZ_hl(:,steps) = cZ_rl(:);
            result.cZ_hr(:,steps) = cZ_rr(:);
            result.dFWZhl_aero(:,steps) = dFWZrl_aero(:);
            result.dFWZhr_aero(:,steps) = dFWZrr_aero(:);
            result.dFWZvl_aero(:,steps) = dFWZfl_aero(:);
            result.dFWZvr_aero(:,steps) = dFWZfr_aero(:);
            result.dFWZhl_x(:,steps) = dFWZrl_x(:);
            result.dFWZhr_x(:,steps) = dFWZrr_x(:);
            result.dFWZvl_x(:,steps) = dFWZfl_x(:);
            result.dFWZvr_x(:,steps) = dFWZfr_x(:);
            result.dFWZhl_y(:,steps) = dFWZrl_y(:);
            result.dFWZhr_y(:,steps) = dFWZrr_y(:);
            result.dFWZvl_y(:,steps) = dFWZfl_y(:);
            result.dFWZvr_y(:,steps) = dFWZfr_y(:);
            result.E_Akku(:,steps) = E_Accu(:);
            result.E_Waerme(:,steps) = E_heat(:);
            result.E_Akku_Reku(:,steps) = E_Accu_Recu(:);
            result.E_ges(:,steps) = E_res(:);
            result.Faero(:,steps) = Faero(:);
            %result.FB = FB*ones(1,length(Track))(:);
            result.FVX(:,steps) = FVX(:);
            result.FVX_hl(:,steps) = FVX_rl(:);
            result.FVX_hr(:,steps) = FVX_rr(:);
            result.FVXre(:,steps) = FVXre(:);
            result.FVXid(:,steps) = FVXid(:);
            result.FR(:,steps) = FR(:);
            result.FL(:,steps) = FL(:);
            result.Fdr(:,steps) = Fdr(:);
            result.FVY(:,steps) = FVY(:);
            result.FWXmax_v(:,steps) = FWXmax_f(:);
            result.FWXmax_h(:,steps) = FWXmax_r(:);
            result.FWXmax_vl(:,steps) = FWXmax_fl(:);
            result.FWXmax_vr(:,steps) = FWXmax_fr(:);
            result.FWXmax_hl(:,steps) = FWXmax_rl(:);
            result.FWXmax_hr(:,steps) = FWXmax_rr(:);
            result.FWYv(:,steps) = FWYf(:);
            result.FWYh(:,steps) = FWYr(:);
            result.FWYmax_v(:,steps) = FWYmax_f(:);
            result.FWYmax_h(:,steps) = FWYmax_r(:);
            result.FWYmax_vl(:,steps) = FWYmax_fl(:);
            result.FWYmax_vr(:,steps) = FWYmax_fr(:);
            result.FWYmax_hl(:,steps) = FWYmax_rl(:);
            result.FWYmax_hr(:,steps) = FWYmax_rr(:);
            result.FWZges(:,steps) = FWZtot(:);
            result.FWZ_hl(:,steps) = FWZ_rl(:);
            result.FWZ_hr(:,steps) = FWZ_rr(:);
            result.FWZ_vl(:,steps) = FWZ_fl(:);
            result.FWZ_vr(:,steps) = FWZ_fr(:);
            result.FWZh(:,steps) = FWZr(:);
            result.FWZv(:,steps) = FWZf(:);
            result.Mi(:,steps) = Mi(:);
            result.ni(:,steps) = ni(:);
            result.P_M(:,steps) = P_M(:);
            result.P_el(:,steps) = P_el(:);
            result.Rdyn_vl(:,steps) = Rdyn_fl(:);
            result.Rdyn_vr(:,steps) = Rdyn_fr(:);
            result.Rdyn_hl(:,steps) = Rdyn_rl(:);
            result.Rdyn_hr(:,steps) = Rdyn_rr(:);
            result.RutschenY_v(:,steps) = slipY_f(:);
            result.RutschenY_h(:,steps) = slipY_r(:);
            result.t(:,steps) = t(:);
            result.TC(:,steps) = TC(:);
            result.Tirelimit(:,steps) = Tirelimit(:);
            result.vAPEXmax(:,steps) = vAPEXmax(:);
            result.vWoBrake(:,steps) = vWoBrake(:);
            result.vV(:,steps) = vV(:);
            result.vRev(:,steps) = vRev(:);
            result.vVYmax(:,steps) = vVYmax(:);
            result.A_Akkuzelle(:,steps) = A_accu_cell(:);
            result.ApexIndizes(:,steps) = ApexIndexes(:);
        %     result.M_eff_inter = M_eff_inter;
            result.motor_eff(:,steps) = motor_eff(:);
            result.gear(:,steps) = gearSelection(:);
                
            if logCellData
                result.Energy_Cellpack(:,steps) = Energy_Cellpack(:);
                result.VirtualCurrent_Cellpack(:,steps) = VirtualCurrent_Cellpack(:);
                result.V_i(:,steps) = V_i(:);
                result.Capacity_Cellpack(:,steps) = Capacity_Cellpack(:);
            end

            % Aero
            result.aero_ph(:,steps) = aero_ph(:);
            result.aero_pv(:,steps) = aero_pv(:);

            % Car Parameters needed to draw result plots and for setup
            % viewer.
            result.P_max(:,steps) = max_power(:);
            result.GAMMA(:,steps) = GAMMA(:); 
            result.m_tot(:,steps) = m_tot(:);
            result.A(:,steps) = A_S(:);
            result.FB(:,steps) = FB(:);
            result.LMUX(:,steps) = TIRparam.LMUX(:);
            result.LMUY(:,steps) = TIRparam.LMUY(:);
            result.c_l(:,steps) = c_l(:);
            result.c_w(:,steps) = c_w(:);
            result.c_l_DRS(:,steps) = c_l_DRS(:);
            result.c_d_DRS(:,steps) = c_d_DRS(:);
            result.downforce_multiplier(:,steps) = downforce_multiplier(:);
            result.drivetrain_eff(:,steps) = drivetrain_eff(:);
            result.num_motors(:,steps) = num_motors(:);
            result.t_L(:,steps) = t_L(:);
            result.rho_L(:,steps) = rho_L(:);
            result.i_G(:,steps) = i_G(:);
            result.z_pinion(:,steps) = z_pinion(:);
            result.z_chaindrive(:,steps) = z_chaindrive(:);
            result.track(:,steps) = track(:);
            result.wheelbase(:,steps) = wheelbase(:);
            result.n_Mmax(:,steps) = n_Mmax(:);
            result.thetaV_X(:,steps) = thetaV_X(:);
            result.thetaV_Y(:,steps) = thetaV_Y(:);
            result.thetaV_Z(:,steps) = thetaV_Z(:);
            result.trq_multiplier(:,steps) = trq_multiplier(:);
            result.m_balast(:,steps) = m_balast(:);
            result.m_driver(:,steps) = m_driver(:);
            result.m_ph(:,steps) = m_ph(:);
            result.p_infl(:,steps) = p_infl(:);
            result.h_COG(:,steps) = h_COG(:);
            result.x_COG(:,steps) = x_COG(:);
            result.h_COG_balast(:,steps) = h_COG_balast(:);
            result.h_COG_driver(:,steps) = h_COG_driver(:);
            result.x_vA(:,steps) = x_vA(:);
            result.DRS_status(:,steps) = DRS_status(:);
            result.i_param(:,steps) = gr(:);
            result.n_shift(:,steps) = n_shift(:);
            result.n_downshift(:,steps) = n_downshift(:);
            result.l_contact_patch_fl(:,steps) = l_contact_patch_fl(:);
            result.l_contact_patch_fr(:,steps) = l_contact_patch_fr(:);
            result.l_contact_patch_rl(:,steps) = l_contact_patch_rl(:);
            result.l_contact_patch_rr(:,steps) = l_contact_patch_rr(:);
            result.kappa_rl(:,steps) = kappa_rl(:);
            result.kappa_rr(:,steps) = kappa_rr(:);
            result.kappa_fl(:,steps) = kappa_fl(:);
            result.kappa_fr(:,steps) = kappa_fr(:);
            result.delta(:,steps) = delta(:);
            result.beta(:,steps) = beta(:);
            result.psi1(:,steps) = psi1(:);
            result.alpha_f(:,steps) = alpha_f(:);
            result.alpha_r(:,steps) = alpha_r(:);
            result.ConstantDownforce(:,steps) = ConstantDownforce(:);
            result.DRS_Radius(:,steps) = DRS_Radius(:);
            result.disciplineID = disciplineID;
            result.numOfLaps = numOfLaps;
            result.lapLength = lapLength;
            result.ptype(:,steps) = ptype(:);

            % sets button progress (progressbar)
            if sensitivityID2 ~= 0 
                currentProg = min(round((size(wbar,2)-2)*(steps/numSteps^2)),size(wbar,2)-2); 
            else
                currentProg = min(round((size(wbar,2)-2)*(steps/numSteps)),size(wbar,2)-2); 
            end
            
            processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 1) = 0.25391; % (royalblue)
            processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 2) = 0.41016;
            processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 3) = 0.87891;
            drawnow; % updates button progress
        end
    end
    
    [~, name, ~] = fileparts(setupFile);
        
        % Saves the results to a .mat file in the same location as the
        % setup file. The name is generated by adding _result to the name.
        savefilename = name + "_result.mat";

        save(path+"/"+savefilename, '-struct','result');
        
        writeToLogfile('File succesfully written', Debug, textAreaHandle);
end