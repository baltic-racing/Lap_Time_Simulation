function Vehiclesim_Endurance_GUI_Version(setupFile, path, TrackFileName, disciplineID, sensitivityID, minValue, stepSize, numSteps, processDataButtonHandle, textAreaHandle, sensitivityID2, minValue2, stepSize2, logCellData, Debug, StartingSpeed, numOfLaps)
    
    
     %% Laptime simulation
%     clear all; clc; close all;
    
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
    
    writeToLogfile('loaded Track!', Debug, textAreaHandle);

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
    FB_1 = FB;

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
    % For Debugging
    % PercentDRS = sum(DRS_status)/i*100;
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
        eta_inv = setup.invertor_eff;           % [-] Inverter efficiency
    else
        eta_inv = 1;
    end
    
    n_shift = setup.n_shift;            % engine rpm when shifting gears [1/min]
    n_downshift = setup.n_downshift;
    t_shift = setup.t_shift;            % time needed to shift gears [s]
    gr = setup.i_param(:);            % Gear ratio for each gear
    num_gears = length(gr);             % highest gear
    gearbox = setup.gearbox;
    
    % Sets gearratio to a constant 1 if no gearbox is present
    if ~gearbox
        gr(1) = 1;
    end
    
    % Gearbox Data (Getriebedaten)
    z_chaindrive = setup.z_chaindrive;      % [-] Number of teeth on sprocket
    z_pinion = setup.z_sprocket;            % [-] Number of teeth on pinion 
    i_P = setup.i_P;
    i_G = z_chaindrive / z_pinion * i_P;          % [-] Gear ratio (Motor to wheel)
    
    % Battery Data (Akkudaten)
    
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
    
    aRev = zeros(1,length(Track)-1);
    aVX = zeros(1,length(Track)-1);
    aVY = zeros(1,length(Track)-1);
    BPPsignal = zeros(1,length(Track)-1);
    cZ_fl = zeros(1,length(Track));
    cZ_fr = zeros(1,length(Track));
    cZ_rl = zeros(1,length(Track));
    cZ_rr = zeros(1,length(Track));
    dFWZrl_aero = zeros(1,length(Track)-1);
    dFWZrr_aero = zeros(1,length(Track)-1);
    dFWZfl_aero = zeros(1,length(Track)-1);
    dFWZfr_aero = zeros(1,length(Track)-1);
    dFWZrl_x = zeros(1,length(Track)-1);
    dFWZrr_x = zeros(1,length(Track)-1);
    dFWZfl_x = zeros(1,length(Track)-1);
    dFWZfr_x = zeros(1,length(Track)-1);
    dFWZrl_y = zeros(1,length(Track)-1);
    dFWZrr_y = zeros(1,length(Track)-1);
    dFWZfl_y = zeros(1,length(Track)-1);
    dFWZfr_y = zeros(1,length(Track)-1);
    
    E_Accu = zeros(1,length(Track));
%     E_ = zeros(1,length(Track));
    E_Accu_Recu = zeros(1,length(Track));
    E_heat = zeros(1,length(Track));
    
    Faero = zeros(1,length(Track)-1);
    FB = FB*ones(1,length(Track));
    FVX = zeros(1,length(Track)-1);
    FVX_f = zeros(1,length(Track)-1);
    FVX_fr = zeros(1,length(Track)-1);
    FVX_fl = zeros(1,length(Track)-1);
    FVX_rl = zeros(1,length(Track)-1);
    FVX_rr = zeros(1,length(Track)-1);
    FVXre = zeros(1,length(Track)-1);
    FVXid = zeros(1,length(Track)-1);
    FR = zeros(1,length(Track)-1);
    FL = zeros(1,length(Track)-1);
    Fdr = zeros(1,length(Track)-1);
    FVY = zeros(1,length(Track)-1);
    FWXmax_f = zeros(1,length(Track));
    FWXmax_r = zeros(1,length(Track));
    FWXmax_fl = zeros(1,length(Track));
    FWXmax_fr = zeros(1,length(Track));
    FWXmax_rl = zeros(1,length(Track));
    FWXmax_rr = zeros(1,length(Track));
    FWYf = zeros(1,length(Track)-1);
    FWYr = zeros(1,length(Track)-1);
    FWYmax_f = zeros(1,length(Track));
    FWYmax_r = zeros(1,length(Track));
    FWYmax_fl = zeros(1,length(Track));
    FWYmax_fr = zeros(1,length(Track));
    FWYmax_rl = zeros(1,length(Track));
    FWYmax_rr = zeros(1,length(Track));
    FWZtot = zeros(1,length(Track));
    FWZ_rl = zeros(1,length(Track));
    FWZ_rr = zeros(1,length(Track));
    FWZ_fl = zeros(1,length(Track));
    FWZ_fr = zeros(1,length(Track));
    FWZr = zeros(1,length(Track)-1);
    FWZf = zeros(1,length(Track)-1);
    Mi = zeros(1,length(Track)-1);
    M_tractive = zeros(1,length(Track)-1);
    ni = zeros(1,length(Track)-1);
    P_M = zeros(1,length(Track)-1);
    P_el = zeros(1,length(Track)-1);
    P_tractive = zeros(1,length(Track)-1);
    Rdyn_fl = zeros(1,length(Track));
    Rdyn_fr = zeros(1,length(Track));
    Rdyn_rl = zeros(1,length(Track));
    Rdyn_rr = zeros(1,length(Track));
    slipY_f = zeros(1,length(Track));
    slipY_r = zeros(1,length(Track));
    t = zeros(1,length(Track)-1);
    
    l_contact_patch_fl = zeros(1,length(Track));
    l_contact_patch_fr = zeros(1,length(Track));
    l_contact_patch_rl = zeros(1,length(Track));
    l_contact_patch_rr = zeros(1,length(Track));
    
    kappa_rl = zeros(1,length(Track));
    kappa_rr = zeros(1,length(Track));
    kappa_fl = zeros(1,length(Track));
    kappa_fr = zeros(1,length(Track));
                
    delta = zeros(1,length(Track));
    beta = zeros(1,length(Track));
    psi1 = zeros(1,length(Track));
    alpha_f = zeros(1,length(Track));
    alpha_r = zeros(1,length(Track));
    
    TC = zeros(1,length(Track));
    TC_front = zeros(1,length(Track));
    ABS = zeros(1,length(Track));
    %DRS_status = zeros(1,length(Track)-1);
    gearSelection = zeros(1,length(Track)-1);
    
    Tirelimit = zeros(1,length(Track));
    vAPEXmax = zeros(1,length(ApexIndexes));
    vV = zeros(1,length(Track)-1);
    vRev = zeros(1,length(Track)-1);
    vVYmax = zeros(1,length(Track)-1);
    A_accu_cell = zeros(1,length(Track));
    P_Mloss = zeros(1, length(Track));
    P_Bh = zeros(1, length(Track));
%     M_eff_inter = zeros(l, length(Track));
    motor_eff = zeros(1, length(Track));
    Capacity_Cellpack = Parrallelcellgroups(1:131,1);
    SOC_Cellpack(1:131,1) = 1;
    Voltage_Cellpack(1:131,1) = 4.2;
    V_i = zeros(1,length(Track));
    VirtualCurrent_Cellpack = zeros(1,length(Track));
    Current_Cellpack_Pointer = zeros(1,length(Track));
    Energy_Cellpack = zeros(1,length(Track));
    Energy_Cellpack_Total = zeros(1,length(Track));
    SOC_Pointer = zeros(131,length(Track));
    Current_Cellpack_Pointer_Voltage = zeros(1,length(Track));

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

            % Axle and wheel loads (static) ((Statische) Achs- und Radlasten)
            FWZtot(1) = FG;                 % [N] Static total axle load (Statische Gesamtachslast)
            FWZr(1) = m_ph/100*FWZtot(1);   % [N] Static rear axle load (Statische Achslast hinten)
            FWZf(1) = FWZtot(1)-FWZr(1);    % [N] Static front axle load (Statische Achslast vorne)
            FWZ_fr(1) = FWZf(1)/2;          % [N] Static front right wheel load (Statische Radlast vorne rechts)  
            FWZ_fl(1) = FWZf(1)/2;          % [N] Static front left wheel load (Statische Radlast vorne links)
            FWZ_rr(1) = FWZr(1)/2;          % [N] Static rear right wheel load (Statische Radlast hinten rechts)
            FWZ_rl(1) = FWZr(1)/2;          % [N] Static rear left wheel load (Statische Radlast hinten links)

            % Interpolated vertical tire stiffness (Vertikale Reifensteifigkeiten interpoliert) 
            [cZ_fl(1), cZ_fr(1), cZ_rl(1), cZ_rr(1)] = vtirestiff(Fz, cZ_tire, FWZ_fl(1), FWZ_fr(1), FWZ_rl(1), FWZ_rr(1));

            % Dynamic tire radii (stationary = static) (Dynamische Reifenradien (im Stand = statisch))
            [Rdyn_fl(1), Rdyn_fr(1), Rdyn_rl(1), Rdyn_rr(1)] = dyn_radii(R0, FWZ_fl(1), FWZ_fr(1), FWZ_rl(1), FWZ_rr(1), cZ_fl(1), cZ_fr(1), cZ_rr(1), cZ_rl(1));

            FWXmax_r(1) = Inf;

            % Static tire loads for calculation (Statische Reifenlasten für Berechnung)
            FWZ_fl_stat = FWZ_fl(1); 
            FWZ_fr_stat = FWZ_fr(1);
            FWZ_rl_stat = FWZ_rl(1);
            FWZ_rr_stat = FWZ_rr(1);

            %% Calculation of the maximum apex speed for all apexes (numerically) (Berechnen der maximalen Kurvengeschwindigkeiten für alle Apexes (numerisch))
            for i = 1:length(ApexIndexes)

                FWYf(i) = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
                FWYr(i) = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
                FWYmax_f(i) = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal übertragbare Querkraft Vorderachse)
                FWYmax_r(i) = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal übertragbare Querkraft Hinterachse)
                vV(i) = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)

                while  FWYf(i) < FWYmax_f(i) && FWYr(i) < FWYmax_r(i) && vV(i) < 30

                    vV(i) = vV(i) + 0.01;   % [m/s] Increaing vehicle speed (Erhöhen der Fahrzeuggeschwindigkeit)
              
                    Faero(i) = aeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i), ConstantDownforce, c_l_DRS, DRS_status(i)); % [N] Aerodynamic force

                    FVY(i) = m_tot*vV(i)^2/R(ApexIndexes(i));    % [N] Centrifugal force (Zentrifugalkraft)

                    aVY(i) = vV(i)^2/R(ApexIndexes(i));  % [m/s²] Lateral acceleration (Querbeschleunigung)

                    % Lateral forces to be applied on front and rear axle (Aufzubringende Querkräfte an Vorder- und Hinterachse)
                    FWYf(i) = lr/wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                    FWYr(i) = lf/wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

                    % Wheel load transfer due to drag forces (Radlastverlagerung in Folge von Aerokräften) 
                    [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = aeroforce_onwheels(Faero(i), aero_ph, aero_pv);

                    % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in Längsrichtung = 0 angenommen)
                    [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = wheelload_longdisp(h_COG, 0, aVX(i), wheelbase); % Loads = 0 assumed

                    % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                    [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = wheelload_latdisp(h_COG, track, lr, lf, wheelbase, FVY(i));

                    % Wheel loads (Radlasten)
                    FWZ_fl(i) = FWZ_fl_stat + dFWZfl_aero(i) + dFWZfl_x(i) + dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                    FWZ_fr(i) = FWZ_fr_stat + dFWZfr_aero(i) + dFWZfr_x(i) + dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                    FWZ_rl(i) = FWZ_rl_stat + dFWZrl_aero(i) + dFWZrl_x(i) + dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                    FWZ_rr(i) = FWZ_rr_stat + dFWZrr_aero(i) + dFWZrr_x(i) + dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)   

                    % Maximum transmissible tire forces in longitudinal direction = 0 assumed (because longitudinal wheel loads = 0 assumed) 

                    % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)    
                    [FWYmax_f(i), FWYmax_r(i)] = lat_tireforces(FWZ_fl(i), FWZ_fr(i),FWZ_rl(i), FWZ_rr(i), GAMMA, TIRparam);

                end

                vAPEXmax(i) = vV(i);   % [m/s] Maximum speed for any apex (Maximalgeschwindigkeit für jede Apex)
            end

            writeToLogfile('caclculated Apex Speeds!', Debug, textAreaHandle);
            
            %% Start/Initial values for first simulation run WITHOUT BRAKES (Startwerte für ersten Simulationslauf OHNE BREMSEN)
            vV(1) = StartingSpeed + 0.001;
            
            gear = 1;

            t(1) = 0;      % [s] Time (Zeit)

            % Supporting variables (Hilfsgrößen)
            z = 1;         % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

            %% Simulation WITHOUT BRAKES (Simulation OHNE BREMSEN)
            NonBrakeApexes = [];

            for i = 1:length(Track)-1
                
                if gearbox
                    if i > 1 && ni(i-1) < n_downshift && gear > 1                          % Shift down
                        gear = gear - 1;
                    elseif i > 1 && ni(i-1) >= n_shift && gear < num_gears                 % Shift one gear up when shifting rpm is reached
                        gear = gear + 1;
                        t_x = t(i)-t(i-1);                                      % step size of the time variable       
                    end
                    ni(i) = vV(i) * 30 / pi * gr(gear) * i_G / Rdyn_rl(i);      % [1/min] calculate rpm with gearbox
                else
                    % Determination of motor speed (no gearbox)
                    ni(i) = vV(i) * 30 / pi * i_G / Rdyn_rl(i);                 % [1/min] Determine current motor speed 
                end
                
                if ni(i) >= n_Mmax % Drehzahlbegrenzer
                    ni(i) = n_Mmax;
                end     

                % Determination of aero forces and motor torque (Bestimmen von Aero-Kräften und Motormoment)      
                Faero(i) = aeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i), ConstantDownforce, c_l_DRS, DRS_status(i)); % [N] Aerodynamic force

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
                    Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Begrenzen des Moments)
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

                P_M(i) = P_M(i) * drivetrain_eff * motor_eff(i) * eta_inv;  % Calculation of motor power after deduction of efficiency of the inverter
                Mi(i) = P_M(i)*(60/ni(i)/2/pi);                       

                % Calculation of tractive forces (Berechnen der Zugkraft)
                if num_motors == 4
                    FVX_fl(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_fl(i);                   % [N] Tractive Force on front left wheel (AWD) 
                    FVX_fr(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_fr(i);                   % [N] Tractive Force on front right wheel (AWD) 
                    FVX_f(i) = FVX_fl(i) + FVX_fr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                    FVX_rl(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_rl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                    FVX_rr(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_rr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                    FVX(i) = FVX_rl(i) + FVX_rr(i);      % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                    if FVX_f(i) > FWXmax_f(i)     % Limiting the tractive force to the traction limit front axle
                        FVX_f(i) = FWXmax_f(i);
                        TC_front(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
                    end
                else
                    FVX_rl(i) = Mi(i)/2*i_G*gr(gear)/Rdyn_rl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                    FVX_rr(i) = Mi(i)/2*i_G*gr(gear)/Rdyn_rr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                    FVX(i) = FVX_rl(i) + FVX_rr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                    FVX_f(i) = 0;
                end

                if FVX(i) > FWXmax_r(i)     % Limiting the tractive force to the traction limit rear axle
                    FVX(i) = FWXmax_r(i);
                    TC(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
                end


                % Driving resistances (Fahrwiderstände) & Vehicle (Fahrzeug)        
               [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = vehicle_resistances_forces(k_R, FWZtot(i), rho_L, vV(i), c_w, A_S, m_tot, R(i), FVX(i), FVX_f(i), 0, c_d_DRS, DRS_status(i));

                if ismember(i,ApexIndexes)
                    if vV(i) > vAPEXmax(z)   % Limiting curve maximum speeds at apexes (Begrenzen auf maximale Kurvengeschwindigkeit in Apexes)
                        vV(i) = vAPEXmax(z);
                    end
                    z = z + 1;
                end

                vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Time (Zeit)

                % Lateral forces to be applied on front and rear axles (Aufzubringende Querkräfte an Vorder- und Hinterachse)
                FWYf(i) = lr/wheelbase*FVY(i);   % [N] Lateral force to be applied on front axle (Aufzubringende Querkraft der Vorderachse)
                FWYr(i) = lf/wheelbase*FVY(i);   % [N] Lateral force to be applied on rear axle (Aufzubringende Querkraft der Hinterachse)

                % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokräften)  
                [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = aeroforce_onwheels(Faero(i), aero_ph, aero_pv);

                % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerungen in Längsrichtung)
                [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = wheelload_longdisp(h_COG, m_tot, aVX(i), wheelbase);

                % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = wheelload_latdisp(h_COG, track, lr, lf, wheelbase, FVY(i));

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
                [FWZr(i+1), FWZf(i+1), FWZtot(i+1)] = axleloads(FWZ_rl(i+1), FWZ_rr(i+1), FWZ_fl(i+1), FWZ_fr(i+1));

                % Vertical tire stiffnesses - for dynamic radii (Vertikale Reifensteifigkeiten)  
                [cZ_fl(i+1), cZ_fr(i+1), cZ_rl(i+1), cZ_rr(i+1)] = vtirestiff(Fz, cZ_tire, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1));

                % Dynamic tire radii (Dynamische Reifenradien)
                [Rdyn_fl(i+1), Rdyn_fr(i+1), Rdyn_rl(i+1), Rdyn_rr(i+1)] = dyn_radii(R0, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1), cZ_fl(i+1), cZ_fr(i+1), cZ_rr(i+1), cZ_rl(i+1));

                % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)
                [FWXmax_fl(i+1), FWXmax_fr(i+1), FWXmax_rl(i+1), FWXmax_rr(i+1), FWXmax_f(i+1), FWXmax_r(i+1)] = longi_tireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam);

                % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)
                [FWYmax_fl(i+1), FWYmax_fr(i+1), FWYmax_rl(i+1), FWYmax_rr(i+1), FWYmax_f(i+1), FWYmax_r(i+1)] = lat_tireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam);

            end

            writeToLogfile('Simulated without brakes!', Debug, textAreaHandle);

            %% BRAKING POINT CALCULATION (BREMSPUNKTBERECHNUNG)
            BrakeIndexes = [];

            k = length(ApexIndexes);
            
            % Reset FB
            FB = FB_1*ones(1,length(Track));

            while k >= 1

                Counter = 0;
                j = ApexIndexes(k);

                while vRev(j-1) < vV(j-1)
                    
                    [~, FL(j-1), Fdr(j-1), FVY(j-1), ~, ~] = vehicle_resistances_forces(k_R, FG, rho_L, vRev(j-1), c_w, A_S, m_tot, R(j-1), FVX(j-1), FVX_f(j-1), FB(j-1), c_d_DRS, DRS_status(j-1));
                    
                    Faero(j-1) = aeroforce(downforce_multiplier, c_l, A_S, rho_L, vRev(j-1), ConstantDownforce, c_l_DRS, DRS_status(j-1));   % [N] Aerodynamic force
                    FR(j-1) = k_R*(FG+Faero(j-1));                                              % [N] Rolling Resistance
                    %FL(j-1) = rho_L*vRev(j-1)^2/2*c_w*A_S;                                      % [N] Air resistance
                    %Fdr(j-1) = FR(j-1)+FL(j-1);                                                 % [N] Overall Resistance 
                    
                    aRev(j-1) = (-Fdr(j-1)-FB(j-1))/m_tot;
                    vRev(j-2) = sqrt(vRev(j-1)^2-2*aRev(j-1)*(s(j-1)-s(j-2)));
                    j = j - 1;
                    Counter = Counter + 1;
                end

                if Counter > 0
                    BrakeIndexes = [BrakeIndexes j-1:ApexIndexes(k)-1];
                end

                if k > 1 && j < ApexIndexes(k-1)
                    NonBrakeApexes = [NonBrakeApexes k-1];
                    k = k - 1;
                end

                k = k - 1;

            end

            %% Start values for simulation WITH BRAKES
            vV(1) = StartingSpeed + 0.001;
            
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
                
                if gearbox
                    if i > 1 && ni(i-1) < n_downshift && gear > 1                          % Shift down
                        gear = gear - 1;  
                    elseif i > 1 && ni(i-1) >= n_shift && gear < num_gears                 % Shift one gear up when shifting rpm is reached
                        gear = gear + 1;
                        t_x = t(i)-t(i-1);                                      % step size of the time variable       
                    end
                    
                    ni(i) = vV(i) * 30 / pi * gr(gear) * i_G / Rdyn_rl(i);     % [1/min] calculate rpm with gearbox
                    n_wheel(i) = vV(i) * 60 / pi / Rdyn_rl(i);
                else
                    % Determination of motor speed (no gearbox)
                    ni(i) = vV(i) * 30 / pi * i_G / Rdyn_rl(i); % [1/min] Determine current motor speed 
                end
                
                if ni(i) >= n_Mmax % Drehzahlbegrenzer
                    ni(i) = n_Mmax;
                end     

               [Faero(i)] = aeroforce(downforce_multiplier, c_l, A_S, rho_L, vV(i), ConstantDownforce, c_l_DRS, DRS_status(i)); % [N] Aerodynamic force
                
               %% Braking
                % Checking if braking is required (Prüfen, ob gebremst werden muss)
                if ismember(i,BrakeIndexes)                       % Initiaion of braking process (Einleiten des Bremsvorgangs)     
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
                    FVX_f(i) = -FB(i)/2;              % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                    FVX_rl(i) = 0;                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                    FVX_rr(i) = 0;                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                    FVX(i) = -FB(i)/2;                % [N] Traction on rear axle (Zugkraft an der Hinterachse)
                    
                    if abs(FVX_f(i)) > FWXmax_f(i)     % Limiting the tractive force to the traction limit rear axle
                        FB(i) = -FWXmax_f(i);
                        ABS(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
                    end     
                    
                    if abs(FVX(i)) > FWXmax_r(i)     % Limiting the tractive force to the traction limit rear axle
                        FB(i) = -FWXmax_r(i);
                        ABS(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
                    end     
                else
                %% Accelerating 
                    Mi(i) = interp1(n,M,ni(i),'linear','extrap'); % [Nm] Motor torque (single motor!)
                    FB(i) = 0;                                    % [N] Braking force
                    P_Bh(i) = 0;                                  % [W] Rear braking power (Bremsleistung hinten)

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

                    % Calculation of tractive forces (Berechnen der Zugkraft)
                    if num_motors == 4
                        FVX_fl(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_fl(i);                   % [N] Tractive Force on front left wheel (AWD) 
                        FVX_fr(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_fr(i);                   % [N] Tractive Force on front right wheel (AWD) 
                        FVX_f(i) = FVX_fl(i) + FVX_fr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                        FVX_rl(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_rl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                        FVX_rr(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_rr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                        FVX(i) = FVX_rl(i) + FVX_rr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                        if FVX_f(i) > FWXmax_f(i)     % Limiting the tractive force to the traction limit front axle
                            FVX_f(i) = FWXmax_f(i);
                            TC_front(i) = 1;              % front Traction control "on" 
                        end
                    elseif num_motors == 2
                        FVX_rl(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_rl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                        FVX_rr(i) = Mi(i)/num_motors*i_G*gr(gear)/Rdyn_rr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                        FVX(i) = FVX_rl(i) + FVX_rr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                        FVX_f(i) = 0;
                    else
                        FVX_rl(i) = Mi(i)/2*i_G*gr(gear)/Rdyn_rl(i);                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                        FVX_rr(i) = Mi(i)/2*i_G*gr(gear)/Rdyn_rr(i);                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                        FVX(i) = FVX_rl(i) + FVX_rr(i);                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                        FVX_f(i) = 0;
                    end

                    if FVX(i) > FWXmax_r(i)     % Limiting the tractive force to the traction limit rear axle
                        FVX(i) = FWXmax_r(i);
                        TC(i) = 1;              % Traction control "on" (Traktionskontrolle "an")
                    end       

                    M_tractive(i) = (FVX(i)+FVX_f(i))/(i_G*gr(gear)/Rdyn_rr(i));            % [Nm] Torque including tractive force
                    P_tractive(i) = M_tractive(i)/(60/ni(i)/2/pi);      % [kW] Motor power required for traction 
                    P_el(i) = (P_tractive(i)/(drivetrain_eff * motor_eff(i) * eta_inv));     % [kW] Motor power including efficiencies
                end


                % Driving resistances (Fahrwiderstände) & Vehicle (Fahrzeug)
                [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = vehicle_resistances_forces(k_R, FWZtot(i), rho_L, vV(i), c_w, A_S, m_tot, R(i), FVX(i), FVX_f(i), FB(i), c_d_DRS, DRS_status(i));

                % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i)));     

                % ToDo Check vehicle speed before applying brakes (ToDo Vor Beginn des Bremsens Geschwindigkeit prüfen)
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
                FT(i) = FVX(i)/2-k_R*FWZ_rl(i+1)-FB(i)/4;   % [N] Umfangskräfte an einem Hinterrad
                l_contact_patch_fl(i) = FWZ_fl(i)/(p_infl*bW);   % [m] Latschlänge vorne links
                l_contact_patch_fr(i) = FWZ_fr(i)/(p_infl*bW);   % [m] Latschlänge vorne rechts
                l_contact_patch_rl(i) = FWZ_rl(i)/(p_infl*bW);   % [m] Latschlänge hinten links
                l_contact_patch_rr(i) = FWZ_rr(i)/(p_infl*bW);   % [m] Latschlänge hinten rechts

                % Längsschlupf nach Carter
                my0 = 2;    % [-] Haftreibungsbeiwert   
                kappa_rl(i) = l_contact_patch_rl(i)/(2*Rdyn_rl(i))*my0*(wheelbase-sqrt(1-FT(i)/(my0*FWZ_rl(i))));
                kappa_rr(i) = l_contact_patch_rr(i)/(2*Rdyn_rr(i))*my0*(wheelbase-sqrt(1-FT(i)/(my0*FWZ_rr(i))));
                kappa_fl(i) = l_contact_patch_fl(i)/(2*Rdyn_fl(i))*my0*(wheelbase-sqrt(1-FT(i)/(my0*FWZ_fl(i))));
                kappa_fr(i) = l_contact_patch_fr(i)/(2*Rdyn_fr(i))*my0*(wheelbase-sqrt(1-FT(i)/(my0*FWZ_fr(i))));
                
                % Lenkwinkel, Schwimmwinkel, Gierrate, Schräglaufwinkel
                if R(i) > 0
                    f = -1;     % [-] Hilfsgröße für Richtung der Winkel
                else
                    f = 1;
                end

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
                 [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = aeroforce_onwheels(Faero(i), aero_ph, aero_pv);

                 % Dynamic wheel load displacements in longitudinal direction (Dynamische Radlastverlagerungen in Längsrichtung)       
                 [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = wheelload_longdisp(h_COG, m_tot, aVX(i), wheelbase);

                 % Dynamic wheel load displacements in lateral direction (Dynamische Radlastverlagerung in Querrichtung)  
                 [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = wheelload_latdisp(h_COG, track, lr, lf, wheelbase, FVY(i));

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
                [FWZr(i+1), FWZf(i+1), FWZtot(i+1)] = axleloads(FWZ_rl(i+1), FWZ_rr(i+1), FWZ_fl(i+1), FWZ_fr(i+1)); 


                % Vertical tire stiffness - for dynamic radii (Vertikale Reifensteifigkeiten)  
                [cZ_fl(i+1), cZ_fr(i+1), cZ_rl(i+1), cZ_rr(i+1)] = vtirestiff(Fz, cZ_tire, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1));

                % Dynamic tire radii (Dynamische Reifenradien)   
                [Rdyn_fl(i+1), Rdyn_fr(i+1), Rdyn_rl(i+1), Rdyn_rr(i+1)] = dyn_radii(R0, FWZ_fl(i+1), FWZ_fr(i+1), FWZ_rl(i+1), FWZ_rr(i+1), cZ_fl(i+1), cZ_fr(i+1), cZ_rr(i+1), cZ_rl(i+1));

                % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)       
                [FWXmax_fl(i+1), FWXmax_fr(i+1), FWXmax_rl(i+1), FWXmax_rr(i+1),FWXmax_f(i+1), FWXmax_r(i+1)] = longi_tireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam);

                % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)      
                [FWYmax_fl(i+1), FWYmax_fr(i+1), FWYmax_rl(i+1), FWYmax_rr(i+1), FWYmax_f(i+1), FWYmax_r(i+1)] = lat_tireforces(FWZ_fl(i), FWZ_fr(i),FWZ_rl(i), FWZ_rr(i), GAMMA, TIRparam);
                
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

            % Car Parameters needed to draw result plots
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

%% Function calls

function [BT_fl, BT_fr, BT_rl, BT_rr] = calculateBrakeTorque(vV, n_wheel, BrakePressure, boreFront, boreRear)
    
    R_m = (R_o + R_i) / 2;

    friction_s = 0.35;
    friciton_k = 0.35;
    
    if n_wheel(i) > 0 
        BT_fl = (friciton_k * BrakePressure * pi * boreFront^2 * R_m * 2) / 4;
        BT_fr = (friciton_k * BrakePressure * pi * boreFront^2 * R_m * 2) / 4;
        BT_rl = (friciton_k * BrakePressure * pi * boreRear^2 * R_m * 2) / 4;
        BT_rr = (friciton_k * BrakePressure * pi * boreRear^2 * R_m * 2) / 4;
    else
        BT_fl = (friciton_s)/4
    end
end

function [dFWZrl_aero, dFWZrr_aero, dFWZfl_aero, dFWZfr_aero] = aeroforce_onwheels(Faero, aero_ph, aero_pv)
%% Individual aerodynamic forces acting on each wheel
%   For calculation of wheel load transfer

        dFWZrl_aero = Faero/2*aero_ph;   % [N] Aerodynamic force on rear left wheel (Aerokraft auf linkes Hinterrad)
        dFWZrr_aero = Faero/2*aero_ph;   % [N] Aerodynamic force on rear right wheel (Aerokraft auf rechtes Hinterrad)
        dFWZfl_aero = Faero/2*aero_pv;   % [N] Aerodynamic force on front left wheel (Aerokraft auf linkes Vorderrad)
        dFWZfr_aero = Faero/2*aero_pv;   % [N] Aerodynamic force on front right wheel (Aerokraft auf rechtes Vorderrad)
            
end

function [dFWZfl_x, dFWZrr_x, dFWZrl_x, dFWZfr_x] = wheelload_longdisp(h_COG, m_ges, aVX, l)
%% Dynamic wheel load displacement in longitudinal direction
%   For calculation of wheel load transfer

        dFWZfl_x = -m_ges*aVX*h_COG/l/2;  % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
        dFWZfr_x = -m_ges*aVX*h_COG/l/2;  % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
        dFWZrl_x = m_ges*aVX*h_COG/l/2;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
        dFWZrr_x = m_ges*aVX*h_COG/l/2;   % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
        
end

function [dFWZfl_y, dFWZfr_y, dFWZrl_y, dFWZrr_y] = wheelload_latdisp(h_COG, B, lh, lv, l, FVY) 
%% Dynamic wheel load displacement in lateral direction
%   For calculation of wheel load transfer

        dFWZfl_y = -h_COG/B*lh/l*FVY;   % [N] Dynamic wheel load transfer to front left wheel (Dynamische Radlastverlagerung linkes Vorderrad)
        dFWZfr_y = h_COG/B*lh/l*FVY;    % [N] Dynamic wheel load transfer to front right wheel (Dynamische Radlastverlagerung rechtes Vorderrad)
        dFWZrl_y = -h_COG/B*lv/l*FVY;   % [N] Dynamic wheel load transfer to rear left wheel (Dynamische Radlastverlagerung linkes Hinterrad)
        dFWZrr_y = h_COG/B*lv/l*FVY;    % [N] Dynamic wheel load transfer to rear right wheel (Dynamische Radlastverlagerung rechtes Hinterrad)
end

function Faero = aeroforce(downforce_multiplier, c_l, A_S, rho_L, vV, ConstantDownforce, c_l_DRS, DRS_status)
%% Aerodynamic force
        if (DRS_status == 1)
            Faero = downforce_multiplier * ((c_l_DRS * A_S * rho_L * (vV)^2 / 2) + ConstantDownforce);  
        else
            Faero = downforce_multiplier * ((c_l * A_S * rho_L * (vV)^2 / 2) + ConstantDownforce);  
        end
end

function [FR, FL, Fdr, FVY, aVX, aVY] = vehicle_resistances_forces(k_R, FWZges, rho_L, vV, c_d, A_S, m_tot, R, FVX, FVX_f, FB, c_d_DRS, DRS_status)
%% Driving resistances and vehicle

        FR = k_R*FWZges;                % [N] Rolling resistance (Rollwiderstand)
        
        if DRS_status == 1
            FL = rho_L*(vV^2/2)*c_d_DRS*A_S;
        else
            FL = rho_L*(vV^2/2)*c_d*A_S;      % [N] Air Resistance/Drag (Luftwiderstand)
        end

        Fdr = FR + FL;                  % [N] Total resistance (Gesamtwiderstand)
        FVY = m_tot*vV^2/R;             % [N] Centrifugal force (Zentrifugalkraft)
        
        aVX = ((FVX+FVX_f)-Fdr)/m_tot;       % [m/s²] Longitudinal acceleration (Längsbeschleunigung)
        aVY = vV^2/R;                   % [m/s²] Lateral acceleration (Querbeschleunigung)
end

function [FWZr, FWZf, FWZges] = axleloads(FWZ_rl, FWZ_rr, FWZ_fl, FWZ_fr)
%% Axle loads calculation
        FWZr = FWZ_rl + FWZ_rr;  % [N] Rear axle load (Achslast hinten)
        FWZf = FWZ_fl + FWZ_fr;  % [N] Front axle load (Achslast vorne)
        FWZges = FWZr + FWZf;    % [N] Total axle and wheel load (Gesamtachs- und Radlast)
        
end

function [Rdyn_fl, Rdyn_fr, Rdyn_rl, Rdyn_rr] = dyn_radii(R0, FWZ_vl, FWZ_vr, FWZ_hl, FWZ_hr, cZ_vl, cZ_vr, cZ_hr, cZ_hl)
%% Dynamic tire radii calculation

        Rdyn_fl = R0 - FWZ_vl/cZ_vl;    % [m] Dynamic front left tire radius (Dynamischer Reifenradius)
        Rdyn_fr = R0 - FWZ_vr/cZ_vr;    % [m] Dynamic front right tire radius (Dynamischer Reifenradius)
        Rdyn_rl = R0 - FWZ_hl/cZ_hl;    % [m] Dynamic rear left tire radius (Dynamischer Reifenradius)
        Rdyn_rr = R0 - FWZ_hr/cZ_hr;    % [m] Dynamic rear right tire radius (Dynamischer Reifenradius)        
end

function [cZ_fl, cZ_fr, cZ_rl, cZ_rr] = vtirestiff(Fz, cZ_tire, FWZ_vl, FWZ_vr, FWZ_hl, FWZ_hr)
%% Vertical tire stiffness

        cZ_fl = interp1(Fz,cZ_tire,FWZ_vl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_fr = interp1(Fz,cZ_tire,FWZ_vr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_rl = interp1(Fz,cZ_tire,FWZ_hl,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)
        cZ_rr = interp1(Fz,cZ_tire,FWZ_hr,'linear','extrap'); % [N/m] Interpolated tire stiffness (Reifensteifigkeit interpoliert)       
end


function [FWYmax_fl, FWYmax_fr, FWYmax_rl, FWYmax_rr, FWYmax_f, FWYmax_r] = lat_tireforces(FWZ_vl, FWZ_vr,FWZ_hl, FWZ_hr, GAMMA, TIRparam)
%% Maximum transmissible tire forces in lateral direction

            FWYmax_fl = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_vl,GAMMA,0,TIRparam))); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
            FWYmax_fr = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_vr,GAMMA,0,TIRparam))); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
            FWYmax_rl = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_hl,GAMMA,0,TIRparam))); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
            FWYmax_rr = max(abs(MF52_Fy_cs(0:0.1:12,FWZ_hr,GAMMA,0,TIRparam))); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)
           
            FWYmax_f = FWYmax_fl + FWYmax_fr;    % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
            FWYmax_r = FWYmax_rl + FWYmax_rr;    % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
end


function [FWXmax_fl, FWXmax_fr, FWXmax_rl, FWXmax_rr, FWXmax_f, FWXmax_r] = longi_tireforces(FWZ_vl, FWZ_vr,FWZ_hl, FWZ_hr, GAMMA, TIRparam)
%% Maximum transmissible tire forces in longitudinal direction
    
        FWXmax_fl = max(abs(MF52_Fx_cs(0,FWZ_vl,GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximum transmissible front left wheel force (Maximal übertragbare Radkraft)
        FWXmax_fr = max(abs(MF52_Fx_cs(0,FWZ_vr,GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximum transmissible front right wheel force (Maximal übertragbare Radkraft)
        FWXmax_rl = max(abs(MF52_Fx_cs(0,FWZ_hl,GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximum transmissible rear left wheel force (Maximal übertragbare Radkraft)
        FWXmax_rr = max(abs(MF52_Fx_cs(0,FWZ_hr,GAMMA,0:0.01:0.2,TIRparam))); % [N] Maximum transmissible rear right wheel force (Maximal übertragbare Radkraft)        
        
        FWXmax_f = FWXmax_fl + FWXmax_fr;                                % [N] Maximum transmissible front axle force (Maximal übertragbare Achskraft)
        FWXmax_r = FWXmax_rl + FWXmax_rr;                                % [N] Maximum transmissible rear axle force (Maximal übertragbare Achskraft)
end


function [x_Track, y_Track, z_Track, s, R, Track, ApexIndexes, lapLength] = loadTrack(TrackFileName, disciplineID, numOfLaps)
%% Loads the track for the simulation    

    % Loads all the track variables from the Trackfile
    load(TrackFileName,'Track')
    
    % Length of the track before adding laps for endurance
    lapLength = length(Track);

    % Checks if AutoX or Endurance is selected (1==AutoX,2==Endurance)
    if disciplineID == 1    % AutoX
        
        x_Track = Track(:,1);           % [m] X-Coordinate of the Track
        y_Track = Track(:,2);           % [m] Y-Coordinate of the Track
        z_Track = Track(:,3);           % [m] Z-Coordinate of the Track
        s = Track(:,4);                 % [m] Track Pathway (Verlauf der Streckenlänge)
        R = Track(:,5);                 % [m] Radius of curves (Kurvenradien)

        % Calls the .m file 'Apexes' to calculate the Apexes of the given track
        ApexIndexes = Apexes(abs(R));  
        
    elseif disciplineID == 2    % Endurance
        x_Track = [];
        y_Track = [];
        z_Track = [];
        s = [];
        R = [];
        
        s_max = max(Track(:,4));
        
        % Expand Track to Endurance Distance
        for i = 1:numOfLaps
            x_Track = [x_Track; Track(:,1)];     % [m] X-Coordinate of the Track
            y_Track = [y_Track; Track(:,2)];     % [m] Y-Coordinate of the Track
            z_Track = [z_Track; Track(:,3)];     % [m] Z-Coordinate of the Track
            s = [s; Track(:,4)+s_max*(i-1)];                 % [m] Track Pathway (Verlauf der Streckenlänge)
            R = [R; Track(:,5)];                 % [m] Radius of curves (Kurvenradien)

            % Calls the .m file 'Apexes' to calculate the Apexes of the given track
            ApexIndexes = Apexes(abs(R));  
        end
        
        % Complete Endurance Track (All Laps)
        Track = [x_Track, y_Track, z_Track, s, R];
    end

     
end