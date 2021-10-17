%% Vehiclesim_Endurance_GUI_Version.m
% Main function of the simulation. 
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function [result] = Vehiclesim_Endurance_GUI_Version(startingParameters, minValue, minValue2, sensitivityID, sensitivityID2)

    if nargin == 5
        startingParameters.minValue = minValue;
        startingParameters.minValue2 = minValue2;
        startingParameters.sensitivityID = sensitivityID;
        startingParameters.sensitivityID2 = sensitivityID2;
    end
    
    %% Load Setup
    % Adds search path of the setup file and then loads the setup. -> Setup
    % can be in any folder this way.
    setup = load(startingParameters.path+"/"+startingParameters.carDatafile, '-mat');

    tic         % Start timing
    
    if startingParameters.textAreaHandle ~= 0
        %% Initalize GUI loading bar
        % Store original button text
        originalButtonText = startingParameters.processDataButtonHandle.Text;
        % When the function ends, return the original button state
        cleanup = onCleanup(@()set(startingParameters.processDataButtonHandle,'Text',originalButtonText,'Icon',''));
        % Change button name to "Processing"
        startingParameters.processDataButtonHandle.Text = 'Simulating...';
        % Put text on top of icon
        startingParameters.processDataButtonHandle.IconAlignment = 'bottom';
        % Create waitbar with same color as button
        wbar = permute(repmat(startingParameters.processDataButtonHandle.BackgroundColor,15,1,200),[1,3,2]);
        % Black frame around waitbar
        wbar([1,end],:,:) = 0;
        wbar(:,[1,end],:) = 0;
        % Load the empty waitbar to the button
        startingParameters.processDataButtonHandle.Icon = wbar;
    end
    
    %% loads the selected track
    [~, ~, ~, s, R, Track, ApexIndexes, lapLength] = loadTrack(startingParameters.TrackFileName, startingParameters.disciplineID, startingParameters.numOfLaps);
    writeToLogfile('loaded Track!', startingParameters.Debug, startingParameters.textAreaHandle);
    
    %% Vehicle Data (Fahrzeugdaten)   
    FG = setup.m_ges*setup.g;           % [N] Force due to weight of the vehicle (Gewichtskraft des Fahrzeugs)
    setup.aero_ph = 1-setup.aero_pv;

    setup.lf = setup.x_cog-setup.x_va;     
%    setup.lf = setup.wheelbase*setup.m_ph/100;                            % [mm] Distance from front axle to CoG (Abstand Vorderachse zu Fahrzeugschwerpunkt)
    setup.lr = setup.wheelbase-setup.lf;                                  % [mm] Distance from rear axle to CoG (Abstand Hinterachse zu Fahrzeugschwerpunkt)

    %% Initalise DRS   
    for i = 1:length(R)
        if R(i) > setup.DRS_Radius && setup.DRS
            DRS_status(i) = 1;
        else
            DRS_status(i) = 0;
        end
    end      
    
    %% Environmental Conditions (Umgebungsbedingungen)
    rho_L = setup.p_L/(setup.R_L*(setup.t_L+273.15));    % [kg/m�] Air density (p_L in bar)

    %% Motor Data (Motordaten)
    n = setup.engine_param(:,1);
    M = setup.engine_param(:,2);   
    
    max_power = setup.p_max * 1000;         % [W] Power Limit in endurance 
    
    if setup.ptype
        eta_inv = setup.invertor_eff;       % [-] Inverter efficiency
    else
        eta_inv = 1;
    end

    gr = setup.i_param(:);                  % Gear ratio for each gear  
    t_x = 0;                                % Initialise the shift dead time
    
    % Sets gearratio to a constant 1 if no gearbox is present
    if ~setup.gearbox
        gr(1) = 1;
    end
    
    i_G = setup.z_chaindrive / setup.z_sprocket * setup.i_P;          % [-] Gear ratio (Motor to wheel)
    
    % HV-Accumulator Data
    ncells_parallel = setup.nZellen_Parallel;       % [-] Number of parallel cell rows 
    capacity_singlecell = setup.capacity_cell;      % [Wh] Capacity of single cell
    number_cells_in_a_row = setup.nZellen_Reihe;    % [-] Number of cells in one row
    capacity_accumulator = number_cells_in_a_row*ncells_parallel*capacity_singlecell; % [Wh] Capacity of total battery
    
    writeToLogfile('loaded Car Data!', startingParameters.Debug, startingParameters.textAreaHandle);

    %% Initialisation of all variables
    
    load('Emraxefficiencydata_interpolated.mat');
    load('CorrectedDischargeInterpolated.mat');
    load('RandomizedCellData.mat');
    load('CellparametersVoltageInterpolation.mat'); 

    %% Tire Model - Magic Tire Formula Model 5.2 (Reifenmodell - Magic Tire Formula Model 5.2)
    
    load('VerticalStiffness_65kPA_IA0.mat');    % Vertical tire stiffness Lookup-Table (Lookup-Table f�r vertikale Reifensteifigkeit)

    % Load tir-file in structure (Tir-File in Struktur laden)
    TIRparam = loadTIR('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

    % Tire Data (Reifendaten)
    R0 = TIRparam.UNLOADED_RADIUS;      % [m] Tire radius - Manufacturing (Fertigungsradius des Reifens)
    bW = TIRparam.WIDTH;                % [m] Tire width (contact area) (Breite des Reifens (Aufstandsbreite))
    p_infl = setup.p_Tire;              % [Pa] Tire pressure (Luftdruck des Reifens)
    GAMMA = setup.camber;               % [�] Camber (Sturz)

    % Scaling factors for grip level (Skalierungsfaktoren f�r Grip-Niveau)
    % 0.75 for optimum tire temperature; 0.6 for low tire temperature (0.75 f�r optimale Reifentemperatur, 0.6 f�r niedrige Reifentemperatur)
    TIRparam.LMUX = setup.LMUX;          % [-] Longitudinal scaling factor (Skalierungsfaktor L�ngsrichtung)
    TIRparam.LMUY = setup.LMUY;          % [-] Lateral scaling factor (Skalierungsfaktor Querrichtung)

    % Additional factors for longitudinal/lateral slip interaction (Zus�tzliche Faktoren f�r Wechselwirkung L�ngs/Querschlupf)
    TIRparam.LXAL = 1.2;    % [-] Influence of slip angle on transmissible longitudinal force (Einfluss Schr�glaufwinkel auf �bertragbare L�ngskraft)
    TIRparam.LYKA = 1.2;    % [-] Influence of longitudinal slip on transmissible lateral force (Einfluss L�ngsschlupf auf �bertragbare Querkraft)
    TIRparam.RBX3 = 0;      % [-] Additional factor Fx for combined slip (Zus�tzlicher Faktor f�r combined slip Fx)

    %% Set progress
    % sets button progress (progressbar)
    if startingParameters.textAreaHandle ~= 0
        currentProg = min(round((size(wbar,2)-2)*(0/startingParameters.numSteps)),size(wbar,2)-2); 
        startingParameters.processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 1) = 0.25391; % (royalblue)
        startingParameters.processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 2) = 0.41016;
        startingParameters.processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 3) = 0.87891;
        drawnow; % updates button progress
    end
    
    %% Sensitivity Analysis
    steps = 0;
    
    % Checks if sensitvity analysis with two variables is started or
    % not, if so the the second loop is activated with the same length
    % as the outer loop.
    if startingParameters.sensitivityID2 ~= 0 
        numSteps2 = startingParameters.numSteps;
    else
        numSteps2 = 1;
    end
    

    %% For used for sensitivity analysis
    for steps1 = 1:startingParameters.numSteps 
        
        switch startingParameters.sensitivityID
            case 1
                setup.A = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 2
                setup.FB = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 3
                setup.TIRparam.LMUX = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 4
                setup.TIRparam.LMUY = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 5
                setup.aero_pv = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 6
                setup.c_l = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 7
                setup.c_w = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 8
                setup.GAMMA = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 9
                setup.downforce_multiplier = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 10
                setup.drivetrain_eff = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 11
                setup.m_balast = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 12
                setup.m_driver = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 13
                setup.m_ges = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 14
                setup.m_ph = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 15
                setup.setup.n_max = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 16
                setup.num_motors = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 17
                setup.p_infl = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 18
                setup.max_power = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 19
                setup.t_L = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
                setup.rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m�] Air density (Luftdichte)
            case 20
                setup.thetaV_X = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 21
                setup.thetaV_Y = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 22
                setup.thetaV_Z = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 23
                setup.track = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 24
                setup.trq_multiplier = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 25
                setup.wheelbase = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 26
                setup.h_cog = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 27
                setup.x_COG = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 28
                setup.h_cog_balast = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 29
                setup.h_cog_driver = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 30
                setup.x_vA = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 31
                setup.z_chaindrive = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
                setup.setup.z_sprocket = setup.z_sprocket;        % [-] Number of teeth on pinion (Z�hnezahl des Ritzels)
                setup.i_G = z_chaindrive / setup.z_sprocket * i_P;     % [-] Gear ratio (Motor to wheel) (�bersetzung Motor zu Rad)
            case 32
                setup.setup.z_sprocket = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
                setup.z_chaindrive = setup.z_chaindrive; % [-] Number of teeth on sprocket (Z�hnezahl des Kettenblatts)
                setup.i_G = z_chaindrive / setup.z_sprocket * i_P;     % [-] Gear ratio (Motor to wheel) (�bersetzung Motor zu Rad)
            case 33       
                setup.i_G = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 34
                setup.ConstantDownforce = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 35
                setup.DRS_Radius = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
        end
        
        
        % Second for loop for 2 variable sensitivity analysis
        for steps2 = 1:numSteps2
            steps = steps + 1;
        
           
            if startingParameters.sensitivityID2 ~= 0 
                switch startingParameters.sensitivityID2
                    case 1
                        setup.A = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 2
                        setup.FB = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 3
                        setup.TIRparam.LMUX = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 4
                        setup.TIRparam.LMUY = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 5
                        setup.aero_pv = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 6
                        setup.c_l = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 7
                        setup.c_w = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 8
                        setup.GAMMA = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 9
                        setup.downforce_multiplier = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 10
                        setup.drivetrain_eff = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 11
                        setup.m_balast = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 12
                        setup.m_driver = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 13
                        setup.m_ges = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 14
                        setup.m_ph = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 15
                        setup.setup.n_max = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 16
                        setup.num_motors = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 17
                        setup.p_infl = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 18
                        setup.max_power = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 19
                        setup.t_L = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                        setup.rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m�] Air density (Luftdichte)
                    case 20
                        setup.thetaV_X = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 21
                        setup.thetaV_Y = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 22
                        setup.thetaV_Z = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 23
                        setup.track = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 24
                        setup.trq_multiplier = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 25
                        setup.wheelbase = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 26
                        setup.h_cog = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 27
                        setup.x_COG = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 28
                        setup.h_cog_balast = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 29
                        setup.h_cog_driver = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 30
                        setup.x_vA = startingParameters.minValue2 + stepSize*(steps2-1);
                    case 31
                        setup.z_chaindrive = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                        setup.setup.z_sprocket = setup.z_sprocket;        % [-] Number of teeth on pinion (Z�hnezahl des Ritzels)
                        setup.i_G = z_chaindrive / setup.z_sprocket;     % [-] Gear ratio (Motor to wheel) (�bersetzung Motor zu Rad)
                    case 32
                        setup.setup.z_sprocket = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                        setup.z_chaindrive = setup.z_chaindrive; % [-] Number of teeth on sprocket (Z�hnezahl des Kettenblatts)
                        setup.i_G = z_chaindrive / setup.z_sprocket;     % [-] Gear ratio (Motor to wheel) (�bersetzung Motor zu Rad)
                    case 33       
                        setup.i_G = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 34
                        setup.ConstantDownforce = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 35
                        setup.DRS_Radius = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                end
            end    
            
            % Initalize simulation start values
            [aRev,aVX,aVY,BPPsignal,cZ_fl,cZ_fr,cZ_rl,cZ_rr,dFWZrl_aero,dFWZrr_aero,dFWZfl_aero,dFWZfr_aero,dFWZrl_x,dFWZrr_x,dFWZfl_x,dFWZfr_x...
            ,dFWZrl_y,dFWZrr_y,dFWZfl_y,dFWZfr_y,E_Accu,E_Accu_Recu,E_heat,Faero,FB,FVX,FVX_f,FVX_fr,FVX_fl,FVX_rl,FVX_rr,FR,FL,Fdr,FVY,...
            FWXmax_f,FWXmax_r,FWXmax_fl,FWXmax_fr,FWXmax_rl,FWXmax_rr,FWYf,FWYr,FWYmax_f,FWYmax_r,FWYmax_fl,FWYmax_fr,FWYmax_rl,FWYmax_rr,FWZtot,FWZ_rl,...
            FWZ_rr,FWZ_fl,FWZ_fr,FWZr,FWZf,Mi,M_tractive,ni,P_M,P_el,P_tractive,Rdyn_fl,Rdyn_fr,Rdyn_rl,Rdyn_rr,slipY_f,slipY_r,t,l_contact_patch_fl,...
            l_contact_patch_fr,l_contact_patch_rl,l_contact_patch_rr,kappa_rl,kappa_rr,kappa_fl,kappa_fr,delta,beta,psi1,alpha_f,alpha_r,TC,TC_front,ABS,...
            gearSelection,Tirelimit,vAPEXmax,vV,vVYmax,A_accu_cell,P_Mloss,P_Bh,motor_eff,Capacity_Cellpack,SOC_Cellpack,Voltage_Cellpack,V_i,VirtualCurrent_Cellpack,...
            Current_Cellpack_Pointer,Energy_Cellpack,Energy_Cellpack_Total,SOC_Pointer,Current_Cellpack_Pointer_Voltage] = initializeStartValues(setup.FB, Track, ApexIndexes);
            
            alpha_fr = zeros(1,length(Track)-1);
            alpha_rl = zeros(1,length(Track)-1);
            alpha_rr = zeros(1,length(Track)-1);
            alpha_rl = zeros(1,length(Track)-1);
        
            vAPEXmax = zeros(1,length(ApexIndexes));
        
            % Axle and wheel loads (static) ((Statische) Achs- und Radlasten)
            FWZtot(1) = FG;                 % [N] Static total axle load (Statische Gesamtachslast)
            FWZr(1) = setup.m_ph/100*FWZtot(1);   % [N] Static rear axle load (Statische Achslast hinten)
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

            % Static tire loads for calculation (Statische Reifenlasten f�r Berechnung)
            FWZ_fl_stat = FWZ_fl(1); 
            FWZ_fr_stat = FWZ_fr(1);
            FWZ_rl_stat = FWZ_rl(1);
            FWZ_rr_stat = FWZ_rr(1);
            
            % Calculate Skidpad Time and Speed
            [t_skidpad, vV_skidpad] = calculateSkidPad(setup.downforce_multiplier, setup.c_l, setup.A, rho_L, setup.ConstantDownforce, setup.c_l_DRS, DRS_status, setup.m_ges, setup.lr, setup.lf, setup.wheelbase, setup.track, setup.aero_ph, setup.aero_pv, setup.h_cog, GAMMA, TIRparam, FWZ_fl_stat, FWZ_fr_stat, FWZ_rl_stat, FWZ_rr_stat);
            
            %% Calculation of the maximum apex speed for all apexes (numerically) (Berechnen der maximalen Kurvengeschwindigkeiten f�r alle Apexes (numerisch))
            for i = 1:length(ApexIndexes)

                FWYf(i) = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
                FWYr(i) = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
                FWYmax_f(i) = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal �bertragbare Querkraft Vorderachse)
                FWYmax_r(i) = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal �bertragbare Querkraft Hinterachse)
                vV(i) = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)

                while  FWYf(i) < FWYmax_f(i) && FWYr(i) < FWYmax_r(i) && vV(i) < 40

                    vV(i) = vV(i) + 0.01;   % [m/s] Increaing vehicle speed (Erh�hen der Fahrzeuggeschwindigkeit)
              
                    Faero(i) = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, rho_L, vV(i), setup.ConstantDownforce, setup.c_l_DRS, DRS_status(ApexIndexes(i))); % [N] Aerodynamic force

                    FVY(i) = setup.m_ges*vV(i)^2/R(ApexIndexes(i));    % [N] Centrifugal force (Zentrifugalkraft)

                    %aVY(i) = vV(i)^2/R(ApexIndexes(i));  % [m/s�] Lateral acceleration (Querbeschleunigung)

                    % Lateral forces to be applied on front and rear axle (Aufzubringende Querkr�fte an Vorder- und Hinterachse)
                    FWYf(i) = setup.lr/setup.wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                    FWYr(i) = setup.lf/setup.wheelbase*abs(FVY(i));   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

                    % Wheel load transfer due to aero forces (Radlastverlagerung in Folge von Aerokr�ften) 
                    [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = calculateAeroforceOnWheels(Faero(i), setup.aero_ph, setup.aero_pv);

                    % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in L�ngsrichtung = 0 angenommen)
                    [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, 0, aVX(i), setup.wheelbase); % Loads = 0 assumed

                    % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                    [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track, setup.lr, setup.lf, setup.wheelbase, FVY(i));

                    % Wheel loads (Radlasten)
                    FWZ_fl(i) = FWZ_fl_stat + dFWZfl_aero(i) + dFWZfl_x(i) + dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                    FWZ_fr(i) = FWZ_fr_stat + dFWZfr_aero(i) + dFWZfr_x(i) + dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                    FWZ_rl(i) = FWZ_rl_stat + dFWZrl_aero(i) + dFWZrl_x(i) + dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                    FWZ_rr(i) = FWZ_rr_stat + dFWZrr_aero(i) + dFWZrr_x(i) + dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)   

                    % Maximum transmissible tire forces in longitudinal direction = 0 assumed (because longitudinal wheel loads = 0 assumed) 
                    
                    % Calculate delta, beta, psi1 and alpha for all wheels
                    % and front / rear
                    [delta(i), beta(i), psi1(i), alpha_f(i), alpha_r(i), alpha_fr(i), alpha_fl(i), alpha_rr(i), alpha_rl(i)] = calculateSteeringData(setup.wheelbase, R(ApexIndexes(i)), setup.lr, setup.lf, vV(i), FWZ_fl(i), FWZ_rl(i));            

                    % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)    
                    [FWYmax_f(i), FWYmax_r(i)] = calculateLatTireforces(FWZ_fl(i), FWZ_fr(i),FWZ_rl(i), FWZ_rr(i), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

                end

                vAPEXmax(i) = vV(i);   % [m/s] Maximum speed for any apex (Maximalgeschwindigkeit f�r jede Apex)
            end

            writeToLogfile('caclculated Apex Speeds!', startingParameters.Debug, startingParameters.textAreaHandle);
            
            %% Start/Initial values for first simulation run WITHOUT BRAKES (Startwerte f�r ersten Simulationslauf OHNE BREMSEN)
            vV(1) = startingParameters.startingSpeed + 0.005;
            
            gear = 1;

            t(1) = 0;      % [s] Time (Zeit)

            % Supporting variables (Hilfsgr��en)
            z = 1;         % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

            %% Simulation WITHOUT BRAKES (Simulation OHNE BREMSEN)
            for i = 1:length(Track)-1            
                
                if i > 1    % if i > 1 use real rpm instead of idle rpm
                    [ni(i), gear, t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, setup.n_max, t_x, ni(i-1), t(i), t(i-1));      % Calculates Gearbox data and rpm
                else
                    [ni(i), gear, t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, setup.n_max, t_x);      % Calculates Gearbox data and rpm
                end 
                    
                % Calculation of aero forces    
                Faero(i) = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, rho_L, vV(i), setup.ConstantDownforce, setup.c_l_DRS, DRS_status(i)); % [N] Aerodynamic force

                % [Nm] Interpolated motor torque (Motormoment interpoliert)
                Mi(i) = interp1(n,M,ni(i),'linear','extrap'); 

                % Pointer for efficiency table (Pointer f�r effizienztabelle)
                rpmpointer = round(ni(i));                          

                if Mi(i) <= 0
                    torquepointer = 1;
                else
                    torquepointer = round(Mi(i));               % Pointer for efficiency table (Pointer f�r effizienztabelle)
                end

                % Motor power & limitation to 80 kW from FS-Rules (for electric cars) (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
                P_M(i) = setup.num_motors * Mi(i) * ni(i) / 60 * 2 * pi;% [W] Total motor power (Gesamt-Motorleistung)
                if setup.ptype && P_M(i) > max_power
                    P_M(i) = max_power;                           % [W] Limited power (Begrenzte Leistung)
                    %Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Begrenzen des Moments)
                end

                if(rpmpointer > setup.n_max)
                    rpmpointer = setup.n_max;
                elseif(rpmpointer < 1)
                    rpmpointer = 1;
                end

                % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
                if setup.ptype
                    motor_eff(i) = M_eff_inter(rpmpointer,torquepointer);
                else
                    motor_eff(i) = 1;
                end

                P_Mloss(i) = P_M(i)*(1-(motor_eff(i)*setup.drivetrain_eff*eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)

                P_M(i) = P_M(i) - P_Mloss(i);  % Calculation of motor power after deduction of efficiency of the inverter
                
                % Calculation Overall Torque with real power
                Mi(i) = P_M(i)*(60/ni(i)/2/pi);  
                
                % Calculate the tractive forces on the wheels
                [FVX_fl(i), FVX_fr(i), FVX_rl(i), FVX_rr(i), FVX(i), FVX_f(i), TC_front(i), TC(i)] = calculateTractiveForces(Mi(i), setup.num_motors, i_G, gr, Rdyn_fl(i), Rdyn_fr(i), Rdyn_rl(i), Rdyn_rr(i), t_x, gear, FWXmax_f(i), FWXmax_r(i), setup.t_shift);

                % Driving resistances (Fahrwiderst�nde) & Vehicle (Fahrzeug)        
                [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = calculateVehicleResistancesForces(setup.k_R, FWZtot(i), rho_L, vV(i), setup.c_w, setup.A, setup.m_ges, R(i), FVX(i), FVX_f(i), setup.c_d_DRS, DRS_status(i), rpmpointer, setup.n_max, 0, 0, 0, 0);

                if ismember(i,ApexIndexes)
                    if vV(i) > vAPEXmax(z)   % Limiting maximum speeds at apexes (Begrenzen auf maximale Kurvengeschwindigkeit in Apexes)
                        vV(i) = vAPEXmax(z);
                    end
                    z = z + 1;
                end

                vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);            % [s] Time (Zeit)

                % Lateral forces to be applied on front and rear axles (Aufzubringende Querkr�fte an Vorder- und Hinterachse)
                FWYf(i) = setup.lr/setup.wheelbase*FVY(i);   % [N] Lateral force to be applied on front axle (Aufzubringende Querkraft der Vorderachse)
                FWYr(i) = setup.lf/setup.wheelbase*FVY(i);   % [N] Lateral force to be applied on rear axle (Aufzubringende Querkraft der Hinterachse)
                
                % Calculate delta, beta, psi1 and alpha for all wheels
                % and front / rear
                [delta(i), beta(i), psi1(i), alpha_f(i), alpha_r(i), alpha_fr(i), alpha_fl(i), alpha_rr(i), alpha_rl(i)] = calculateSteeringData(setup.wheelbase, R(i), setup.lr, setup.lf, vV(i), FWZ_fl(i), FWZ_rl(i));       

                % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokr�ften)  
                [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = calculateAeroforceOnWheels(Faero(i), setup.aero_ph, setup.aero_pv);

                % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerungen in L�ngsrichtung)
                [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, aVX(i), setup.wheelbase);

                % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track, setup.lr, setup.lf, setup.wheelbase, FVY(i));

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

                % Maximum transmissible tire forces in longitudinal direction (Maximal �bertragbare Reifenkr�fte in L�ngsrichtung)
                [FWXmax_fl(i+1), FWXmax_fr(i+1), FWXmax_rl(i+1), FWXmax_rr(i+1), FWXmax_f(i+1), FWXmax_r(i+1)] = calculateLongiTireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

                % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)
                [FWYmax_fl(i+1), FWYmax_fr(i+1), FWYmax_rl(i+1), FWYmax_rr(i+1), FWYmax_f(i+1), FWYmax_r(i+1)] = calculateLatTireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

            end
            
            vWoBrake = vV;                                      % Save velocity without braking for log file.
            avXWoBrake = aVX;                                   % Save longitudinal acceleration for log file.
            avYWoBrake = aVY;                                   % Save lateral acceleration for log file.

            writeToLogfile('Simulated without brakes!', startingParameters.Debug, startingParameters.textAreaHandle);

            %% BRAKING POINT CALCULATION (BREMSPUNKTBERECHNUNG)
            [BrakeIndexes, NonBrakeApexes, vRev] = calculateBrakepoints(FB, Track, ApexIndexes, vAPEXmax, setup.m_ges, setup.downforce_multiplier, setup.c_l, setup.c_w, setup.A, rho_L, setup.ConstantDownforce, setup.c_l_DRS, DRS_status, setup.aero_ph, setup.aero_pv, vV, setup.k_R, FG, setup.h_cog, setup.wheelbase, setup.track, setup.lr, setup.lf, GAMMA, TIRparam, FWZ_fl_stat, FWZ_fr_stat, FWZ_rl_stat, FWZ_rr_stat, R, s, setup.brakeBias_setup, startingParameters.brakeFunction);

            %% Start values for simulation WITH BRAKES
            vV(1) = startingParameters.startingSpeed + 0.005;
            
            t_x = 0;        % Reset Shift time
            
            gear = 1;
            
            t(1) = 0;       % [s] Time (Zeit)

            E_Accu(1) = 0;       % [J] Energy consumed by battery (Verbrauchte Energie Akku)
            E_heat(1) = 0; 
            E_Accu_Recu(i) = 0;  % [J] Energy recuperated by battery (Rekuperierte Energie Akku)

            % Supporting variables (Hilfsgr��en)
            z = 1;          % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

            %% SIMULATION WITH BRAKES (SIMULATION MIT BREMSEN)
            for i = 1:length(Track)-1

                % Checking at which apex the vehicle is (�berpr�fen, vor welcher Apex das Auto ist)
                if ismember(i,ApexIndexes)  
                    z = z + 1;
                end

                % Saving the current gear for the result file
                gearSelection(i) = gear;
                
                if i > 1    % if i > 1 use real rpm instead of idle rpm
                    [ni(i), gear, t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, setup.n_max, t_x, ni(i-1), t(i), t(i-1));      % Calculates Gearbox data and rpm
                else
                    [ni(i), gear, t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, vV(i), gr, gear, Rdyn_rl(i), Rdyn_rr(i), i_G, setup.n_max, t_x);      % Calculates Gearbox data and rpm
                end 

                [Faero(i)] = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, rho_L, vV(i), setup.ConstantDownforce, setup.c_l_DRS, DRS_status(i)); % [N] Aerodynamic force
                
               %% Braking
                % Checking if braking is required (Pr�fen, ob gebremst werden muss)
                if ismember(i,BrakeIndexes) && not(ismember(z,NonBrakeApexes))                      % Initiaion of braking process (Einleiten des Bremsvorgangs)     
                    % Braking 
                    Mi(i) = 0;                                    % [Nm] Motor torque (Motormoment)
                    BPPsignal(i) = 1;                             % [-] Brake signal (Bremssignal)
                    P_Bh(i) = FWZr(i)/FWZtot(i)*FB(i)*vV(i);      % [W] Rear braking power for recuperation (Bremsleistung hinten f�r Rekuperation)

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
                    
                    [FB_fl(i), FB_fr(i), FB_rl(i), FB_rr(i), ~, BrakeBias(i), ABS(i)] = calculateDeceleration(FB(i), setup.m_ges, Fdr(i), FWXmax_fl(i), FWXmax_fr(i), FWXmax_rl(i), FWXmax_rr(i), setup.brakeBias_setup);
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
                    rpmpointer = round(ni(i));                   % Pointer for efficiency table (Pointer f�r effizienztabelle)

                    if Mi(i) <= 0
                        torquepointer = 1;
                    else
                        torquepointer = round(Mi(i));               % Pointer for efficiency table (Pointer f�r effizienztabelle)
                    end 

                    P_M(i) = setup.num_motors * Mi(i) * ni(i) / 60 * 2 * pi; % [W] Total motor power (P_el!)

                    % Limiting the maximal power when using an electric
                    % drivetrain
                    if setup.ptype && P_M(i) > max_power
                        P_M(i) = max_power;                           % [W] Limited power (P_el!)
                        %Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Total Motor Torque!)
                    end

                    % Checks if the rpmpointer is higher than the maximum rpm and
                    % adjusts it if needed
                    if(rpmpointer > setup.n_max)
                        rpmpointer = setup.n_max;
                    elseif(rpmpointer < 1)
                        rpmpointer = 1;
                    end

                    % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
                    if setup.ptype
                        motor_eff(i) = M_eff_inter(rpmpointer,torquepointer);
                    else
                        motor_eff(i) = 1;
                    end

                    P_Mloss(i) = P_M(i)*(1-(motor_eff(i)*setup.drivetrain_eff*eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)

                    %P_el(i) = P_M(i);

                    % Calculation of motor power after deduction of efficiency of the inverter (ALL MOTORS!)
                    P_M(i) = P_M(i) - P_Mloss(i);  

                    % Calculation Overall Torque with real power
                    Mi(i) = P_M(i)*(60/ni(i)/2/pi);    

                    % Calculate the tractive forces on the wheels
                    [FVX_fl(i), FVX_fr(i), FVX_rl(i), FVX_rr(i), FVX(i), FVX_f(i), TC_front(i), TC(i)] = calculateTractiveForces(Mi(i), setup.num_motors, i_G, gr, Rdyn_fl(i), Rdyn_fr(i), Rdyn_rl(i), Rdyn_rr(i), t_x, gear, FWXmax_f(i), FWXmax_r(i), setup.t_shift);

                    M_tractive(i) = (FVX(i)+FVX_f(i))/(i_G*gr(gear)/Rdyn_rr(i));            % [Nm] Torque including tractive force
                    P_tractive(i) = M_tractive(i)/(60/ni(i)/2/pi);      % [kW] Motor power required for traction 
                    P_el(i) = (P_tractive(i)/(setup.drivetrain_eff * motor_eff(i) * eta_inv));     % [kW] Motor power including efficiencies
                end

                % Driving resistances (Fahrwiderst�nde) & Vehicle (Fahrzeug)
                [FR(i), FL(i), Fdr(i), FVY(i), aVX(i), aVY(i)] = calculateVehicleResistancesForces(setup.k_R, FWZtot(i), rho_L, vV(i), setup.c_w, setup.A, setup.m_ges, R(i), FVX(i), FVX_f(i), setup.c_d_DRS, DRS_status(i), rpmpointer, setup.n_max, FB_fl(i), FB_fr(i), FB_rl(i), FB_rr(i));

                % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                vV(i+1) = sqrt(vV(i)^2+2*aVX(i)*(s(i+1)-s(i)));     

                % ToDo Check vehicle speed before applying brakes 
                % Limit Braking before Apex if car is allready slower than
                % needed
                if ismember(i,BrakeIndexes) && vV(i+1) < vAPEXmax(z)   % Begrenzen der Geschwindigkeit auf ApexGeschwindigkeit (Bremst solange bis Geschwindigkeiten gleich)
                    vV(i+1) = vAPEXmax(z);                         % [m/s] Total vehicle speed 
                end
                
                if vV(i+1) < min(vAPEXmax) && ismember(i,BrakeIndexes)
                    vV(i+1) = min(vAPEXmax);
                elseif vV(i+1) > vWoBrake(i+1)
                    vV(i+1) = vWoBrake(i+1);
                end

                t(i+1) = t(i)+(s(i+1)-s(i))/vV(i+1);                % [s] Time (Zeit)

                % Battery energy capacity (Energiemenge Akku)
                if (t(i+1)-t(i)) * P_Bh(i) > 16.5*10^3 * (t(i+1)-t(i)) 
                    E_Accu_Recu(i+1) = E_Accu_Recu(i) + 16.5*10^3 * (t(i+1)-t(i)); % [J] 
                else
                    E_Accu_Recu(i+1) = E_Accu_Recu(i) + (t(i+1)-t(i)) * P_Bh(i); % [J]
                end
                
                % Latschl�ngen f�r L�ngsschlupfberechnung nach Carter
                FU_fl(i) = FVX(i)/2-setup.k_R*FWZ_fl(i+1)-FB_fl(i);   % [N] Umfangskr�fte an einem Hinterrad
                FU_fr(i) = FVX(i)/2-setup.k_R*FWZ_fr(i+1)-FB_fr(i);   % [N] Umfangskr�fte an einem Hinterrad
                FU_rl(i) = FVX(i)/2-setup.k_R*FWZ_rl(i+1)-FB_rl(i);   % [N] Umfangskr�fte an einem Hinterrad
                FU_rr(i) = FVX(i)/2-setup.k_R*FWZ_rr(i+1)-FB_rr(i);   % [N] Umfangskr�fte an einem Hinterrad
                
                l_contact_patch_fl(i) = FWZ_fl(i)/(p_infl*bW);   % [m] Latschl�nge vorne links
                l_contact_patch_fr(i) = FWZ_fr(i)/(p_infl*bW);   % [m] Latschl�nge vorne rechts
                l_contact_patch_rl(i) = FWZ_rl(i)/(p_infl*bW);   % [m] Latschl�nge hinten links
                l_contact_patch_rr(i) = FWZ_rr(i)/(p_infl*bW);   % [m] Latschl�nge hinten rechts

                % L�ngsschlupf nach Carter
                my0 = 2;    % [-] Haftreibungsbeiwert   
                kappa_fl(i) = l_contact_patch_fl(i)/(2*Rdyn_fl(i))*my0*(setup.wheelbase-sqrt(1-FU_fl(i)/(my0*FWZ_fl(i))));
                kappa_fr(i) = l_contact_patch_fr(i)/(2*Rdyn_fr(i))*my0*(setup.wheelbase-sqrt(1-FU_fr(i)/(my0*FWZ_fr(i))));
                kappa_rl(i) = l_contact_patch_rl(i)/(2*Rdyn_rl(i))*my0*(setup.wheelbase-sqrt(1-FU_rl(i)/(my0*FWZ_rl(i))));
                kappa_rr(i) = l_contact_patch_rr(i)/(2*Rdyn_rr(i))*my0*(setup.wheelbase-sqrt(1-FU_rr(i)/(my0*FWZ_rr(i))));
                
                % Calculate delta, beta, psi1 and alpha for all wheels
                % and front / rear
                [delta(i), beta(i), psi1(i), alpha_f(i), alpha_r(i), alpha_fr(i), alpha_fl(i), alpha_rr(i), alpha_rl(i)] = calculateSteeringData(setup.wheelbase, R(i), setup.lr, setup.lf, vV(i), FWZ_fl(i), FWZ_rl(i));  

                E_Accu(i+1) = E_Accu(i) + (t(i+1)-t(i)) * P_el(i); % [J] 
                E_heat(i+1) = E_heat(i) + (P_el(i)-P_tractive(i)) * (t(i+1)-t(i));
                %E_Waerme(i+1) = P_el(i)*(1-M_eff_inter(i))+(E_Akku(i) + (t(i+1)-t(i)))*(1-eta_inv);    % [J] Motor losses + 5% for inverter losses; Drivetrain losses with heat (Motorverluste + 5% flat f�r inverterverlsute % Drivetrain Verluste auch W�rme)


                % Lateral forces on front and rear axle (Querkr�fte an Vorder- und Hinterachse)
                FWYf(i) = setup.lr/setup.wheelbase*FVY(i);   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                FWYr(i) = setup.lf/setup.wheelbase*FVY(i);   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

                if FWYf(i) > FWYmax_f(i)
                   slipY_f(i) = 1;  
                end

                 if FWYr(i) > FWYmax_r(i)
                   slipY_r(i) = 1;  
                end

                 % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokr�ften)        
                 [dFWZrl_aero(i), dFWZrr_aero(i), dFWZfl_aero(i), dFWZfr_aero(i)] = calculateAeroforceOnWheels(Faero(i), setup.aero_ph, setup.aero_pv);

                 % Dynamic wheel load displacements in longitudinal direction (Dynamische Radlastverlagerungen in L�ngsrichtung)       
                 [dFWZfl_x(i), dFWZfr_x(i), dFWZrl_x(i), dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, aVX(i), setup.wheelbase);

                 % Dynamic wheel load displacements in lateral direction (Dynamische Radlastverlagerung in Querrichtung)  
                 [dFWZfl_y(i), dFWZfr_y(i), dFWZrl_y(i), dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track, setup.lr, setup.lf, setup.wheelbase, FVY(i));

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

                % Maximum transmissible tire forces in longitudinal direction (Maximal �bertragbare Reifenkr�fte in L�ngsrichtung)       
                [FWXmax_fl(i+1), FWXmax_fr(i+1), FWXmax_rl(i+1), FWXmax_rr(i+1),FWXmax_f(i+1), FWXmax_r(i+1)] = calculateLongiTireforces(FWZ_fl(i+1), FWZ_fr(i+1),FWZ_rl(i+1), FWZ_rr(i+1), GAMMA, TIRparam, alpha_f(i), alpha_r(i));

                % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)      
                [FWYmax_fl(i+1), FWYmax_fr(i+1), FWYmax_rl(i+1), FWYmax_rr(i+1), FWYmax_f(i+1), FWYmax_r(i+1)] = calculateLatTireforces(FWZ_fl(i), FWZ_fr(i),FWZ_rl(i), FWZ_rr(i), GAMMA, TIRparam, alpha_f(i), alpha_r(i));
                
                % Maximum cornering velocity 
                vVYmax(i+1) = sqrt((FWYmax_f(i+1)+FWYmax_r(i+1))*R(i+1)/setup.m_ges); % Calculating maximum possible lateral velocity with given Tire forces [m/s] (inaccuaracy because tire force is based on aero force)

%                 %Akkustr�me
%                 V_i(i) = sum(Voltage_Cellpack(:,i));
% 
%                 % Battery Currents (Akkustr�me)
%                 A_accu_cell(i) = P_el(i) / V_i(i) / ncells_parallel;  
% 
%                 Current_Cellpack_Pointer(i) = P_M(i) / V_i(i) *10 ; %Strombelastung eines %er Parrallel Paketes in 0,1A parameter f�r die berechnung der korrigierten belastung mit h�hren verlusten durch h�here zellstr�me
%                 if Current_Cellpack_Pointer(i) <= 1
%                     Current_Cellpack_Pointer(i)=1;
%                 end
% 
%                 if Current_Cellpack_Pointer(i) >= 1500 %begrenzen des max Zellstromes auf 30A pro Zelle im 5er parralelverbund also 150A
%                     Current_Cellpack_Pointer(i)=1500;
%                 end
% 
%                 VirtualCurrent_Cellpack(i) = CorrectedDischargeInterpolated(1,round(Current_Cellpack_Pointer(i))); %Berechnung der Virtuell h�heren zellstr�me basierend auf den h�heren verlsuten durch h�here Str�me
% 
%                 Energy_Cellpack(i) = (VirtualCurrent_Cellpack(i)*(t(i+1)-t(i))) - ((P_Bh(i)/V_i(i))*(t(i+1)-t(i))) ; %Energieverbrauch in As f�r ein 5erpacket an akkuzellen -> Akkustrom zum zeitpunkt i mal Zeitdifferenz zwischen i und i+1
%                 Energy_Cellpack_Total(i+1) = Energy_Cellpack_Total(i) + Energy_Cellpack(i); % �ber Endurance Run Integrierte Energieverbrauch in As f�r ein 5erpacket an akkuzellen
% 
%                 Capacity_Cellpack(1:131,i+1) =  Capacity_Cellpack(1:131,i)- Energy_Cellpack(i); 
% 
%                 SOC_Cellpack(1:131,i+1) = Capacity_Cellpack(1:131,i)./Capacity_Cellpack(1:131,1); %Berechnung des SOC f�r den n�chsten tick basierend auf der aktuellen cellcapacity und der im n�chsten tick
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
%                 %% @Lukas bitte pr�fen if, wegen sonst auftretender Fehler bei langen Strecken!
%                 if (SOC_Pointer(1:131,i+1)>1001) 
%                     Voltage_Cellpack(1:131,i+1) = Voltage_inter(Current_Cellpack_Pointer_Voltage(1,i),1001);    
%                 else
%                     Voltage_Cellpack(1:131,i+1) = Voltage_inter(Current_Cellpack_Pointer_Voltage(1,i),SOC_Pointer(1:131,i+1));
%                 end
                
             
            end
            
            % Conversion of battery energy capacity (Umrechnen der Energiemengen des Akkus)
            E_Accu = E_Accu/(3.6e6);                % [kWh] Energy consumed by battery per lap (Verbrauchte Akku-Energie je Runde)
            E_heat = E_heat/(3.6e6);                % [kWh] 3.6e6 Joule conversion (3.6e6 Umrechnung Joule)
            E_Accu_Recu = E_Accu_Recu/(3.6e6);      % [kWh] Energy recuperated by batter per lap (Rekuperierte Akku-Energie je Runde)
            E_res = E_Accu - E_Accu_Recu;           % [kWh] Resulting energy consumption per lap (Resultierender Verbrauch je Runde)

            %%  Output of the values (Ausgabe der Werte)
            tEnd = toc;
            
            writeToLogfile(['Simulated with brakes! ' num2str(t(end)) ' s'], startingParameters.Debug, startingParameters.textAreaHandle);         

            t_tot = t(end);

            %% Writing the results to Mat File (Schreiben der Ergebnisse in Mat File)
            
            %% skidpad data
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
            result.motor_eff(:,steps) = motor_eff(:);
            result.gear(:,steps) = gearSelection(:);
            
            %% Data without braking 
            result.vWoBrake(:,steps) = vWoBrake(:);
            result.avXWoBrake(:,steps) = avXWoBrake(:);
            result.avYWoBrake(:,steps) = avYWoBrake(:);
                
            %% Cell Data
            if startingParameters.logCellData
                result.Energy_Cellpack(:,steps) = Energy_Cellpack(:);
                result.VirtualCurrent_Cellpack(:,steps) = VirtualCurrent_Cellpack(:);
                result.V_i(:,steps) = V_i(:);
                result.Capacity_Cellpack(:,steps) = Capacity_Cellpack(:);
            end

            % Car Parameters needed to draw result plots and for setup
            % viewer.
            result.rho_L(:,steps) = rho_L(:);
            result.i_G(:,steps) = i_G(:);

            result.DRS_status(:,steps) = DRS_status(:);

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
            result.alpha_rr(:,steps) = alpha_rr(:);
            result.alpha_rl(:,steps) = alpha_rl(:);

            result.lapLength = lapLength;

            result = catstruct(result, setup, startingParameters);                              % Combine Result and Setup to an single struct
            result.processDataButtonHandle = 0;
            result.textAreaHandle = 0;

            % sets button progress (progressbar)
            if startingParameters.sensitivityID2 ~= 0 && startingParameters.textAreaHandle ~= 0
                currentProg = min(round((size(wbar,2)-2)*(steps/startingParameters.numSteps^2)),size(wbar,2)-2); 
            elseif startingParameters.textAreaHandle ~= 0
                currentProg = min(round((size(wbar,2)-2)*(steps/startingParameters.numSteps)),size(wbar,2)-2); 
            end
            
            if startingParameters.textAreaHandle ~= 0
                startingParameters.processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 1) = 0.25391; % (royalblue)
                startingParameters.processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 2) = 0.41016;
                startingParameters.processDataButtonHandle.Icon(2:end-1, 2:currentProg+1, 3) = 0.87891;
                drawnow; % updates button progress
            end
        end
    end
    
    if startingParameters.textAreaHandle ~= 0
        [~, name, ~] = fileparts(startingParameters.carDatafile);
        
        % Saves the results to a .mat file in the same location as the
        % setup file. The name is generated by adding _result to the name.
        savefilename = name + "_result.mat";

        save(startingParameters.path+"/"+savefilename, '-struct','result');
        
        writeToLogfile('File succesfully written', startingParameters.Debug, startingParameters.textAreaHandle);
    end
end