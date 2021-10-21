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
    [~, ~, ~, s, R, Track, sim.ApexIndexes, lapLength] = loadTrack(startingParameters.TrackFileName, startingParameters.disciplineID, startingParameters.numOfLaps);
    writeToLogfile('loaded Track!', startingParameters.Debug, startingParameters.textAreaHandle);
    
    %% Vehicle Data (Fahrzeugdaten)   
    sim.FG = setup.m_ges*setup.g;           % [N] Force due to weight of the vehicle (Gewichtskraft des Fahrzeugs)
    setup.aero_ph = 1-setup.aero_pv;

    setup.lf = setup.x_cog-setup.x_va;     
%    setup.lf = setup.wheelbase*setup.m_ph/100;                            % [mm] Distance from front axle to CoG (Abstand Vorderachse zu Fahrzeugschwerpunkt)
    setup.lr = setup.wheelbase-setup.lf;                                  % [mm] Distance from rear axle to CoG (Abstand Hinterachse zu Fahrzeugschwerpunkt)

    %% Initalise DRS   
    for i = 1:length(R)
        if R(i) > setup.DRS_Radius && setup.DRS
            sim.DRS_status(i) = 1;
        else
            sim.DRS_status(i) = 0;
        end
    end      
    
    %% Environmental Conditions (Umgebungsbedingungen)
    sim.rho_L = setup.p_L/(setup.R_L*(setup.t_L+273.15));    % [kg/m�] Air density (p_L in bar)

    %% Motor Data (Motordaten)
    sim.n = setup.engine_param(:,1);
    sim.M = setup.engine_param(:,2);   
    
    %setup.p_max * 1000;         % [W] Power Limit in endurance 
    
    if setup.ptype
        sim.eta_inv = setup.invertor_eff;       % [-] Inverter efficiency
    else
        sim.eta_inv = 1;
    end

    sim.gr = setup.i_param(:);                  % Gear ratio for each gear  
    sim.t_x = 0;                                % Initialise the shift dead time
    
    % Sets gearratio to a constant 1 if no gearbox is present
    if ~setup.gearbox
        sim.gr(1) = 1;
    end
    
    sim.i_G = setup.z_chaindrive / setup.z_sprocket * setup.i_P;          % [-] Gear ratio (Motor to wheel)
    
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
    TIRparam.LXAL = setup.LXAL;    % [-] Influence of slip angle on transmissible longitudinal force (Einfluss Schr�glaufwinkel auf �bertragbare L�ngskraft)
    TIRparam.LYKA = setup.LYKA;    % [-] Influence of longitudinal slip on transmissible lateral force (Einfluss L�ngsschlupf auf �bertragbare Querkraft)
    TIRparam.RBX3 = setup.RBX3;      % [-] Additional factor Fx for combined slip (Zus�tzlicher Faktor f�r combined slip Fx)

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
                setup.p_max = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
            case 19
                setup.t_L = startingParameters.minValue + startingParameters.stepSize*(steps1-1);
                sim.rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m�] Air density (Luftdichte)
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
                        setup.p_max = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                    case 19
                        setup.t_L = startingParameters.minValue2 + startingParameters.stepSize2*(steps2-1);
                        sim.rho_L = p_L*10^5/(R_L*(t_L+273.15));    % [kg/m�] Air density (Luftdichte)
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
            sim = initializeStartValues(sim, setup.FB, Track, sim.ApexIndexes);

            sim.FB = setup.FB * ones(1,length(Track)-1);
            
            sim.alpha_fr = zeros(1,length(Track)-1);
            sim.alpha_rl = zeros(1,length(Track)-1);
            sim.alpha_rr = zeros(1,length(Track)-1);
            sim.alpha_rl = zeros(1,length(Track)-1);
        
            sim.vAPEXmax = zeros(1,length(sim.ApexIndexes));
        
            % Axle and wheel loads (static) ((Statische) Achs- und Radlasten)
            sim.FWZtot(1) = sim.FG;                 % [N] Static total axle load (Statische Gesamtachslast)
            sim.FWZr(1) = setup.m_ph/100*sim.FWZtot(1);   % [N] Static rear axle load (Statische Achslast hinten)
            sim.FWZf(1) = sim.FWZtot(1)-sim.FWZr(1);    % [N] Static front axle load (Statische Achslast vorne)
            sim.FWZ_fr(1) = sim.FWZf(1)/2;          % [N] Static front right wheel load (Statische Radlast vorne rechts)  
            sim.FWZ_fl(1) = sim.FWZf(1)/2;          % [N] Static front left wheel load (Statische Radlast vorne links)
            sim.FWZ_rr(1) = sim.FWZr(1)/2;          % [N] Static rear right wheel load (Statische Radlast hinten rechts)
            sim.FWZ_rl(1) = sim.FWZr(1)/2;          % [N] Static rear left wheel load (Statische Radlast hinten links)

            % Interpolated vertical tire stiffness (Vertikale Reifensteifigkeiten interpoliert) 
            [sim.cZ_fl(1), sim.cZ_fr(1), sim.cZ_rl(1), sim.cZ_rr(1)] = calculateVtirestiff(Fz, cZ_tire, sim.FWZ_fl(1), sim.FWZ_fr(1), sim.FWZ_rl(1), sim.FWZ_rr(1));

            % Dynamic tire radii (stationary = static) (Dynamische Reifenradien (im Stand = statisch))
            [sim.Rdyn_fl(1), sim.Rdyn_fr(1), sim.Rdyn_rl(1), sim.Rdyn_rr(1)] = calculateDynRadii(R0, sim.FWZ_fl(1), sim.FWZ_fr(1), sim.FWZ_rl(1), sim.FWZ_rr(1), sim.cZ_fl(1), sim.cZ_fr(1), sim.cZ_rr(1), sim.cZ_rl(1));

            sim.FWXmax_r(1) = Inf;

            % Static tire loads for calculation (Statische Reifenlasten f�r Berechnung)
            sim.FWZ_fl_stat = sim.FWZ_fl(1); 
            sim.FWZ_fr_stat = sim.FWZ_fr(1);
            sim.FWZ_rl_stat = sim.FWZ_rl(1);
            sim.FWZ_rr_stat = sim.FWZ_rr(1);
            
            % Calculate Skidpad Time and Speed
            [sim.t_skidpad, sim.vV_skidpad] = calculateSkidPad(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status(i), setup.m_ges, setup.lr, setup.lf, setup.wheelbase, setup.track, setup.aero_ph, setup.aero_pv, setup.h_cog, GAMMA, TIRparam, sim.FWZ_fl_stat, sim.FWZ_fr_stat, sim.FWZ_rl_stat, sim.FWZ_rr_stat);
            
            %% Calculation of the maximum apex speed for all apexes (numerically) (Berechnen der maximalen Kurvengeschwindigkeiten f�r alle Apexes (numerisch))
            for i = 1:length(sim.ApexIndexes)

                sim.FWYf(i) = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
                sim.FWYr(i) = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
                sim.FWYmax_f(i) = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal �bertragbare Querkraft Vorderachse)
                sim.FWYmax_r(i) = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal �bertragbare Querkraft Hinterachse)
                sim.vV(i) = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)

                while  sim.FWYf(i) < sim.FWYmax_f(i) && sim.FWYr(i) < sim.FWYmax_r(i) && sim.vV(i) < 40

                    sim.vV(i) = sim.vV(i) + 0.01;   % [m/s] Increaing vehicle speed (Erh�hen der Fahrzeuggeschwindigkeit)
              
                    sim.Faero(i) = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, sim.vV(i), setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status(sim.ApexIndexes(i))); % [N] Aerodynamic force

                    sim.FVY(i) = setup.m_ges*sim.vV(i)^2/R(sim.ApexIndexes(i));    % [N] Centrifugal force (Zentrifugalkraft)

                    %aVY(i) = vV(i)^2/R(sim.ApexIndexes(i));  % [m/s�] Lateral acceleration (Querbeschleunigung)

                    % Lateral forces to be applied on front and rear axle (Aufzubringende Querkr�fte an Vorder- und Hinterachse)
                    sim.FWYf(i) = setup.lr/setup.wheelbase*abs(sim.FVY(i));   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                    sim.FWYr(i) = setup.lf/setup.wheelbase*abs(sim.FVY(i));   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

                    % Wheel load transfer due to aero forces (Radlastverlagerung in Folge von Aerokr�ften) 
                    [sim.dFWZrl_aero(i), sim.dFWZrr_aero(i), sim.dFWZfl_aero(i), sim.dFWZfr_aero(i)] = calculateAeroforceOnWheels(sim.Faero(i), setup.aero_ph, setup.aero_pv);

                    % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in L�ngsrichtung = 0 angenommen)
                    [sim.dFWZfl_x(i), sim.dFWZfr_x(i), sim.dFWZrl_x(i), sim.dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, 0, 0, setup.wheelbase); % Loads = 0 assumed

                    % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                    [sim.dFWZfl_y(i), sim.dFWZfr_y(i), sim.dFWZrl_y(i), sim.dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track, setup.lr, setup.lf, setup.wheelbase, sim.FVY(i));

                    % Wheel loads (Radlasten)
                    sim.FWZ_fl(i) = sim.FWZ_fl_stat + sim.dFWZfl_aero(i) + sim.dFWZfl_x(i) + sim.dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                    sim.FWZ_fr(i) = sim.FWZ_fr_stat + sim.dFWZfr_aero(i) + sim.dFWZfr_x(i) + sim.dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                    sim.FWZ_rl(i) = sim.FWZ_rl_stat + sim.dFWZrl_aero(i) + sim.dFWZrl_x(i) + sim.dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                    sim.FWZ_rr(i) = sim.FWZ_rr_stat + sim.dFWZrr_aero(i) + sim.dFWZrr_x(i) + sim.dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts) 
                    


                    % Maximum transmissible tire forces in longitudinal direction = 0 assumed (because longitudinal wheel loads = 0 assumed) 
                    
                    % Calculate delta, beta, psi1 and alpha for all wheels
                    % and front / rear
                    [sim.delta(i), sim.beta(i), sim.psi1(i), sim.alpha_f(i), sim.alpha_r(i), sim.alpha_fr(i), sim.alpha_fl(i), sim.alpha_rr(i), sim.alpha_rl(i)] = calculateSteeringData(setup.wheelbase, R(sim.ApexIndexes(i)), setup.lr, setup.lf, sim.vV(i), sim.FWZ_fl(i), sim.FWZ_rl(i),sim.FWZ_fr(i),sim.FWZ_rr(i));            

                    % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)    
                    [sim.FWYmax_f(i), sim.FWYmax_r(i)] = calculateLatTireforces(sim.FWZ_fl(i), sim.FWZ_fr(i), sim.FWZ_rl(i), sim.FWZ_rr(i), GAMMA, TIRparam, sim.alpha_f(i), sim.alpha_r(i));

                end

                sim.vAPEXmax(i) = sim.vV(i);   % [m/s] Maximum speed for any apex (Maximalgeschwindigkeit f�r jede Apex)
            end

            writeToLogfile('caclculated Apex Speeds!', startingParameters.Debug, startingParameters.textAreaHandle);
            
            %% Start/Initial values for first simulation run WITHOUT BRAKES (Startwerte f�r ersten Simulationslauf OHNE BREMSEN)
            sim.vV(1) = startingParameters.startingSpeed + 0.005;
            
            sim.gear = 1;

            sim.t(1) = 0;      % [s] Time (Zeit)

            % Supporting variables (Hilfsgr��en)
            sim.z = 1;         % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)
            
            sim.FWXmax_f(1) = 0;
            sim.FWXmax_f(1) = 0;

            %% Simulation WITHOUT BRAKES (Simulation OHNE BREMSEN)
            for i = 1:length(Track)-1            
                
                if i > 1    % if i > 1 use real rpm instead of idle rpm
                    [sim.ni(i), sim.gear, sim.t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, sim.vV(i), sim.gr, sim.gear, sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.i_G, setup.n_max, sim.t_x, sim.ni(i-1), sim.t(i), sim.t(i-1));      % Calculates Gearbox data and rpm
                else
                    [sim.ni(i), sim.gear, sim.t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, sim.vV(i), sim.gr, sim.gear, sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.i_G, setup.n_max, sim.t_x);      % Calculates Gearbox data and rpm
                end 
                    
                % Calculation of aero forces    
                sim.Faero(i) = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, sim.vV(i), setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status(i)); % [N] Aerodynamic force

                % [Nm] Interpolated motor torque (Motormoment interpoliert)
                sim.Mi(i) = interp1(sim.n,sim.M,sim.ni(i),'linear','extrap'); 

                % Pointer for efficiency table (Pointer f�r effizienztabelle)
                sim.rpmpointer = round(sim.ni(i));                          

                if sim.Mi(i) <= 0
                    sim.torquepointer = 1;
                else
                    sim.torquepointer = round(sim.Mi(i));               % Pointer for efficiency table (Pointer f�r effizienztabelle)
                end

                % Motor power & limitation to 80 kW from FS-Rules (for electric cars) (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
                sim.P_M(i) = setup.num_motors * sim.Mi(i) * sim.ni(i) / 60 * 2 * pi;% [W] Total motor power (Gesamt-Motorleistung)
                if setup.ptype && sim.P_M(i) > setup.p_max
                    sim.P_M(i) = setup.p_max;                           % [W] Limited power (Begrenzte Leistung)
                    %Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Begrenzen des Moments)
                end

                if(sim.rpmpointer > setup.n_max)
                    sim.rpmpointer = setup.n_max;
                elseif(sim.rpmpointer < 1)
                    sim.rpmpointer = 1;
                end

                % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
                if setup.ptype
                    sim.motor_eff(i) = M_eff_inter(sim.rpmpointer,sim.torquepointer);
                else
                    sim.motor_eff(i) = 1;
                end

                sim.P_Mloss(i) = sim.P_M(i)*(1-(sim.motor_eff(i)*setup.drivetrain_eff*sim.eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)

                sim.P_M(i) = sim.P_M(i) - sim.P_Mloss(i);  % Calculation of motor power after deduction of efficiency of the inverter
                
                % Calculation Overall Torque with real power
                sim.Mi(i) = sim.P_M(i)*(60/sim.ni(i)/2/pi);  
                
                % Calculate the tractive forces on the wheels
                [sim.FVX_fl(i), sim.FVX_fr(i), sim.FVX_rl(i), sim.FVX_rr(i), sim.FVX(i), sim.FVX_f(i), sim.TC_front(i), sim.TC(i)] = calculateTractiveForces(sim.Mi(i), setup.num_motors, sim.i_G, sim.gr, sim.Rdyn_fl(i), sim.Rdyn_fr(i), sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.t_x, sim.gear, sim.FWXmax_f(i), sim.FWXmax_r(i), setup.t_shift);

                % Driving resistances (Fahrwiderst�nde) & Vehicle (Fahrzeug)        
                [sim.FR(i), sim.FL(i), sim.Fdr(i), sim.FVY(i), sim.aVX(i), sim.aVY(i)] = calculateVehicleResistancesForces(setup.k_R, sim.FWZtot(i), sim.rho_L, sim.vV(i), setup.c_w, setup.A, setup.m_ges, R(i), sim.FVX(i), sim.FVX_f(i), setup.c_d_DRS, sim.DRS_status(i), sim.rpmpointer, setup.n_max, 0, 0, 0, 0);

                if ismember(i,sim.ApexIndexes)
                    if sim.vV(i) > sim.vAPEXmax(sim.z)   % Limiting maximum speeds at apexes (Begrenzen auf maximale Kurvengeschwindigkeit in Apexes)
                        sim.vV(i) = sim.vAPEXmax(sim.z);
                    end
                    sim.z = sim.z + 1;
                end

                sim.vV(i+1) = sqrt(sim.vV(i)^2+2*sim.aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                sim.t(i+1) = sim.t(i)+(s(i+1)-s(i))/sim.vV(i+1);            % [s] Time (Zeit)

                % Lateral forces to be applied on front and rear axles (Aufzubringende Querkr�fte an Vorder- und Hinterachse)
                sim.FWYf(i) = setup.lr/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied on front axle (Aufzubringende Querkraft der Vorderachse)
                sim.FWYr(i) = setup.lf/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied on rear axle (Aufzubringende Querkraft der Hinterachse)
                
                % Calculate delta, beta, psi1 and alpha for all wheels
                % and front / rear
                [sim.delta(i), sim.beta(i), sim.psi1(i), sim.alpha_f(i), sim.alpha_r(i), sim.alpha_fr(i), sim.alpha_fl(i), sim.alpha_rr(i), sim.alpha_rl(i)] = calculateSteeringData(setup.wheelbase, R(i), setup.lr, setup.lf, sim.vV(i), sim.FWZ_fl(i), sim.FWZ_rl(i), sim.FWZ_fr(i), sim.FWZ_rr(i));       

                % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokr�ften)  
                [sim.dFWZrl_aero(i), sim.dFWZrr_aero(i), sim.dFWZfl_aero(i), sim.dFWZfr_aero(i)] = calculateAeroforceOnWheels(sim.Faero(i), setup.aero_ph, setup.aero_pv);

                % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerungen in L�ngsrichtung)
                [sim.dFWZfl_x(i), sim.dFWZfr_x(i), sim.dFWZrl_x(i), sim.dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, sim.aVX(i), setup.wheelbase);

                % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                [sim.dFWZfl_y(i), sim.dFWZfr_y(i), sim.dFWZrl_y(i), sim.dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track, setup.lr, setup.lf, setup.wheelbase, sim.FVY(i));

                % Wheel loads (Radlasten)
                sim.FWZ_fl(i+1) = sim.FWZ_fl(1) + sim.dFWZfl_aero(i) + sim.dFWZfl_x(i) + sim.dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                sim.FWZ_fr(i+1) = sim.FWZ_fr(1) + sim.dFWZfr_aero(i) + sim.dFWZfr_x(i) + sim.dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                sim.FWZ_rl(i+1) = sim.FWZ_rl(1) + sim.dFWZrl_aero(i) + sim.dFWZrl_x(i) + sim.dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                sim.FWZ_rr(i+1) = sim.FWZ_rr(1) + sim.dFWZrr_aero(i) + sim.dFWZrr_x(i) + sim.dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)

                % Limiting the wheel loads to (almost) zero (Begrenzen der Radlasten auf (quasi) Null)
                if sim.FWZ_fl(i+1) < 0
                    sim.FWZ_fl(i+1) = 0.001;
                end
                if sim.FWZ_fr(i+1) < 0
                    sim.FWZ_fr(i+1) = 0.001;
                end
                if sim.FWZ_rl(i+1) < 0
                    sim.FWZ_rl(i+1) = 0.001;
                end
                if sim.FWZ_rr(i+1) < 0
                    sim.FWZ_rr(i+1) = 0.001;
                end

                % Axle loads - for dynamic radii (Achslasten)
                [sim.FWZr(i+1), sim.FWZf(i+1), sim.FWZtot(i+1)] = calculateAxleloads(sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), sim.FWZ_fl(i+1), sim.FWZ_fr(i+1));

                % Vertical tire stiffnesses - for dynamic radii (Vertikale Reifensteifigkeiten)  
                [sim.cZ_fl(i+1), sim.cZ_fr(i+1), sim.cZ_rl(i+1), sim.cZ_rr(i+1)] = calculateVtirestiff(Fz, cZ_tire, sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1));

                % Dynamic tire radii (Dynamische Reifenradien)
                [sim.Rdyn_fl(i+1), sim.Rdyn_fr(i+1), sim.Rdyn_rl(i+1), sim.Rdyn_rr(i+1)] = calculateDynRadii(R0, sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), sim.cZ_fl(i+1), sim.cZ_fr(i+1), sim.cZ_rr(i+1), sim.cZ_rl(i+1));

                % Maximum transmissible tire forces in longitudinal direction (Maximal �bertragbare Reifenkr�fte in L�ngsrichtung)
                [sim.FWXmax_fl(i+1), sim.FWXmax_fr(i+1), sim.FWXmax_rl(i+1), sim.FWXmax_rr(i+1), sim.FWXmax_f(i+1), sim.FWXmax_r(i+1)] = calculateLongiTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i), sim.alpha_r(i));

                % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)
                [sim.FWYmax_fl(i+1), sim.FWYmax_fr(i+1), sim.FWYmax_rl(i+1), sim.FWYmax_rr(i+1), sim.FWYmax_f(i+1), sim.FWYmax_r(i+1)] = calculateLatTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i), sim.alpha_r(i));

            end
            
            sim.vWoBrake = sim.vV;                                      % Save velocity without braking for log file.
            sim.avXWoBrake = sim.aVX;                                   % Save longitudinal acceleration for log file.
            sim.avYWoBrake = sim.aVY;                                   % Save lateral acceleration for log file.

            writeToLogfile('Simulated without brakes!', startingParameters.Debug, startingParameters.textAreaHandle);

            %% BRAKING POINT CALCULATION (BREMSPUNKTBERECHNUNG)
            [sim.BrakeIndexes, sim.NonBrakeApexes, sim.vRev] = calculateBrakepoints(setup.FB, Track, sim.ApexIndexes, sim.vAPEXmax, setup.m_ges, setup.downforce_multiplier, setup.c_l, setup.c_w, setup.A, sim.rho_L, setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status, setup.aero_ph, setup.aero_pv, sim.vV, setup.k_R, sim.FG, setup.h_cog, setup.wheelbase, setup.track, setup.lr, setup.lf, GAMMA, TIRparam, sim.FWZ_fl_stat, sim.FWZ_fr_stat, sim.FWZ_rl_stat, sim.FWZ_rr_stat, R, s, setup.brakeBias_setup, startingParameters.brakeFunction);

            %% Start values for simulation WITH BRAKES
            sim.vV(1) = startingParameters.startingSpeed + 0.005;
            
            sim.t_x = 0;        % Reset Shift time
            
            sim.gear = 1;
            
            sim.t(1) = 0;       % [s] Time (Zeit)

            sim.E_Accu(1) = 0;       % [J] Energy consumed by battery (Verbrauchte Energie Akku)
            sim.E_heat(1) = 0; 
            sim.E_Accu_Recu(i) = 0;  % [J] Energy recuperated by battery (Rekuperierte Energie Akku)

            % Supporting variables (Hilfsgr��en)
            sim.z = 1;          % [-] Determination of the upcoming apex (Bestimmung der anstehenden Apex)

            %% SIMULATION WITH BRAKES (SIMULATION MIT BREMSEN)
            for i = 1:length(Track)-1

                % Checking at which apex the vehicle is (�berpr�fen, vor welcher Apex das Auto ist)
                if ismember(i,sim.ApexIndexes)  
                    sim.z = sim.z + 1;
                end

                % Saving the current gear for the result file
                sim.gearSelection(i) = sim.gear;
                
                if i > 1    % if i > 1 use real rpm instead of idle rpm
                    [sim.ni(i), sim.gear, sim.t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, sim.vV(i), sim.gr, sim.gear, sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.i_G, setup.n_max, sim.t_x, sim.ni(i-1), sim.t(i), sim.t(i-1));      % Calculates Gearbox data and rpm
                else
                    [sim.ni(i), sim.gear, sim.t_x] = calculateGearbox(setup.gearbox, setup.idleRPM, setup.n_shift, setup.n_downshift, sim.vV(i), sim.gr, sim.gear, sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.i_G, setup.n_max, sim.t_x);      % Calculates Gearbox data and rpm
                end 

                [sim.Faero(i)] = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, sim.vV(i), setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status(i)); % [N] Aerodynamic force
                
               %% Braking
                % Checking if braking is required (Pr�fen, ob gebremst werden muss)
                if ismember(i,sim.BrakeIndexes) && not(ismember(sim.z,sim.NonBrakeApexes))                      % Initiaion of braking process (Einleiten des Bremsvorgangs)     
                    % Braking 
                    sim.Mi(i) = 0;                                    % [Nm] Motor torque (Motormoment)
                    sim.BPPsignal(i) = 1;                             % [-] Brake signal (Bremssignal)
                    sim.P_Bh(i) = sim.FWZr(i)/sim.FWZtot(i)*sim.FB(i)*sim.vV(i);      % [W] Rear braking power for recuperation (Bremsleistung hinten f�r Rekuperation)

                    sim.P_M(i) = 0;
                    sim.Mi(i) = 0;
                    sim.motor_eff(i) = 0;
                    sim.P_Mloss(i) = 0;

                    sim.M_tractive(i) = 0;
                    sim.P_tractive(i) = 0;
                    sim.P_el(i) = 0;

                    sim.FVX_fl(i) = 0;                   % [N] Tractive Force on front left wheel (AWD) 
                    sim.FVX_fr(i) = 0;                   % [N] Tractive Force on front right wheel (AWD) 
                    sim.FVX_f(i) = 0;              % [N] Traction on rear axle (Zugkraft an der Hinterachse)

                    sim.FVX_rl(i) = 0;                   % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
                    sim.FVX_rr(i) = 0;                   % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
                    sim.FVX(i) = 0;                % [N] Traction on rear axle (Zugkraft an der Hinterachse)              
                    
                    [sim.FB_fl(i), sim.FB_fr(i), sim.FB_rl(i), sim.FB_rr(i), ~, sim.BrakeBias(i), sim.ABS(i)] = calculateDeceleration(sim.FB(i), setup.m_ges, sim.Fdr(i), sim.FWXmax_fl(i), sim.FWXmax_fr(i), sim.FWXmax_rl(i), sim.FWXmax_rr(i), setup.brakeBias_setup);
                else
                %% Accelerating 
                
                    sim.Mi(i) = interp1(sim.n,sim.M,sim.ni(i),'linear','extrap'); % [Nm] Motor torque (single motor!)
                    sim.FB(i) = 0;                                    % [N] Braking force
                    sim.P_Bh(i) = 0;                                  % [W] Rear braking power (Bremsleistung hinten)
                    sim.FB_fl(i) = 0;    
                    sim.FB_fr(i) = 0;    
                    sim.FB_rl(i) = 0;    
                    sim.FB_rr(i) = 0;    

                    % Motor power & limitation to 80 kW from FS-Rules (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
                    sim.rpmpointer = round(sim.ni(i));                   % Pointer for efficiency table (Pointer f�r effizienztabelle)

                    if sim.Mi(i) <= 0
                        sim.torquepointer = 1;
                    else
                        sim.torquepointer = round(sim.Mi(i));               % Pointer for efficiency table (Pointer f�r effizienztabelle)
                    end 

                    sim.P_M(i) = setup.num_motors * sim.Mi(i) * sim.ni(i) / 60 * 2 * pi; % [W] Total motor power (P_el!)

                    % Limiting the maximal power when using an electric
                    % drivetrain
                    if setup.ptype && sim.P_M(i) > setup.p_max
                        sim.P_M(i) = setup.p_max;                           % [W] Limited power (P_el!)
                        %Mi(i) = P_M(i)*60/ni(i)/2/pi;                 % [Nm] Limiting the torque (Total Motor Torque!)
                    end

                    % Checks if the rpmpointer is higher than the maximum rpm and
                    % adjusts it if needed
                    if(sim.rpmpointer > setup.n_max)
                        sim.rpmpointer = setup.n_max;
                    elseif(sim.rpmpointer < 1)
                        sim.rpmpointer = 1;
                    end

                    % Motor efficiency at given speed and torque (Motor Effizienz bei Drehzahl und Moment)
                    if setup.ptype
                        sim.motor_eff(i) = M_eff_inter(sim.rpmpointer,sim.torquepointer);
                    else
                        sim.motor_eff(i) = 1;
                    end

                    sim.P_Mloss(i) = sim.P_M(i)*(1-(sim.motor_eff(i)*setup.drivetrain_eff*sim.eta_inv)); % Calculation of power loss (berechnung der Verlustleistung)

                    %P_el(i) = P_M(i);

                    % Calculation of motor power after deduction of efficiency of the inverter (ALL MOTORS!)
                    sim.P_M(i) = sim.P_M(i) - sim.P_Mloss(i);  

                    % Calculation Overall Torque with real power
                    sim.Mi(i) = sim.P_M(i)*(60/sim.ni(i)/2/pi); 

                    % Calculate the tractive forces on the wheels
                    [sim.FVX_fl(i), sim.FVX_fr(i), sim.FVX_rl(i), sim.FVX_rr(i), sim.FVX(i), sim.FVX_f(i), sim.TC_front(i), sim.TC(i)] = calculateTractiveForces(sim.Mi(i), setup.num_motors, sim.i_G, sim.gr, sim.Rdyn_fl(i), sim.Rdyn_fr(i), sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.t_x, sim.gear, sim.FWXmax_f(i), sim.FWXmax_r(i), setup.t_shift);

                    sim.M_tractive(i) = (sim.FVX(i)+sim.FVX_f(i))/(sim.i_G*sim.gr(sim.gear)/sim.Rdyn_rr(i));        % [Nm] Torque including tractive force
                    sim.P_tractive(i) = sim.M_tractive(i)/(60/sim.ni(i)/2/pi);                                      % [kW] Motor power required for tractive forces
                    sim.P_el(i) = (sim.P_tractive(i)/(setup.drivetrain_eff * sim.motor_eff(i) * sim.eta_inv));      % [kW] Motor power required for tractive forces including efficiencies
                end

                % Driving resistances (Fahrwiderst�nde) & Vehicle (Fahrzeug)
                [sim.FR(i), sim.FL(i), sim.Fdr(i), sim.FVY(i), sim.aVX(i), sim.aVY(i)] = calculateVehicleResistancesForces(setup.k_R, sim.FWZtot(i), sim.rho_L, sim.vV(i), setup.c_w, setup.A, setup.m_ges, R(i), sim.FVX(i), sim.FVX_f(i), setup.c_d_DRS, sim.DRS_status(i), sim.rpmpointer, setup.n_max, sim.FB_fl(i), sim.FB_fr(i), sim.FB_rl(i), sim.FB_rr(i));

                % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
                sim.vV(i+1) = sqrt(sim.vV(i)^2+2*sim.aVX(i)*(s(i+1)-s(i)));     

                % ToDo Check vehicle speed before applying brakes 
                % Limit Braking before Apex if car is allready slower than
                % needed
                if ismember(i,sim.BrakeIndexes) && sim.vV(i+1) < sim.vAPEXmax(sim.z)   % Begrenzen der Geschwindigkeit auf ApexGeschwindigkeit (Bremst solange bis Geschwindigkeiten gleich)
                    sim.vV(i+1) = sim.vAPEXmax(sim.z);                         % [m/s] Total vehicle speed 
                end
                
                if ismember(i,sim.BrakeIndexes) && sim.vV(i+1) < min(sim.vAPEXmax)
                    sim.vV(i+1) = min(sim.vAPEXmax);
                elseif sim.vV(i+1) > sim.vWoBrake(i+1)
                    sim.vV(i+1) = sim.vWoBrake(i+1);
                end 

                sim.t(i+1) = sim.t(i)+(s(i+1)-s(i))/sim.vV(i+1);                % [s] Time (Zeit)

                % Battery energy capacity (Energiemenge Akku)
                if (sim.t(i+1)-sim.t(i)) * sim.P_Bh(i) > 16.5*10^3 * (sim.t(i+1)-sim.t(i)) 
                    sim.E_Accu_Recu(i+1) = sim.E_Accu_Recu(i) + 16.5*10^3 * (sim.t(i+1)-sim.t(i)); % [J] 
                else
                    sim.E_Accu_Recu(i+1) = sim.E_Accu_Recu(i) + (sim.t(i+1)-sim.t(i)) * sim.P_Bh(i); % [J]
                end
                
                % Latschl�ngen f�r L�ngsschlupfberechnung nach Carter
                % (WRONG USE FVX for each wheel wrong formula for front wheels!)
                sim.FU_fl(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_fl(i+1)-sim.FB_fl(i);   % [N] Umfangskr�fte an einem Hinterrad
                sim.FU_fr(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_fr(i+1)-sim.FB_fr(i);   % [N] Umfangskr�fte an einem Hinterrad
                sim.FU_rl(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_rl(i+1)-sim.FB_rl(i);   % [N] Umfangskr�fte an einem Hinterrad
                sim.FU_rr(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_rr(i+1)-sim.FB_rr(i);   % [N] Umfangskr�fte an einem Hinterrad
                
                sim.l_contact_patch_fl(i) = sim.FWZ_fl(i)/(p_infl*bW);   % [m] Latschl�nge vorne links
                sim.l_contact_patch_fr(i) = sim.FWZ_fr(i)/(p_infl*bW);   % [m] Latschl�nge vorne rechts
                sim.l_contact_patch_rl(i) = sim.FWZ_rl(i)/(p_infl*bW);   % [m] Latschl�nge hinten links
                sim.l_contact_patch_rr(i) = sim.FWZ_rr(i)/(p_infl*bW);   % [m] Latschl�nge hinten rechts

                % L�ngsschlupf nach Carter
                sim.my0 = 2;    % [-] Haftreibungsbeiwert   
                sim.kappa_fl(i) = sim.l_contact_patch_fl(i)/(2*sim.Rdyn_fl(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_fl(i)/(sim.my0*sim.FWZ_fl(i))));
                sim.kappa_fr(i) = sim.l_contact_patch_fr(i)/(2*sim.Rdyn_fr(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_fr(i)/(sim.my0*sim.FWZ_fr(i))));
                sim.kappa_rl(i) = sim.l_contact_patch_rl(i)/(2*sim.Rdyn_rl(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_rl(i)/(sim.my0*sim.FWZ_rl(i))));
                sim.kappa_rr(i) = sim.l_contact_patch_rr(i)/(2*sim.Rdyn_rr(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_rr(i)/(sim.my0*sim.FWZ_rr(i))));
                
                % Calculate delta, beta, psi1 and alpha for all wheels
                % and front / rear
                [sim.delta(i), sim.beta(i), sim.psi1(i), sim.alpha_f(i), sim.alpha_r(i), sim.alpha_fr(i), sim.alpha_fl(i), sim.alpha_rr(i), sim.alpha_rl(i)] = calculateSteeringData(setup.wheelbase, R(i), setup.lr, setup.lf, sim.vV(i), sim.FWZ_fl(i), sim.FWZ_rl(i),sim.FWZ_fr(i),sim.FWZ_rr(i));  

                sim.E_Accu(i+1) = sim.E_Accu(i) + (sim.t(i+1)-sim.t(i)) * sim.P_el(i); % [J] 
                sim.E_heat(i+1) = sim.E_heat(i) + (sim.P_el(i)-sim.P_tractive(i)) * (sim.t(i+1)-sim.t(i));
                %E_Waerme(i+1) = P_el(i)*(1-M_eff_inter(i))+(E_Akku(i) + (t(i+1)-t(i)))*(1-eta_inv);    % [J] Motor losses + 5% for inverter losses; Drivetrain losses with heat (Motorverluste + 5% flat f�r inverterverlsute % Drivetrain Verluste auch W�rme)


                % Lateral forces on front and rear axle (Querkr�fte an Vorder- und Hinterachse)
                sim.FWYf(i) = setup.lr/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
                sim.FWYr(i) = setup.lf/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

                if sim.FWYf(i) > sim.FWYmax_f(i)
                   sim.slipY_f(i) = 1;  
                end

                 if sim.FWYr(i) > sim.FWYmax_r(i)
                   sim.slipY_r(i) = 1;  
                end

                 % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokr�ften)        
                 [sim.dFWZrl_aero(i), sim.dFWZrr_aero(i), sim.dFWZfl_aero(i), sim.dFWZfr_aero(i)] = calculateAeroforceOnWheels(sim.Faero(i), setup.aero_ph, setup.aero_pv);

                 % Dynamic wheel load displacements in longitudinal direction (Dynamische Radlastverlagerungen in L�ngsrichtung)       
                 [sim.dFWZfl_x(i), sim.dFWZfr_x(i), sim.dFWZrl_x(i), sim.dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, sim.aVX(i), setup.wheelbase);

                 % Dynamic wheel load displacements in lateral direction (Dynamische Radlastverlagerung in Querrichtung)  
                 [sim.dFWZfl_y(i), sim.dFWZfr_y(i), sim.dFWZrl_y(i), sim.dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track, setup.lr, setup.lf, setup.wheelbase, sim.FVY(i));

                 % Wheel Loads (Radlasten)
                 sim.FWZ_fl(i+1) = sim.FWZ_fl(1) + sim.dFWZfl_aero(i) + sim.dFWZfl_x(i) + sim.dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
                 sim.FWZ_fr(i+1) = sim.FWZ_fr(1) + sim.dFWZfr_aero(i) + sim.dFWZfr_x(i) + sim.dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
                 sim.FWZ_rl(i+1) = sim.FWZ_rl(1) + sim.dFWZrl_aero(i) + sim.dFWZrl_x(i) + sim.dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
                 sim.FWZ_rr(i+1) = sim.FWZ_rr(1) + sim.dFWZrr_aero(i) + sim.dFWZrr_x(i) + sim.dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)

                 % Limiting the wheel loads to (almost) zero (Begrenzen der Radlasten auf (quasi) Null)
                 if sim.FWZ_fl(i+1) < 0
                     sim.FWZ_fl(i+1) = 0.001;
                 end
                 if sim.FWZ_fr(i+1) < 0
                     sim.FWZ_fr(i+1) = 0.001;
                 end
                 if sim.FWZ_rl(i+1) < 0
                     sim.FWZ_rl(i+1) = 0.001;
                 end
                 if sim.FWZ_rr(i+1) < 0
                     sim.FWZ_rr(i+1) = 0.001;
                 end

                % Axle loads - for dynamic radii (Achslasten)      
                [sim.FWZr(i+1), sim.FWZf(i+1), sim.FWZtot(i+1)] = calculateAxleloads(sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), sim.FWZ_fl(i+1), sim.FWZ_fr(i+1)); 

                % Vertical tire stiffness - for dynamic radii (Vertikale Reifensteifigkeiten)  
                [sim.cZ_fl(i+1), sim.cZ_fr(i+1), sim.cZ_rl(i+1), sim.cZ_rr(i+1)] = calculateVtirestiff(Fz, cZ_tire, sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1));

                % Dynamic tire radii (Dynamische Reifenradien)   
                [sim.Rdyn_fl(i+1), sim.Rdyn_fr(i+1), sim.Rdyn_rl(i+1), sim.Rdyn_rr(i+1)] = calculateDynRadii(R0, sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), sim.cZ_fl(i+1), sim.cZ_fr(i+1), sim.cZ_rr(i+1), sim.cZ_rl(i+1));

                % Maximum transmissible tire forces in longitudinal direction (Maximal �bertragbare Reifenkr�fte in L�ngsrichtung)       
                [sim.FWXmax_fl(i+1), sim.FWXmax_fr(i+1), sim.FWXmax_rl(i+1), sim.FWXmax_rr(i+1), sim.FWXmax_f(i+1), sim.FWXmax_r(i+1)] = calculateLongiTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i), sim.alpha_r(i));

                % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)      
                [sim.FWYmax_fl(i+1), sim.FWYmax_fr(i+1), sim.FWYmax_rl(i+1), sim.FWYmax_rr(i+1), sim.FWYmax_f(i+1), sim.FWYmax_r(i+1)] = calculateLatTireforces(sim.FWZ_fl(i), sim.FWZ_fr(i), sim.FWZ_rl(i), sim.FWZ_rr(i), GAMMA, TIRparam, sim.alpha_f(i), sim.alpha_r(i));
                
                % Maximum cornering velocity 
                sim.vVYmax(i+1) = sqrt((sim.FWYmax_f(i+1)+sim.FWYmax_r(i+1))*R(i+1)/setup.m_ges); % Calculating maximum possible lateral velocity with given Tire forces [m/s] (inaccuaracy because tire force is based on aero force)

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
            sim.E_Accu = sim.E_Accu/(3.6e6);                % [kWh] Energy consumed by battery per lap (Verbrauchte Akku-Energie je Runde)
            sim.E_heat = sim.E_heat/(3.6e6);                % [kWh] 3.6e6 Joule conversion (3.6e6 Umrechnung Joule)
            sim.E_Accu_Recu = sim.E_Accu_Recu/(3.6e6);      % [kWh] Energy recuperated by batter per lap (Rekuperierte Akku-Energie je Runde)
            sim.E_res = sim.E_Accu - sim.E_Accu_Recu;           % [kWh] Resulting energy consumption per lap (Resultierender Verbrauch je Runde)

            %%  Output of the values (Ausgabe der Werte)
            sim.tEnd = toc;
            
            writeToLogfile(['Simulated with brakes! ' num2str(sim.t(end)) ' s'], startingParameters.Debug, startingParameters.textAreaHandle);         

            sim.t_tot = sim.t(end);

            %% Writing the results to Mat File (Schreiben der Ergebnisse in Mat File)              
            result.Track = Track;
            result.lapLength = lapLength;
            
            result = catstruct(result, setup, startingParameters, sim);                              % Combine Result and Setup to an single struct
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