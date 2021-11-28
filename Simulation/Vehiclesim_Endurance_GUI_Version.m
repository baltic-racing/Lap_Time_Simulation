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

    Accumulator = [];

    %% Initialisation of all variable
    load('Emraxefficiencydata_interpolated.mat');
    load('CorrectedDischargeInterpolated.mat');
    load('RandomizedCellData.mat');
    load('CellparametersVoltageInterpolation.mat'); 
    
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

    setup.lf = setup.x_cog-setup.x_va;                                  % [mm] Distance from front axle to CoG  
    setup.lr = setup.wheelbase-setup.lf;                                % [mm] Distance from rear axle to CoG

    %% Initalise DRS   
    for i = 1:length(R)
        if abs(R(i)) > setup.DRS_Radius && setup.DRS
            sim.DRS_status(i) = 1;
        else
            sim.DRS_status(i) = 0;
        end
    end     
    
    %% Environmental Conditions (Umgebungsbedingungen)
    sim.rho_L = setup.p_L/(setup.R_L*(setup.t_L+273.15));    % [kg/m^2] Air density (p_L in bar)

    %% Motor Data (Motordaten)
    sim.n = setup.engine_param(:,1);
    sim.M = setup.engine_param(:,2);
    
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
    setup.capacity_accumulator = setup.nZellen_Reihe*setup.nZellen_Parallel*setup.capacity_cell; % [Wh] Capacity of total battery
    
    writeToLogfile('loaded Car Data!', startingParameters.Debug, startingParameters.textAreaHandle);

    %% Tire Model - Magic Tire Formula Model 5.2 (Reifenmodell - Magic Tire Formula Model 5.2)
    
    load('VerticalStiffness_65kPA_IA0.mat');    % Vertical tire stiffness Lookup-Table

    % Load tir-file in structure (Tir-File in Struktur laden)
    TIRparam = loadTIR('C19_CONTINENTAL_FORMULASTUDENT_205_470_R13_65kPa.tir');

    % Tire Data (Reifendaten)
    R0 = TIRparam.UNLOADED_RADIUS;      % [m] Unloaded Tire radius (no wheel loads)
    bW = TIRparam.WIDTH;                % [m] Tire width (contact area) (Breite des Reifens (Aufstandsbreite))
    p_infl = setup.p_Tire;              % [Pa] Tire pressure (Luftdruck des Reifens)
    GAMMA = setup.camber;               % [°] Camber (Sturz)

    % Scaling factors for tire grip level
    % 0.75 for optimum tire temperature; 0.6 for low tire temperature
    TIRparam.LMUX = setup.LMUX;          % [-] Longitudinal scaling factor
    TIRparam.LMUY = setup.LMUY;          % [-] Lateral scaling factor

    % Additional factors for longitudinal/lateral slip interaction
    TIRparam.LXAL = setup.LXAL;         % [-] Influence of slip angle on transmissible longitudinal force
    TIRparam.LYKA = setup.LYKA;         % [-] Influence of longitudinal slip on transmissible lateral force
    TIRparam.RBX3 = setup.RBX3;         % [-] Additional factor Fx for combined slip

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

    % Adjust Variables used for sensitivity Analysis
    if startingParameters.sensitivityID ~= "0"
        setup.(startingParameters.sensitivityID) = startingParameters.minValue;
    end

    if startingParameters.sensitivityID ~= "0"
        setup.(startingParameters.sensitivityID2) = startingParameters.minValue2;  
    end
            
    % Initalize simulation start values
    sim = initializeStartValues(sim, setup.FB, Track, sim.ApexIndexes);

    sim.FB = setup.FB * ones(1,length(Track)-1);

    sim.vAPEXmax = zeros(1,length(sim.ApexIndexes));

    % Axle and wheel loads (static) ((Statische) Achs- und Radlasten)
    sim.FWZtot(1) = sim.FG;                 % [N] Static total axle load (Statische Gesamtachslast)

    sim.m_ph = (setup.x_cog - setup.x_va) / setup.wheelbase;
    sim.FWZf_stat = (1 - (setup.x_cog - setup.x_va) / setup.wheelbase) * sim.FG;
    sim.FWZr_stat = sim.FG - sim.FWZf_stat;    % [N] Static front axle load (Statische Achslast vorne)

    % Static tire loads for calculation
    sim.FWZ_fl_stat = sim.FWZf_stat / 2; 
    sim.FWZ_fr_stat = sim.FWZf_stat / 2;
    sim.FWZ_rl_stat = sim.FWZr_stat / 2; 
    sim.FWZ_rr_stat = sim.FWZr_stat / 2;  

    sim.FWZ_fr(1) = sim.FWZ_fl_stat;         % [N] Static front right wheel load (Statische Radlast vorne rechts)  
    sim.FWZ_fl(1) = sim.FWZ_fr_stat;          % [N] Static front left wheel load (Statische Radlast vorne links)
    sim.FWZ_rr(1) = sim.FWZ_rl_stat;          % [N] Static rear right wheel load (Statische Radlast hinten rechts)
    sim.FWZ_rl(1) = sim.FWZ_rr_stat;          % [N] Static rear left wheel load (Statische Radlast hinten links)

    % Interpolated vertical tire stiffness (Vertikale Reifensteifigkeiten interpoliert) 
    [sim.cZ_fl(1), sim.cZ_fr(1), sim.cZ_rl(1), sim.cZ_rr(1)] = calculateVtirestiff(Fz, cZ_tire, sim.FWZ_fl(1), sim.FWZ_fr(1), sim.FWZ_rl(1), sim.FWZ_rr(1));

    % Dynamic tire radii (stationary = static) (Dynamische Reifenradien (im Stand = statisch))
    [sim.Rdyn_fl(1), sim.Rdyn_fr(1), sim.Rdyn_rl(1), sim.Rdyn_rr(1)] = calculateDynRadii(R0, sim.FWZ_fl(1), sim.FWZ_fr(1), sim.FWZ_rl(1), sim.FWZ_rr(1), sim.cZ_fl(1), sim.cZ_fr(1), sim.cZ_rr(1), sim.cZ_rl(1));

    sim.FWXmax_r(1) = Inf;

    % Calculate Skidpad Time and Speed
    [sim.t_skidpad, sim.vV_skidpad] = calculateSkidPad(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, setup.ConstantDownforce, setup.c_l_DRS, 0, setup.m_ges, setup.lr, setup.lf, setup.wheelbase, setup.track_f, setup.track_r, setup.aero_ph, setup.aero_pv, setup.h_cog, GAMMA, TIRparam, sim.FWZ_fl_stat, sim.FWZ_fr_stat, sim.FWZ_rl_stat, sim.FWZ_rr_stat);
    
    %% Calculation of the maximum apex speed for all apexes (numerically) (Berechnen der maximalen Kurvengeschwindigkeiten f�r alle Apexes (numerisch))
    for i = 1:length(sim.ApexIndexes)

        sim.FWYf(i) = 0;            % [N] Start/Initial value of front axle lateral force (Startwert Querkraft Vorderachse)
        sim.FWYr(i) = 0;            % [N] Start/Initial value of rear axle lateral force (Startwert Querkraft Hinterachse)
        sim.FWYmax_f(i) = 0.1;      % [N] Start/Initial value of maximum transmissible front axle lateral force (Startwert maximal �bertragbare Querkraft Vorderachse)
        sim.FWYmax_r(i) = 0.1;      % [N] Start/Initial value of maximum transmissible rear axle lateral force (Startwert maximal �bertragbare Querkraft Hinterachse)
        sim.vV(i) = 0;              % [m/s] Start/Initial value of vehicle speed (Startwert Fahrzeuggeschwindigkeit)

        %while  sim.FWYf(i) < sim.FWYmax_f(i) && sim.FWYr(i) < sim.FWYmax_r(i) && sim.vV(i) < 50 % ToDo: replace with max speed in last gear
        while  abs(sim.FWYf(i))+abs(sim.FWYr(i)) < abs(sim.FWYmax_f(i))+abs(sim.FWYmax_r(i)) && abs(sim.vV(i)) < 50

            sim.vV(i) = sim.vV(i) + 0.01;   % [m/s] Increaing vehicle speed (Erhoehen der Fahrzeuggeschwindigkeit)
      
            sim.Faero(i) = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, sim.vV(i), setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status(sim.ApexIndexes(i))); % [N] Aerodynamic force

            sim.FVY(i) = setup.m_ges*sim.vV(i)^2/R(sim.ApexIndexes(i));    % [N] Centrifugal force (Zentrifugalkraft)

            sim.aVY(i) = sim.vV(i)^2/R(sim.ApexIndexes(i));  % [m/s^2] Lateral acceleration (Querbeschleunigung)

            % Lateral forces to be applied on front and rear axle (Aufzubringende Querkr�fte an Vorder- und Hinterachse)
            sim.FWYf(i) = setup.lr/setup.wheelbase*abs(sim.FVY(i));   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
            sim.FWYr(i) = setup.lf/setup.wheelbase*abs(sim.FVY(i));   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

            % Wheel load transfer due to aero forces (Radlastverlagerung in Folge von Aerokr�ften) 
            [sim.dFWZrl_aero(i), sim.dFWZrr_aero(i), sim.dFWZfl_aero(i), sim.dFWZfr_aero(i)] = calculateAeroforceOnWheels(sim.Faero(i), setup.aero_ph, setup.aero_pv);

            % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in L�ngsrichtung = 0 angenommen)
            [sim.dFWZfl_x(i), sim.dFWZfr_x(i), sim.dFWZrl_x(i), sim.dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, 0, setup.wheelbase); % Loads = 0 assumed

            % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
            [sim.dFWZfl_y(i), sim.dFWZfr_y(i), sim.dFWZrl_y(i), sim.dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track_f, setup.track_r, setup.lr, setup.lf, setup.wheelbase, sim.FVY(i));

            % Wheel loads (Radlasten)
            sim.FWZ_fl(i) = sim.FWZ_fl_stat + sim.dFWZfl_aero(i) + sim.dFWZfl_x(i) + sim.dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
            sim.FWZ_fr(i) = sim.FWZ_fr_stat + sim.dFWZfr_aero(i) + sim.dFWZfr_x(i) + sim.dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
            sim.FWZ_rl(i) = sim.FWZ_rl_stat + sim.dFWZrl_aero(i) + sim.dFWZrl_x(i) + sim.dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
            sim.FWZ_rr(i) = sim.FWZ_rr_stat + sim.dFWZrr_aero(i) + sim.dFWZrr_x(i) + sim.dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts) 
            
            % Maximum transmissible tire forces in longitudinal direction = 0 assumed (because longitudinal wheel loads = 0 assumed) 
            
            % Calculate delta, beta, psi1 and alpha for all wheels
            % and front / rear
            [sim.delta(i), sim.beta(i), sim.psi1(i), sim.alpha_f(i), sim.alpha_r(i), sim.alpha_fr(i), sim.alpha_fl(i), sim.alpha_rr(i), sim.alpha_rl(i), sim.delta_fl(i), sim.delta_fr(i), sim.delta_sw(i), sim.ackermann(i), sim.ackermannPercent(i)] = calculateSteeringData(setup.wheelbase, R(sim.ApexIndexes(i)), setup.lr, setup.lf, sim.vV(i), sim.FWZ_fl(i), sim.FWZ_rl(i), sim.FWZ_fr(i), sim.FWZ_rr(i), setup.track_f);         

            % Maximum transmissible tire forces in lateral direction   
            [sim.FWYmax_f(i), sim.FWYmax_r(i)] = calculateLatTireforces(sim.FWZ_fl(i), sim.FWZ_fr(i), sim.FWZ_rl(i), sim.FWZ_rr(i), GAMMA, TIRparam, sim.alpha_f(i), sim.alpha_r(i));

        end

        sim.vAPEXmax(i) = sim.vV(i);   % [m/s] Maximum speed for any apex
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

    %% Simulation WITHOUT BRAKES
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

        % Pointer for efficiency table
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
        [sim.FVX_fl(i), sim.FVX_fr(i), sim.FVX_rl(i), sim.FVX_rr(i), sim.FVX(i), sim.FVX_f(i), sim.TC_front(i), sim.TC(i)] = calculateTractiveForces(sim.Mi(i), setup.num_motors, sim.i_G, sim.gr, sim.Rdyn_fl(i), sim.Rdyn_fr(i), sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.t_x, sim.gear, sim.FWXmax_f(i), sim.FWXmax_r(i), setup.t_shift, sim.FWYmax_f(i), sim.FWYmax_r(i), sim.FWYf(i), sim.FWYr(i));

        % Driving resistances      
        [sim.FR(i), sim.FL(i), sim.Fdr(i), sim.FVY(i), sim.aVX(i), sim.aVY(i)] = calculateVehicleResistancesForces(setup.k_R, sim.FWZtot(i), sim.rho_L, sim.vV(i), setup.c_w, setup.A, setup.m_ges, R(i), sim.FVX(i), sim.FVX_f(i), setup.c_d_DRS, sim.DRS_status(i), sim.rpmpointer, setup.n_max, 0, 0, 0, 0);

        if ismember(i,sim.ApexIndexes)
            if sim.vV(i) > sim.vAPEXmax(sim.z)   % Limiting maximum speeds at apexes (Begrenzen auf maximale Kurvengeschwindigkeit in Apexes)
                sim.vV(i) = sim.vAPEXmax(sim.z);
                sim.ni(i) = sim.vV(i) * 30 / pi * sim.gr(sim.gear) * sim.i_G / ((sim.Rdyn_rr(i)+sim.Rdyn_rl(i))/2);
            end
            sim.z = sim.z + 1;
        end
        
        if sim.ni(i) == setup.n_max && sim.aVX(i) == 0
            sim.vV(i+1) = (setup.n_max * pi * ((sim.Rdyn_rr(i)+sim.Rdyn_rl(i))/2)) / (30 * sim.gr(sim.gear) * sim.i_G);
            %sim.vV(i+1) = (setup.n_max * pi * R0) / (30 * sim.gr(sim.gear) * sim.i_G);
        else
            sim.vV(i+1) = sqrt(sim.vV(i)^2+2*sim.aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
        end
        
        sim.t(i+1) = sim.t(i)+(s(i+1)-s(i))/sim.vV(i+1);            % [s] Time (Zeit)

        % Lateral forces to be applied on front and rear axles (Aufzubringende Querkr�fte an Vorder- und Hinterachse)
        sim.FWYf(i) = setup.lr/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied on front axle (Aufzubringende Querkraft der Vorderachse)
        sim.FWYr(i) = setup.lf/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied on rear axle (Aufzubringende Querkraft der Hinterachse)                                  

        % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokr�ften)  
        [sim.dFWZrl_aero(i), sim.dFWZrr_aero(i), sim.dFWZfl_aero(i), sim.dFWZfr_aero(i)] = calculateAeroforceOnWheels(sim.Faero(i), setup.aero_ph, setup.aero_pv);

        % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerungen in L�ngsrichtung)
        [sim.dFWZfl_x(i), sim.dFWZfr_x(i), sim.dFWZrl_x(i), sim.dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, sim.aVX(i), setup.wheelbase);

        % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
        [sim.dFWZfl_y(i), sim.dFWZfr_y(i), sim.dFWZrl_y(i), sim.dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track_f, setup.track_r, setup.lr, setup.lf, setup.wheelbase, sim.FVY(i));

        % Wheel loads (Radlasten)
        sim.FWZ_fl(i) = sim.FWZ_fl_stat + sim.dFWZfl_aero(i) + sim.dFWZfl_x(i) + sim.dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
        sim.FWZ_fr(i) = sim.FWZ_fr_stat + sim.dFWZfr_aero(i) + sim.dFWZfr_x(i) + sim.dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
        sim.FWZ_rl(i) = sim.FWZ_rl_stat + sim.dFWZrl_aero(i) + sim.dFWZrl_x(i) + sim.dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
        sim.FWZ_rr(i) = sim.FWZ_rr_stat + sim.dFWZrr_aero(i) + sim.dFWZrr_x(i) + sim.dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)

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
        
        % Calculate delta, beta, psi1 and alpha for all wheels
        % and front / rear
        [sim.delta(i+1), sim.beta(i+1), sim.psi1(i+1), sim.alpha_f(i+1), sim.alpha_r(i+1), sim.alpha_fr(i+1), sim.alpha_fl(i+1), sim.alpha_rr(i+1), sim.alpha_rl(i+1), sim.delta_fl(i+1), sim.delta_fr(i+1), sim.delta_sw(i+1), sim.ackermann(i+1), sim.ackermannPercent(i+1)] = calculateSteeringData(setup.wheelbase, R(i+1), setup.lr, setup.lf, sim.vV(i+1), sim.FWZ_fl(i+1), sim.FWZ_rl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rr(i+1), setup.track_f);  

        % Vertical tire stiffnesses - for dynamic radii (Vertikale Reifensteifigkeiten)  
        [sim.cZ_fl(i+1), sim.cZ_fr(i+1), sim.cZ_rl(i+1), sim.cZ_rr(i+1)] = calculateVtirestiff(Fz, cZ_tire, sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1));

        % Dynamic tire radii (Dynamische Reifenradien)
        [sim.Rdyn_fl(i+1), sim.Rdyn_fr(i+1), sim.Rdyn_rl(i+1), sim.Rdyn_rr(i+1)] = calculateDynRadii(R0, sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), sim.cZ_fl(i+1), sim.cZ_fr(i+1), sim.cZ_rr(i+1), sim.cZ_rl(i+1));

        % Maximum transmissible tire forces in longitudinal direction (Maximal �bertragbare Reifenkr�fte in L�ngsrichtung)
        [sim.FWXmax_fl(i+1), sim.FWXmax_fr(i+1), sim.FWXmax_rl(i+1), sim.FWXmax_rr(i+1), sim.FWXmax_f(i+1), sim.FWXmax_r(i+1)] = calculateLongiTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i+1), sim.alpha_r(i+1));

        % Maximum transmissible tire forces in lateral direction (Maximal �bertragbare Reifenkr�fte in Querrichtung)
        [sim.FWYmax_fl(i+1), sim.FWYmax_fr(i+1), sim.FWYmax_rl(i+1), sim.FWYmax_rr(i+1), sim.FWYmax_f(i+1), sim.FWYmax_r(i+1)] = calculateLatTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i+1), sim.alpha_r(i+1));
        
        % Maximum cornering velocity 
        sim.vVYmax(i+1) = sqrt((abs(sim.FWYmax_f(i+1))+abs(sim.FWYmax_r(i+1)))*abs(R(i+1))/setup.m_ges); % Calculating maximum possible lateral velocity with given Tire forces [m/s] (inaccuaracy because tire force is based on aero force)
        
%         if (sim.vV(i+1) > sim.vVYmax(i+1))
%             sim.vV(i+1) = sim.vVYmax(i+1);
%         end
    end
    
    sim.vWoBrake = sim.vV;                                      % Save velocity without braking for log file.
    sim.avXWoBrake = sim.aVX;                                   % Save longitudinal acceleration for log file.
    sim.avYWoBrake = sim.aVY;                                   % Save lateral acceleration for log file.

    writeToLogfile('Simulated without brakes!', startingParameters.Debug, startingParameters.textAreaHandle);

    %% BRAKING POINT CALCULATION (BREMSPUNKTBERECHNUNG)
    [sim.BrakeIndexes, sim.NonBrakeApexes, sim.vRev] = calculateBrakepoints(setup.FB, Track, sim.ApexIndexes, sim.vAPEXmax, setup.m_ges, setup.downforce_multiplier, setup.c_l, setup.c_w, setup.A, sim.rho_L, setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status, setup.aero_ph, setup.aero_pv, sim.vV, setup.k_R, sim.FG, setup.h_cog, setup.wheelbase, setup.track_f, setup.track_r, setup.lr, setup.lf, GAMMA, TIRparam, sim.FWZ_fl_stat, sim.FWZ_fr_stat, sim.FWZ_rl_stat, sim.FWZ_rr_stat, R, s, setup.brakeBias_setup, startingParameters.brakeFunction);

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

        % Checking at which apex the vehicle is
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
   
        %% Braking
        % Checking if braking is required
        if ismember(i,sim.BrakeIndexes) && not(ismember(sim.z,sim.NonBrakeApexes))                      % Initiaion of braking process (Einleiten des Bremsvorgangs)     
            
            [sim.Faero(i)] = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, sim.vV(i), setup.ConstantDownforce, setup.c_l_DRS, 0); % [N] Aerodynamic force
            
            sim.DRS_status(i) = 0;
            
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

            sim.FVX_fl(i) = 0;                  % [N] Tractive Force on front left wheel (AWD) 
            sim.FVX_fr(i) = 0;                  % [N] Tractive Force on front right wheel (AWD) 
            sim.FVX_f(i) = 0;                   % [N] Traction on rear axle (Zugkraft an der Hinterachse)

            sim.FVX_rl(i) = 0;                  % [N] Traction on rear left wheel (Zugkraft an linkem Hinterrad)
            sim.FVX_rr(i) = 0;                  % [N] Traction on rear right wheel (Zugkraft an rechtem Hinterrad)
            sim.FVX(i) = 0;                     % [N] Traction on rear axle (Zugkraft an der Hinterachse)              
            
            [sim.FB_fl(i), sim.FB_fr(i), sim.FB_rl(i), sim.FB_rr(i), ~, sim.BrakeBias(i), sim.ABS(i)] = calculateDeceleration(sim.FB(i), setup.m_ges, sim.Fdr(i), sim.FWXmax_fl(i), sim.FWXmax_fr(i), sim.FWXmax_rl(i), sim.FWXmax_rr(i), setup.brakeBias_setup);
        else
        %% Accelerating 
        
            [sim.Faero(i)] = calculateAeroforce(setup.downforce_multiplier, setup.c_l, setup.A, sim.rho_L, sim.vV(i), setup.ConstantDownforce, setup.c_l_DRS, sim.DRS_status(i)); % [N] Aerodynamic force
        
            sim.Mi(i) = interp1(sim.n,sim.M,sim.ni(i),'linear','extrap');   % [Nm] Motor torque (single motor!)
            sim.FB(i) = 0;                                                  % [N] Braking force
            sim.P_Bh(i) = 0;                                                % [W] Rear braking power (Bremsleistung hinten)
            sim.FB_fl(i) = 0;                                               % [N] Brake Force front left
            sim.FB_fr(i) = 0;                                               % [N] Brake Force front right
            sim.FB_rl(i) = 0;                                               % [N] Brake Force rear left
            sim.FB_rr(i) = 0;                                               % [N] Brake Force rear right

            % Motor power & limitation to 80 kW from FS-Rules (Motorleistung & Begrenzung auf 80 kW aus FS-Rules)
            sim.rpmpointer = round(sim.ni(i));                              % Pointer for efficiency table 

            if sim.Mi(i) <= 0
                sim.torquepointer = 1;
            else
                sim.torquepointer = round(sim.Mi(i));                       % Pointer for efficiency table (Pointer f�r effizienztabelle)
            end 

            sim.P_M(i) = setup.num_motors * sim.Mi(i) * sim.ni(i) / 60 * 2 * pi; % [W] Total motor power (P_el!)

            % Limiting the maximal power when using an electric
            % drivetrain
            if setup.ptype && sim.P_M(i) > setup.p_max
                sim.P_M(i) = setup.p_max;                                       % [W] Limited power (P_el!)
                %Mi(i) = P_M(i)*60/ni(i)/2/pi;                  % [Nm] Limiting the torque (Total Motor Torque!)
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
            [sim.FVX_fl(i), sim.FVX_fr(i), sim.FVX_rl(i), sim.FVX_rr(i), sim.FVX(i), sim.FVX_f(i), sim.TC_front(i), sim.TC(i)] = calculateTractiveForces(sim.Mi(i), setup.num_motors, sim.i_G, sim.gr, sim.Rdyn_fl(i), sim.Rdyn_fr(i), sim.Rdyn_rl(i), sim.Rdyn_rr(i), sim.t_x, sim.gear, sim.FWXmax_f(i), sim.FWXmax_r(i), setup.t_shift, sim.FWYmax_f(i), sim.FWYmax_r(i), sim.FWYf(i), sim.FWYr(i));

            sim.M_tractive(i) = (sim.FVX(i)+sim.FVX_f(i))/(sim.i_G*sim.gr(sim.gear)/sim.Rdyn_rr(i));        % [Nm] Torque including tractive force
            sim.P_tractive(i) = sim.M_tractive(i)/(60/sim.ni(i)/2/pi);                                      % [kW] Motor power required for tractive forces
            sim.P_el(i) = (sim.P_tractive(i)/(setup.drivetrain_eff * sim.motor_eff(i) * sim.eta_inv));      % [kW] Motor power required for tractive forces including efficiencies
        end

        % Driving resistances (Fahrwiderstaende) & Vehicle (Fahrzeug)
        [sim.FR(i), sim.FL(i), sim.Fdr(i), sim.FVY(i), sim.aVX(i), sim.aVY(i)] = calculateVehicleResistancesForces(setup.k_R, sim.FWZtot(i), sim.rho_L, sim.vV(i), setup.c_w, setup.A, setup.m_ges, R(i), sim.FVX(i), sim.FVX_f(i), setup.c_d_DRS, sim.DRS_status(i), sim.rpmpointer, setup.n_max, sim.FB_fl(i), sim.FB_fr(i), sim.FB_rl(i), sim.FB_rr(i));    
        
        % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
        if sim.ni(i) == setup.n_max && sim.aVX(i) == 0
            sim.vV(i+1) = (setup.n_max * pi * ((sim.Rdyn_rr(i)+sim.Rdyn_rl(i))/2)) / (30 * sim.gr(sim.gear) * sim.i_G);
            %sim.vV(i+1) = (setup.n_max * pi * R0) / (30 * sim.gr(sim.gear) * sim.i_G);
        else
            sim.vV(i+1) = sqrt(sim.vV(i)^2+2*sim.aVX(i)*(s(i+1)-s(i))); % [m/s] Total vehicle speed (Gesamt-Fahrzeuggeschwindigkeit)
        end

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
        
        % Latschlaengen fuer Laengsschlupfberechnung nach Carter
        % (WRONG USE FVX for each wheel wrong formula for front wheels!)
        sim.FU_fl(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_fl(i+1)-sim.FB_fl(i);   % [N] Umfangskraefte an einem Hinterrad
        sim.FU_fr(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_fr(i+1)-sim.FB_fr(i);   % [N] Umfangskraefte an einem Hinterrad
        sim.FU_rl(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_rl(i+1)-sim.FB_rl(i);   % [N] Umfangskraefte an einem Hinterrad
        sim.FU_rr(i) = sim.FVX(i)/2-setup.k_R*sim.FWZ_rr(i+1)-sim.FB_rr(i);   % [N] Umfangskraefte an einem Hinterrad
        
        sim.l_contact_patch_fl(i) = sim.FWZ_fl(i)/(p_infl*bW);   % [m] Latschlaenge vorne links
        sim.l_contact_patch_fr(i) = sim.FWZ_fr(i)/(p_infl*bW);   % [m] Latschlaenge vorne rechts
        sim.l_contact_patch_rl(i) = sim.FWZ_rl(i)/(p_infl*bW);   % [m] Latschlaenge hinten links
        sim.l_contact_patch_rr(i) = sim.FWZ_rr(i)/(p_infl*bW);   % [m] Latschlaenge hinten rechts

        % Laengsschlupf nach Carter
        sim.my0 = 2;    % [-] Haftreibungsbeiwert   
        sim.kappa_fl(i) = sim.l_contact_patch_fl(i)/(2*sim.Rdyn_fl(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_fl(i)/(sim.my0*sim.FWZ_fl(i))));
        sim.kappa_fr(i) = sim.l_contact_patch_fr(i)/(2*sim.Rdyn_fr(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_fr(i)/(sim.my0*sim.FWZ_fr(i))));
        sim.kappa_rl(i) = sim.l_contact_patch_rl(i)/(2*sim.Rdyn_rl(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_rl(i)/(sim.my0*sim.FWZ_rl(i))));
        sim.kappa_rr(i) = sim.l_contact_patch_rr(i)/(2*sim.Rdyn_rr(i))*sim.my0*(setup.wheelbase-sqrt(1-sim.FU_rr(i)/(sim.my0*sim.FWZ_rr(i))));

        sim.E_Accu(i+1) = sim.E_Accu(i) + (sim.t(i+1)-sim.t(i)) * sim.P_el(i); % [J] 
        sim.E_heat(i+1) = sim.E_heat(i) + (sim.P_el(i)-sim.P_tractive(i)) * (sim.t(i+1)-sim.t(i));
        %E_Waerme(i+1) = P_el(i)*(1-M_eff_inter(i))+(E_Akku(i) + (t(i+1)-t(i)))*(1-eta_inv);    % [J] Motor losses + 5% for inverter losses; Drivetrain losses with heat (Motorverluste + 5% flat fuer inverterverlsute % Drivetrain Verluste auch Waerme)


        % Lateral forces on front and rear axle (Querkraefte an Vorder- und Hinterachse)
        sim.FWYf(i) = setup.lr/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied to the front axle (Aufzubringende Querkraft der Vorderachse)
        sim.FWYr(i) = setup.lf/setup.wheelbase*sim.FVY(i);   % [N] Lateral force to be applied to the rear axle (Aufzubringende Querkraft der Hinterachse)

        if sim.FWYf(i) > sim.FWYmax_f(i)
           sim.slipY_f(i) = 1;  
        end

         if sim.FWYr(i) > sim.FWYmax_r(i)
           sim.slipY_r(i) = 1;  
        end

        % Wheel load transfer due to aerodynamic forces (Radlastverlagerung in Folge von Aerokraeften)        
        [sim.dFWZrl_aero(i), sim.dFWZrr_aero(i), sim.dFWZfl_aero(i), sim.dFWZfr_aero(i)] = calculateAeroforceOnWheels(sim.Faero(i), setup.aero_ph, setup.aero_pv);

        % Dynamic wheel load displacements in longitudinal direction (Dynamische Radlastverlagerungen in Laengsrichtung)       
        [sim.dFWZfl_x(i), sim.dFWZfr_x(i), sim.dFWZrl_x(i), sim.dFWZrr_x(i)] = calculateWheelloadLongDisp(setup.h_cog, setup.m_ges, sim.aVX(i), setup.wheelbase);

        % Dynamic wheel load displacements in lateral direction (Dynamische Radlastverlagerung in Querrichtung)  
        [sim.dFWZfl_y(i), sim.dFWZfr_y(i), sim.dFWZrl_y(i), sim.dFWZrr_y(i)] = calculateWheelloadLatDisp(setup.h_cog, setup.track_f, setup.track_r, setup.lr, setup.lf, setup.wheelbase, sim.FVY(i));

        % Wheel loads (Radlasten)
        sim.FWZ_fl(i) = sim.FWZ_fl_stat + sim.dFWZfl_aero(i) + sim.dFWZfl_x(i) + sim.dFWZfl_y(i); % [N] Front left wheel load (Radlast vorne links)
        sim.FWZ_fr(i) = sim.FWZ_fr_stat + sim.dFWZfr_aero(i) + sim.dFWZfr_x(i) + sim.dFWZfr_y(i); % [N] Front right wheel load (Radlast vorne rechts)
        sim.FWZ_rl(i) = sim.FWZ_rl_stat + sim.dFWZrl_aero(i) + sim.dFWZrl_x(i) + sim.dFWZrl_y(i); % [N] Rear left wheel load (Radlast hinten links)
        sim.FWZ_rr(i) = sim.FWZ_rr_stat + sim.dFWZrr_aero(i) + sim.dFWZrr_x(i) + sim.dFWZrr_y(i); % [N] Rear right wheel load (Radlast hinten rechts)

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

        % Calculate delta, beta, psi1 and alpha for all wheels
        % and front / rear
        [sim.delta(i+1), sim.beta(i+1), sim.psi1(i+1), sim.alpha_f(i+1), sim.alpha_r(i+1), sim.alpha_fr(i+1), sim.alpha_fl(i+1), sim.alpha_rr(i+1), sim.alpha_rl(i+1), sim.delta_fl(i+1), sim.delta_fr(i+1), sim.delta_sw(i+1), sim.ackermann(i+1), sim.ackermannPercent(i+1)] = calculateSteeringData(setup.wheelbase, R(i+1), setup.lr, setup.lf, sim.vV(i+1), sim.FWZ_fl(i+1), sim.FWZ_rl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rr(i+1), setup.track_f);  
        
        % Maximum transmissible tire forces in longitudinal direction   
        [sim.FWXmax_fl(i+1), sim.FWXmax_fr(i+1), sim.FWXmax_rl(i+1), sim.FWXmax_rr(i+1), sim.FWXmax_f(i+1), sim.FWXmax_r(i+1)] = calculateLongiTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i+1), sim.alpha_r(i+1));

        % Maximum transmissible tire forces in lateral direction
        [sim.FWYmax_fl(i+1), sim.FWYmax_fr(i+1), sim.FWYmax_rl(i+1), sim.FWYmax_rr(i+1), sim.FWYmax_f(i+1), sim.FWYmax_r(i+1)] = calculateLatTireforces(sim.FWZ_fl(i+1), sim.FWZ_fr(i+1), sim.FWZ_rl(i+1), sim.FWZ_rr(i+1), GAMMA, TIRparam, sim.alpha_f(i+1), sim.alpha_r(i+1));
        
        % Maximum cornering velocity
        sim.vVYmax(i+1) = sqrt((abs(sim.FWYmax_f(i+1))+abs(sim.FWYmax_r(i+1)))*abs(R(i+1))/setup.m_ges); % Calculating maximum possible lateral velocity with given Tire forces [m/s] (inaccuaracy because tire force is based on aero force)
        
%         if (sim.vV(i+1) > sim.vVYmax(i+1))
%             sim.vV(i+1) = sim.vVYmax(i+1);
%         end

        [sim.rollMoment_f(i), sim.rollMoment_r(i), sim.rollAngleChassis(i)] = calculateRollMoment(sim.aVY(i), setup.m_ges, setup.h_cog, setup.m_ph, setup.h_rc_f, setup.h_rc_r, setup.g, setup.x_va, setup.x_cog, setup.wheelbase);

        sim.V_i(i) = sum(sim.Voltage_Cellpack(:,i));
        
        % Battery Currents (Akkustroeme)
        sim.A_accu_cell(i) = sim.P_el(i) / sim.V_i(i) / setup.nZellen_Parallel;  
        
        sim.Current_Cellpack_Pointer(i) = sim.P_M(i) / sim.V_i(i) * 10 ; %Strombelastung eines %er Parrallel Paketes in 0,1A parameter fuer die berechnung der korrigierten belastung mit hoehren verlusten durch hoehere zellstroeme
        if sim.Current_Cellpack_Pointer(i) <= 1
            sim.Current_Cellpack_Pointer(i)=1;
        end
        
        if sim.Current_Cellpack_Pointer(i) >= 1500 %begrenzen des max Zellstromes auf 30A pro Zelle im 5er parralelverbund also 150A
            sim.Current_Cellpack_Pointer(i)=1500;
        end
        
        sim.VirtualCurrent_Cellpack(i) = CorrectedDischargeInterpolated(1,round(sim.Current_Cellpack_Pointer(i))); %Berechnung der Virtuell hoeheren zellstr�me basierend auf den h�heren verlsuten durch h�here Str�me
        
        sim.Energy_Cellpack(i) = (sim.VirtualCurrent_Cellpack(i)*(sim.t(i+1)-sim.t(i))) - ((sim.P_Bh(i) / sim.V_i(i)) * (sim.t(i+1)-sim.t(i))) ; %Energieverbrauch in As f�r ein 5erpacket an akkuzellen -> Akkustrom zum zeitpunkt i mal Zeitdifferenz zwischen i und i+1
        sim.Energy_Cellpack_Total(i+1) = sim.Energy_Cellpack_Total(i) + sim.Energy_Cellpack(i); % �ber Endurance Run Integrierte Energieverbrauch in As f�r ein 5erpacket an akkuzellen
        
        sim.Capacity_Cellpack(1:131,i+1) =  sim.Capacity_Cellpack(1:131,i) - sim.Energy_Cellpack(i); 
        
        sim.SOC_Cellpack(1:131,i+1) = sim.Capacity_Cellpack(1:131,i) ./ sim.Capacity_Cellpack(1:131,1); %Berechnung des SOC f�r den n�chsten tick basierend auf der aktuellen cellcapacity und der im n�chsten tick
        
        sim.SOC_Pointer(1:131,i+1) = round(sim.SOC_Cellpack(1:131,i+1)*1000);
        sim.Current_Cellpack_Pointer_Voltage(1,i+1) = round(sim.Current_Cellpack_Pointer(i)/5);
        
        if sim.Current_Cellpack_Pointer_Voltage(i) <= 3
            sim.Current_Cellpack_Pointer_Voltage(i) = 3;
        end
        
        if size(Track,1) < sim.SOC_Pointer(1:131,i+1)
            sim.SOC_Pointer(1:131,i+1) = size(Track,1); 
        end
        
        if size(Track,1) < sim.Current_Cellpack_Pointer_Voltage(1,i)
            sim.Current_Cellpack_Pointer_Voltage(1,i) = size(Track,1);
        end
        
        try % Catch Exception when sim is empty 
            sim.Voltage_Cellpack(1:131,i+1) = Voltage_inter(sim.Current_Cellpack_Pointer_Voltage(1,i),sim.SOC_Pointer(1,i+1));
        catch
        end

        %[Accumulator(i)] = calculateAccumulatorData(Accumulator, sim.P_el(i), sim.P_M(i), sim.P_Bh(i), sim.t, i);  
    end
    
    % Conversion of battery energy capacity (Umrechnen der Energiemengen des Akkus)
    sim.E_Accu = sim.E_Accu / (3.6e6);                % [kWh] Energy consumed by battery per lap (Verbrauchte Akku-Energie je Runde)
    sim.E_heat = sim.E_heat / (3.6e6);                % [kWh] 3.6e6 Joule conversion (3.6e6 Umrechnung Joule)
    sim.E_Accu_Recu = sim.E_Accu_Recu / (3.6e6);      % [kWh] Energy recuperated by batter per lap (Rekuperierte Akku-Energie je Runde)
    sim.E_res = sim.E_Accu - sim.E_Accu_Recu;         % [kWh] Resulting energy consumption per lap (Resultierender Verbrauch je Runde)

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
    if startingParameters.sensitivityID2 ~= "0" && startingParameters.textAreaHandle ~= 0
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
    
    if startingParameters.textAreaHandle ~= 0
        [~, name, ~] = fileparts(startingParameters.carDatafile);
        
        % Saves the results to a .mat file in the same location as the
        % setup file. The name is generated by adding _result to the name.
        savefilename = name + "_result.mat";

        save(startingParameters.path+"/"+savefilename, '-struct','result');
        
        writeToLogfile('File succesfully written', startingParameters.Debug, startingParameters.textAreaHandle);
    end
end
