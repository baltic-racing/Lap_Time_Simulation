%% LoadData.m
% Loads the Data of the app from an setup .mat file
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function loadSetup(app, file)

    try
        % Loads all the data from the .mat file
        setup = load(file, '-mat'); 

        %% Chassis Parameters
        app.m_ges = setup.m_ges;
        app.h_cog = setup.h_cog;
        app.x_cog = setup.x_cog;
        app.x_va = setup.x_va;
        app.m_driver = setup.m_driver;
        app.h_cog_driver = setup.h_cog_driver;
        app.x_cog_driver = setup.x_cog_driver;
        app.m_ballast = setup.m_ballast;
        app.h_cog_ballast = setup.h_cog_ballast;
        app.x_cog_ballast = setup.x_cog_ballast;
        app.thetaV_X = setup.thetaV_X;
        app.thetaV_Y = setup.thetaV_Y;
        app.thetaV_Z = setup.thetaV_Z;      

        %% Suspension Parameters
        app.wheelbase = setup.wheelbase;
        app.track = setup.track;
        app.J_Tire = setup.J_Tire;
        app.p_Tire = setup.p_Tire;
        app.LMUX = setup.LMUX;
        app.LMUY = setup.LMUY;
        app.k_R = setup.k_R;
        app.FB = setup.FB;
        app.camber = setup.camber;
        app.m_ph = setup.m_ph;
        app.brakeBias = setup.brakeBias_setup;

        %% Drivetrain Parameters
        app.ptype = setup.ptype;
        app.p_max = setup.p_max;
        app.n_max = setup.n_max;
        app.drivetrain_eff = setup.drivetrain_eff;
        app.invertor_eff = setup.invertor_eff;
        app.z_chaindrive = setup.z_chaindrive;
        app.z_sprocket = setup.z_sprocket;
        app.trq_multiplier = setup.trq_multiplier;
        app.engine_param = setup.engine_param;
        app.num_motors = setup.num_motors;
        app.gearbox = setup.gearbox;
        app.i_P = setup.i_P;

        if app.gearbox
            setup.i_param = app.i_param;
            setup.n_shift = app.n_shift;
            setup.n_downshift = app.n_downshift;
            setup.t_shift = app.t_shift;
        else
            setup.i_param = [0 0];
            setup.n_shift = 0;
            setup.n_downshift = 0;
            setup.t_shift = 0;
        end

        app.i_param = setup.i_param;
        app.n_shift = setup.n_shift;
        app.t_shift = setup.t_shift;
        app.idleRPM = setup.idleRPM;

        %% Aerodynamic Parameters
        app.c_w = setup.c_w;
        app.c_l = setup.c_l;
        app.A = setup.A;
        app.downforce_data = setup.downforce_data;
        app.downforce_multiplier = setup.downforce_multiplier;
        app.aero_pv = setup.aero_pv;
        app.DRS = setup.DRS;
        app.c_d_DRS = setup.c_d_DRS;
        app.c_l_DRS = setup.c_l_DRS;
        app.ConstantDownforce = setup.ConstantDownforce;
        app.DRS_Radius = setup.DRS_Radius;

        %% Accumulator Parameters
        app.V_i = setup.V_i;
        app.Energy_i = setup.Energy_i;
        app.nZellen_Parallel = setup.nZellen_Parallel;
        app.nZellen_Reihe = setup.nZellen_Reihe;
        app.capacity_cell = setup.capacity_cell;  

        %% Environment Parameters
        app.t_L = setup.t_L;
        app.p_L = setup.p_L;
        app.R_L = setup.R_L;
        app.g = setup.g;
    catch error
        % Write Error Message to log file
        writeToLogfile(error.message);
    end
end