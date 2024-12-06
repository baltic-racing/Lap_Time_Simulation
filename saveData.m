%% SaveData.m
% Saves the Data of the app to an setup .mat file
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function saveData(app)

    try
        % Saves all data to an .mat file
        %% Chassis Parameters
        setup.m_ges = app.m_ges;
        setup.h_cog = app.h_cog;
        setup.x_cog = app.x_cog;
        setup.x_va = app.x_va;
        setup.m_driver = app.m_driver;
        setup.h_cog_driver = app.h_cog_driver;
        setup.x_cog_driver = app.x_cog_driver;
        setup.m_ballast = app.m_ballast;
        setup.h_cog_ballast = app.h_cog_ballast;
        setup.x_cog_ballast = app.x_cog_ballast;
        setup.thetaV_X = app.thetaV_X;
        setup.thetaV_Y = app.thetaV_Y;
        setup.thetaV_Z = app.thetaV_Z;      

        %% Suspension Parameters
        setup.wheelbase = app.wheelbase;
        setup.track_f = app.track_f;
        setup.track_r = app.track_r;
        setup.J_Tire = app.J_Tire;
        setup.p_Tire = app.p_Tire;
        setup.LMUX = app.LMUX;
        setup.LMUY = app.LMUY;
        setup.LXAL = app.LXAL;
        setup.LYKA = app.LYKA;
        setup.RBX3 = app.RBX3;
        setup.k_R = app.k_R;
        setup.FB = app.FB;
        setup.camber = app.camber;
        setup.m_ph = app.m_ph;
        setup.brakeBias_setup = app.brakeBias;
        setup.h_rc_f = app.h_rc_f;
        setup.h_rc_r = app.h_rc_r;
        setup.tirFile = convertCharsToStrings(app.tirFile);
        setup.KinematicPointsFront = app.KinematicPointsFront;
        setup.KinematicPointsRear = app.KinematicPointsRear;
        %setup.CHAS_LowForFront = app.CHAS_LowForFront;


        %% Drivetrain Parameters
        setup.ptype = app.ptype;
        setup.p_max = app.p_max;
        setup.n_max = app.n_max;
        setup.drivetrain_eff = app.drivetrain_eff;
        setup.invertor_eff = app.invertor_eff; 
        setup.z_chaindrive = app.z_chaindrive;
        setup.z_sprocket = app.z_sprocket;
        setup.trq_multiplier = app.trq_multiplier;
        setup.engine_param = app.engine_param;
        setup.num_motors = app.num_motors;  
        setup.gearbox = app.gearbox;
        setup.i_P = app.i_P;

        if app.gearbox
            setup.i_param = app.i_param;
            setup.n_shift = app.n_shift;
            setup.n_downshift = app.n_downshift;
            setup.t_shift = app.t_shift;
        else
            setup.i_param = [];
            setup.n_shift = [];
            setup.n_downshift = [];
            setup.t_shift = [];
        end
        
        setup.idleRPM = app.idleRPM;

        %% Aerodynamic Parameters
        setup.c_w = app.c_w;
        setup.c_l = app.c_l;
        setup.A = app.A;
        setup.downforce_data = app.downforce_data;
        setup.downforce_multiplier = app.downforce_multiplier;
        setup.aero_pv = app.aero_pv;
        setup.DRS = app.DRS;
        setup.c_d_DRS = app.c_d_DRS;
        setup.c_l_DRS = app.c_l_DRS;
        setup.ConstantDownforce = app.ConstantDownforce;
        setup.DRS_Radius = app.DRS_Radius;

        %% Accumulator Parameters
        setup.V_i = app.V_i;
        setup.Energy_i = app.Energy_i;
        setup.nZellen_Parallel = app.nZellen_Parallel;
        setup.nZellen_Reihe = app.nZellen_Reihe;
        setup.capacity_cell = app.capacity_cell;             

        %% Environment Parameters
        setup.t_L = app.t_L;
        setup.p_L = app.p_L;
        setup.R_L = app.R_L;
        setup.g = app.g;

        % Opens a save file dialog
        [filename, pathname] = uiputfile('*.mat','Save Car Variables As','CarSetup.mat');

        newfilename = fullfile(pathname, filename);

        % Saves file with the given parameters
        save(newfilename, '-struct','setup');
    catch error
        % Write Error Message to log file
        writeToLogfile(error.message);
    end
end