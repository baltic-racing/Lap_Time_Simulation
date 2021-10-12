function PlotSensitivity(resultFile, xparam, yparam, dots, zparam, resx, resy)
    
    result1 = load(resultFile, '-mat');

    x_axis = {};
    y_axis = {};
    xlabelText = '';
    ylabelText = '';
    
    try 
        % Checks if the track file is 4 or 5 columns and read it acordingly          
        [~, columns] = size(result1(1).Track);
    catch 
        result1 = result1.result;

        % Checks if the track file is 4 or 5 columns and read it acordingly          
        [~, columns] = size(result1(1).Track);
    end
    
    switch xparam  
        case 1
            x_axis = [result1.t_ges];
            xlabelText = 'lap time [s]';
        case 2
            x_axis = [result1.A];
            xlabelText = 'A [m²]';
        case 3
            x_axis = [result1.FB];
            xlabelText = 'Maximum Braking Force [N]';
        case 4
            x_axis = [result1.LMUX];
            xlabelText = 'Longitudinal scaling factor [-]';
        case 5
            x_axis = [result1.LMUY];
            xlabelText = 'Lateral scaling factor [-]';
        case 6
            x_axis = [result1.aero_pv];
            xlabelText = 'Aerodynamic Force on front axle [-]';
        case 7
            x_axis = [result1.c_l];
            xlabelText = 'c_l value [-]';
        case 8
            x_axis = [result1.c_w];
            xlabelText = 'c_d value [-]';
        case 9
            x_axis = [result1.camber];
            xlabelText = 'Camber [°]';
        case 10
            x_axis = [result1.downforce_multiplier];
            xlabelText = 'downforce multiplier [-]';
        case 11
            x_axis = [result1.drivetrain_eff];
            xlabelText = 'Overall efficiency of powertrain [-]';
        case 12
            x_axis = [result1.m_balast];
            xlabelText = 'mass of ballast [kg]';
        case 13
            x_axis = [result1.m_driver];
            xlabelText = 'mass of driver [kg]';
        case 14
            x_axis = [result1.m_tot];
            xlabelText = 'total mass [kg]';
        case 15
            x_axis = [result1.m_ph];
            xlabelText = 'Percentage of rear axle wheel load [%]';
        case 16
            x_axis = [result1.n_max];
            xlabelText = 'Maximum RPM of the motor [1/min]';
        case 17
            x_axis = [result1.num_motors];
            xlabelText = 'number of E-motors [-]';
        case 18
            x_axis = [result1.p_Tire];
            xlabelText = 'Tire pressure [Pa]';
        case 19
            x_axis = [result1.P_max];
            xlabelText = 'Power Limit in endurance [W]';
        case 20
            x_axis = [result1.t_L];
            xlabelText = 'Ambient air temperature [°C]';
        case 21
            x_axis = [result1.thetaV_X];
            xlabelText = 'Moment of Inertia of vehicle about X-Axis [kg*m²]';
        case 22
            x_axis = [result1.thetaV_Y];
            xlabelText = 'Moment of Inertia of vehicle about Y-Axis [kg*m²]';
        case 23
            x_axis = [result1.thetaV_Z];
            xlabelText = 'Moment of Inertia of vehicle about Z-Axis [kg*m²]';
        case 24
            x_axis = [result1.track];
            xlabelText = 'track [-]';
        case 25
            x_axis = [result1.trq_multiplier];
            xlabelText = 'Torque multiplier [-]';
        case 26
            x_axis = [result1.wheelbase];
            xlabelText = 'wheelbase [mm]';
        case 27
            x_axis = [result1.h_COG];
            xlabelText = 'Height of vehicles Center of Gravity(COG)[mm]';
        case 28
            x_axis = [result1.x_cog];
            xlabelText = 'x-Coordinate of vehicles COG [mm]';
        case 29
            x_axis = [result1.x_cog_balast];
            xlabelText = 'x-Coordinate COG [mm]';
        case 30
            x_axis = [result1.x_cog_driver];
            xlabelText = 'x-Coordinate COG [mm]';
        case 31
            x_axis = [result1.x_va];
            xlabelText = 'x-Coordinate of the front axle [mm]';
        case 32
            x_axis = [result1.z_chaindrive];
            xlabelText = 'Number of teeth on sprocket [-]';
        case 33
            x_axis = [result1.z_sprocket];
            xlabelText = 'Number of teeth on pinion [-]';
        case 34
            x_axis = [result1.i_G];
            xlabelText = 'Gear Ratio [-]';
        case 35
            x_axis = [result1.ConstantDownforce];
            xlabelText = 'Constant Downforce [N]';
        case 36
            x_axis = [result1.DRS_Radius];
            xlabelText = 'DRS Radius [m]';
        case 37
            x_axis = [result1.vV_skidpad];
            xlabelText = 'vV_skidpad [m/s]';
        case 38
            x_axis = [result1.t_skidpad];
            xlabelText = 't_skidpad [m/s]';   
        case 39
            x_axis = [result1.c_l]./[result1.c_w];
            xlabelText = 'Aero Efficiency c_l/c_w [-]';
    end
    
    switch yparam  
        case 1
            y_axis = [result1.t_ges];
            ylabelText = 'lap time [s]';
        case 2
            y_axis = [result1.A];
            ylabelText = 'A [m²]';
        case 3
            y_axis = [result1.FB];
            ylabelText = 'Maximum Braking Force [N]';
        case 4
            y_axis = [result1.LMUX];
            ylabelText = 'Longitudinal scaling factor [-]';
        case 5
            y_axis = [result1.LMUY];
            ylabelText = 'Lateral scaling factor [-]';
        case 6
            y_axis = [result1.aero_pv];
            ylabelText = 'Aerodynamic Force on front axle [-]';
        case 7
            y_axis = [result1.c_l];
            ylabelText = 'c_l value [-]';
        case 8
            y_axis = [result1.c_w];
            ylabelText = 'c_d value [-]';
        case 9
            y_axis = [result1.camber];
            ylabelText = 'Camber [°]';
        case 10
            y_axis = [result1.downforce_multiplier];
            ylabelText = 'downforce multiplier [-]';
        case 11
            y_axis = [result1.drivetrain_eff];
            ylabelText = 'Overall efficiency of powertrain [-]';
        case 12
            y_axis = [result1.m_balast];
            ylabelText = 'mass of ballast [kg]';
        case 13
            y_axis = [result1.m_driver];
            ylabelText = 'mass of driver [kg]';
        case 14
            y_axis = [result1.m_tot];
            ylabelText = 'total mass [kg]';
        case 15
            y_axis = [result1.m_ph];
            ylabelText = 'Percentage of rear axle wheel load [%]';
        case 16
            y_axis = [result1.n_max];
            ylabelText = 'Maximum RPM of the motor [1/min]';
        case 17
            y_axis = [result1.num_motors];
            ylabelText = 'number of E-motors [-]';
        case 18
            y_axis = [result1.p_Tire];
            ylabelText = 'Tire pressure [Pa]';
        case 19
            y_axis = [result1.P_max];
            ylabelText = 'Power Limit in endurance [W]';
        case 20
            y_axis = [result1.t_L];
            ylabelText = 'Ambient air temperature [°C]';
        case 21
            y_axis = [result1.thetaV_X];
            ylabelText = 'Moment of Inertia of vehicle about X-Axis [kg*m²]';
        case 22
            y_axis = [result1.thetaV_Y];
            ylabelText = 'Moment of Inertia of vehicle about Y-Axis [kg*m²]';
        case 23
            y_axis = [result1.thetaV_Z];
            ylabelText = 'Moment of Inertia of vehicle about Z-Axis [kg*m²]';
        case 24
            y_axis = [result1.track];
            ylabelText = 'track [-]';
        case 25
            y_axis = [result1.trq_multiplier];
            ylabelText = 'Torque multiplier [-]';
        case 26
            y_axis = [result1.wheelbase];
            ylabelText = 'wheelbase [mm]';
        case 27
            y_axis = [result1.h_cog];
            ylabelText = 'Height of vehicles Center of Gravity(COG)[mm]';
        case 28
            y_axis = [result1.x_cog];
            ylabelText = 'x-Coordinate of vehicles COG [mm]';
        case 29
            y_axis = [result1.x_cog_balast];
            ylabelText = 'x-Coordinate of vehicles COG [mm]';
        case 30
            y_axis = [result1.x_cog_driver];
            ylabelText = 'x-Coordinate of vehicles COG [mm]';
        case 31
            y_axis = [result1.x_va];
            ylabelText = 'x-Coordinate of the front axle [mm]';
        case 32
            y_axis = [result1.z_chaindrive];
            ylabelText = 'Number of teeth on sprocket [-]';
        case 33
            y_axis = [result1.z_sprocket];
            ylabelText = 'Number of teeth on pinion [-]';
        case 34
            y_axis = [result1.i_G];
            ylabelText = 'Gear Ratio [-]';
        case 35
            y_axis = [result1.ConstantDownforce];
            ylabelText = 'Constant Downforce [N]';
        case 36
            y_axis = [result1.DRS_Radius];
            ylabelText = 'DRS Radius [m]';
        case 37
            y_axis = [result1.vV_skidpad];
            ylabelText = 'vV_skidpad [m/s]';
        case 38
            y_axis = [result1.t_skidpad];
            ylabelText = 't_skidpad [m/s]';   
        case 39
            y_axis = sort([result1.c_l]./[result1.c_w]);
            ylabelText = 'Aero Efficiency c_l/c_w [-]';
    end
    
    if nargin == 7
        switch zparam  
            case 1
                z_axis = [result1.t_ges];
                zlabelText = 'lap time [s]';
            case 2
                z_axis = [result1.A];
                zlabelText = 'A [m²]';
            case 3
                z_axis = [result1.FB];
                zlabelText = 'Maximum Braking Force [N]';
            case 4
                z_axis = [result1.LMUX];
                zlabelText = 'Longitudinal scaling factor [-]';
            case 5
                z_axis = [result1.LMUY];
                zlabelText = 'Lateral scaling factor [-]';
            case 6
                z_axis = [result1.aero_pv];
                zlabelText = 'Aerodynamic Force on front axle [-]';
            case 7
                z_axis = [result1.c_l];
                zlabelText = 'c_l value [-]';
            case 8
                z_axis = [result1.c_w];
                zlabelText = 'c_d value [-]';
            case 9
                z_axis = [result1.camber];
                zlabelText = 'Camber [°]';
            case 10
                z_axis = [result1.downforce_multiplier];
                zlabelText = 'downforce multiplier [-]';
            case 11
                z_axis = [result1.drivetrain_eff];
                zlabelText = 'Overall efficiency of powertrain [-]';
            case 12
                z_axis = [result1.m_balast];
                zlabelText = 'mass of ballast [kg]';
            case 13
                z_axis = [result1.m_driver];
                zlabelText = 'mass of driver [kg]';
            case 14
                z_axis = [result1.m_tot];
                zlabelText = 'total mass [kg]';
            case 15
                z_axis = [result1.m_ph];
                zlabelText = 'Percentage of rear axle wheel load [%]';
            case 16
                z_axis = [result1.n_max];
                zlabelText = 'Maximum RPM of the motor [1/min]';
            case 17
                z_axis = [result1.num_motors];
                zlabelText = 'number of E-motors [-]';
            case 18
                z_axis = [result1.p_Tire];
                zlabelText = 'Tire pressure [Pa]';
            case 19
                z_axis = [result1.P_max];
                zlabelText = 'Power Limit in endurance [W]';
            case 20
                z_axis = [result1.t_L];
                zlabelText = 'Ambient air temperature [°C]';
            case 21
                z_axis = [result1.thetaV_X];
                zlabelText = 'Moment of Inertia of vehicle about X-Axis [kg*m²]';
            case 22
                z_axis = [result1.thetaV_Y];
                zlabelText = 'Moment of Inertia of vehicle about Y-Axis [kg*m²]';
            case 23
                z_axis = [result1.thetaV_Z];
                zlabelText = 'Moment of Inertia of vehicle about Z-Axis [kg*m²]';
            case 24
                z_axis = [result1.track];
                zlabelText = 'track [-]';
            case 25
                z_axis = [result1.trq_multiplier];
                zlabelText = 'Torque multiplier [-]';
            case 26
                z_axis = [result1.wheelbase];
                zlabelText = 'wheelbase [mm]';
            case 27
                z_axis = [result1.h_cog];
                zlabelText = 'Height of vehicles Center of Gravity(COG)[mm]';
            case 28
                z_axis = [result1.x_cog];
                zlabelText = 'x-Coordinate of vehicles COG [mm]';
            case 29
                z_axis = [result1.x_cog_balast];
                zlabelText = 'x-Coordinate of vehicles COG [mm]';
            case 30
                z_axis = [result1.x_cog_driver];
                zlabelText = 'x-Coordinate of vehicles COG [mm]';
            case 31
                z_axis = [result1.x_va];
                zlabelText = 'x-Coordinate of the front axle [mm]';
            case 32
                z_axis = [result1.z_chaindrive];
                zlabelText = 'Number of teeth on sprocket [-]';
            case 33
                z_axis = [result1.z_sprocket];
                zlabelText = 'Number of teeth on pinion [-]';
            case 34
                z_axis = [result1.i_G];
                zlabelText = 'Gear Ratio [-]';
            case 35
                z_axis = [result1.ConstantDownforce];
                zlabelText = 'Constant Downforce [N]';
            case 36
                z_axis = [result1.DRS_Radius];
                zlabelText = 'DRS Radius [m]';
            case 37
                z_axis = [result1.vV_skidpad];
                zlabelText = 'vV_skidpad [m/s]';
            case 38
                z_axis = [result1.t_skidpad];
                zlabelText = 't_skidpad [m/s]';   
            case 39
                z_axis = sort([result1.c_l]./[result1.c_w]);
                zlabelText = 'Aero Efficiency c_l/c_w [-]';
        end
        
        if dots
            plot3(x_axis,y_axis,z_axis,'-o','MarkerIndices',1:length(y_axis))
        else
%             [X,Y,Z] = meshgrid(x_axis,y_axis,z_axis);
%             mesh(X,Y,Z);

            % z wurde in Abh. von x und y gemessen
            % Festlegen des Gitters in x und y Koordinate
            rangeScaleFactor = 0.1;
            rangeX = min(x_axis):resx:max(x_axis);
            rangeY = min(y_axis):resy:max(y_axis);
            
%             rangeX  = -10:0.1:10;
%             rangeY  = -10:0.1:10;

            [X,Y]=meshgrid(rangeX,rangeY);
            % Interpolation der Messwerte Z an den Gitterpunkten X,Y
            Z=griddata(x_axis,y_axis,z_axis,X,Y,'cubic');
            % Plot als Fläche
            surf(X,Y,Z)
            hold off
            grid on 
            
            %plot3(x_axis,y_axis,z_axis)
        end
        
        zlabel(zlabelText,'FontSize',10)
    else
        if dots
            plot(x_axis,y_axis,'-o','MarkerIndices',1:length(y_axis))
        else
            plot(x_axis,y_axis)
        end
    end
    
    
    
    title('Sensitvity Analysis','FontSize',12)
    xlabel(xlabelText,'FontSize',10)
    ylabel(ylabelText,'FontSize',10)
    grid on
    box on
end