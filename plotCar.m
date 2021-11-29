%% plotCar.m
% Plots the car suspension and calculates kinematic data such as Roll Centers, Instant Centers and so on.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function plotCar(app, CallingApp, AxesObject)
    
    % Define where to plot the car
    if nargin == 2  % If only two input arguments are given use the standard object for plotting.
        AxesObject = app.UIAxes;
    end

    % cla old plot
    cla(AxesObject);

    %% Set viewing angle
    view(AxesObject, -37.5, 30);

    %% initlize variables
    wheelbase = app.WheelbasemmEditField.Value;

    x_fa = CallingApp.x_va;
    x_ra = x_fa + wheelbase;
    
    %% additional arguments
    zOffset = 43;
    A_ArmWidth = 3;
    TieRodWidth = 3;
    
    x_cog = CallingApp.x_cog;
    y_cog = 0;
    z_cog = CallingApp.h_cog;
    
    track_f = app.TrackfrontmmEditField.Value;
    track_r = app.TrackrearmmEditField.Value;
    
%     h_rc_f = app.RollCenterheightfrontmmEditField.Value;
%     h_rc_r = app.RollCenterheightrearmmEditField.Value;
    
    tire_radius = app.TireradiusmmEditField.Value;
    tire_width = app.TirewidthmmEditField.Value;    
   
    % hold UIAxes to plot other points and lines
    hold(AxesObject,'on');
    
    %% CoG
    if (app.DrawCoGCheckBox.Value)       
        % Draw CoG
        %scatter3(AxesObject,x_cog, y_cog, z_cog);
        plot3(AxesObject,[x_cog; x_cog - 100], [0; 0], [z_cog; z_cog], 'Color', 'r');
        plot3(AxesObject,[x_cog; x_cog], [0; 100], [z_cog; z_cog], 'Color', 'g');
        plot3(AxesObject,[x_cog; x_cog], [0; 0], [z_cog; z_cog + 100], 'Color', 'b');
    end 
    
    
    %% Wheel Center axis
    if (app.DrawWheelaxisCheckBox.Value)
        plot3(AxesObject,[x_fa; x_fa],[-track_f/2; track_f/2],[tire_radius tire_radius])
        plot3(AxesObject,[x_fa+wheelbase; x_fa+wheelbase],[-track_r/2; track_r/2],[tire_radius tire_radius])
    end
    
    %% Plot tires
    teta=-pi:0.01:pi;
    z=tire_radius*cos(teta)+tire_radius;
    
    % x_front
    x_f=tire_radius*sin(teta)+x_fa;
    
    % x_rear
    x_r=tire_radius*sin(teta)+x_fa+wheelbase;
               
    for i = 0:2:tire_width
        % Draw front tires
        if (app.DrawFrontTiresCheckBox.Value)
            % FR
            plot3(AxesObject,x_f,ones(1,numel(x_f))*track_f/2+i/2,z, 'k')
            plot3(AxesObject,x_f,ones(1,numel(x_f))*track_f/2-i/2,z, 'k')
        
            % FL
            plot3(AxesObject,x_f,ones(1,numel(x_f))*-track_f/2+i/2,z, 'k')
            plot3(AxesObject,x_f,ones(1,numel(x_f))*-track_f/2-i/2,z, 'k')
        end

        % Draw rear tires
        if (app.DrawRearTiresCheckBox.Value)
            % RR
            plot3(AxesObject,x_r,ones(1,numel(x_r))*track_r/2+i/2,z, 'k')
            plot3(AxesObject,x_r,ones(1,numel(x_r))*track_r/2-i/2,z, 'k')
        
            % RL
            plot3(AxesObject,x_r,ones(1,numel(x_r))*-track_r/2+i/2,z, 'k')
            plot3(AxesObject,x_r,ones(1,numel(x_r))*-track_r/2-i/2,z, 'k')
        end
    end
    
    axis(AxesObject, "equal");
    axis(AxesObject, 'off');

    %% Plot A-Arms
    %% Front Suspension
%     CHAS_LowForFront = [388.914, -156.064, 67.743];
%     CHAS_LowAftFront = [599.771, -153.063, 76.107];
%     UPRI_LowPntFront = [479.410, -588.899, 73.000];
% 
%     CHAS_UppForFront = [395.756, -259.077, 245.545];
%     CHAS_UppAftFront = [617.179, -287.046, 245.545];
%     UPRI_UppPntFront = [504.310, -570.000, 295.000];
% 
%     CHAS_TiePntFront = [350.172, -226.166, 72.400];
%     UPRI_TiePntFront = [435.850, -630.210, 78.825];

    CHAS_LowForFront = app.CHAS_LowForFront + [0, 0, zOffset];
    CHAS_LowAftFront = app.CHAS_LowAftFront + [0, 0, zOffset];
    UPRI_LowPntFront = app.UPRI_LowPntFront + [0, 0, zOffset];

    CHAS_UppForFront = app.CHAS_UppForFront + [0, 0, zOffset];
    CHAS_UppAftFront = app.CHAS_UppAftFront + [0, 0, zOffset];
    UPRI_UppPntFront = app.UPRI_UppPntFront + [0, 0, zOffset];

    CHAS_TiePntFront = app.CHAS_TiePntFront + [0, 0, zOffset];
    UPRI_TiePntFront = app.UPRI_TiePntFront + [0, 0, zOffset];

    %% Rear Suspension
%     CHAS_LowForRear = [1876.610, 66.501, 118.500];
%     CHAS_LowAftRear = [2151.485, 45.485, 69.946];
%     UPRI_LowPntRear = [2012.500, 588.422, 73.000];
% 
%     CHAS_UppForRear = [1881.000, 159.103, 199.991];
%     CHAS_UppAftRear = [2134.104, 153.990, 204.054];
%     UPRI_UppPntRear = [2012.500, 568.050, 283.000];
% 
%     CHAS_TiePntRear = [2171.012, 130.926, 203.734];
%     UPRI_TiePntRear = [2085.775, 567.088, 283.240];

    CHAS_LowForRear = app.CHAS_LowForRear + [0, 0, zOffset];
    CHAS_LowAftRear = app.CHAS_LowAftRear + [0, 0, zOffset];
    UPRI_LowPntRear = app.UPRI_LowPntRear + [0, 0, zOffset];

    CHAS_UppForRear = app.CHAS_UppForRear + [0, 0, zOffset];
    CHAS_UppAftRear = app.CHAS_UppAftRear + [0, 0, zOffset];
    UPRI_UppPntRear = app.UPRI_UppPntRear + [0, 0, zOffset];

    CHAS_TiePntRear = app.CHAS_TiePntRear + [0, 0, zOffset];
    UPRI_TiePntRear = app.UPRI_TiePntRear + [0, 0, zOffset];

    %% Plot Front Suspension
    if (app.DrawFrontSuspensionCheckBox.Value)
        % Plot fl lower A-arm
        plot3(AxesObject,[CHAS_LowForFront(1); UPRI_LowPntFront(1); CHAS_LowAftFront(1)],[CHAS_LowForFront(2); UPRI_LowPntFront(2); CHAS_LowAftFront(2)],[CHAS_LowForFront(3); UPRI_LowPntFront(3); CHAS_LowAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_LowForFront(1); UPRI_LowPntFront(1); CHAS_LowAftFront(1)],[-CHAS_LowForFront(2); -UPRI_LowPntFront(2); -CHAS_LowAftFront(2)],[CHAS_LowForFront(3); UPRI_LowPntFront(3); CHAS_LowAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        %scatter3(AxesObject,CHAS_LowFor(1),CHAS_LowFor(2),CHAS_LowFor(3)); % CHAS_LowFor
    
        % Plot fl upper A-arm
        plot3(AxesObject,[CHAS_UppForFront(1); UPRI_UppPntFront(1); CHAS_UppAftFront(1)],[CHAS_UppForFront(2); UPRI_UppPntFront(2); CHAS_UppAftFront(2)],[CHAS_UppForFront(3); UPRI_UppPntFront(3); CHAS_UppAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_UppForFront(1); UPRI_UppPntFront(1); CHAS_UppAftFront(1)],[-CHAS_UppForFront(2); -UPRI_UppPntFront(2); -CHAS_UppAftFront(2)],[CHAS_UppForFront(3); UPRI_UppPntFront(3); CHAS_UppAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    
        % Plot Tie rod
        plot3(AxesObject,[CHAS_TiePntFront(1); UPRI_TiePntFront(1)],[CHAS_TiePntFront(2); UPRI_TiePntFront(2)],[CHAS_TiePntFront(3); UPRI_TiePntFront(3)],'Color','blue','LineWidth',TieRodWidth);
        plot3(AxesObject,[CHAS_TiePntFront(1); UPRI_TiePntFront(1)],[-CHAS_TiePntFront(2); -UPRI_TiePntFront(2)],[CHAS_TiePntFront(3); UPRI_TiePntFront(3)],'Color','blue','LineWidth',TieRodWidth);   

        %scatter3(AxesObject,CHAS_TiePnt(1),CHAS_TiePnt(2),CHAS_TiePnt(3)); % CHAS_TiePnt
    end

    if (app.DrawFrontUprightCheckBox.Value)
        % Plot front uprights
        fill3(AxesObject,[UPRI_LowPntFront(1), UPRI_UppPntFront(1), UPRI_TiePntFront(1)],[UPRI_LowPntFront(2), UPRI_UppPntFront(2), UPRI_TiePntFront(2)],[UPRI_LowPntFront(3), UPRI_UppPntFront(3), UPRI_TiePntFront(3)],[0.9290 0.6940 0.1250]);
        fill3(AxesObject,[UPRI_LowPntFront(1), UPRI_UppPntFront(1), UPRI_TiePntFront(1)],[-UPRI_LowPntFront(2), -UPRI_UppPntFront(2), -UPRI_TiePntFront(2)],[UPRI_LowPntFront(3), UPRI_UppPntFront(3), UPRI_TiePntFront(3)],[0.9290 0.6940 0.1250]);
    end

    if (app.DrawRearSuspensionCheckBox.Value)
        % Plot fl lower A-arm
        plot3(AxesObject,[CHAS_LowForRear(1); UPRI_LowPntRear(1); CHAS_LowAftRear(1)],[CHAS_LowForRear(2); UPRI_LowPntRear(2); CHAS_LowAftRear(2)],[CHAS_LowForRear(3); UPRI_LowPntRear(3); CHAS_LowAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_LowForRear(1); UPRI_LowPntRear(1); CHAS_LowAftRear(1)],[-CHAS_LowForRear(2); -UPRI_LowPntRear(2); -CHAS_LowAftRear(2)],[CHAS_LowForRear(3); UPRI_LowPntRear(3); CHAS_LowAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        %scatter3(AxesObject,CHAS_LowFor(1),CHAS_LowFor(2),CHAS_LowFor(3)); % CHAS_LowFor
    
        % Plot fl upper A-arm
        plot3(AxesObject,[CHAS_UppForRear(1); UPRI_UppPntRear(1); CHAS_UppAftRear(1)],[CHAS_UppForRear(2); UPRI_UppPntRear(2); CHAS_UppAftRear(2)],[CHAS_UppForRear(3); UPRI_UppPntRear(3); CHAS_UppAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_UppForRear(1); UPRI_UppPntRear(1); CHAS_UppAftRear(1)],[-CHAS_UppForRear(2); -UPRI_UppPntRear(2); -CHAS_UppAftRear(2)],[CHAS_UppForRear(3); UPRI_UppPntRear(3); CHAS_UppAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    
        % Plot Tie rod
        plot3(AxesObject,[CHAS_TiePntRear(1); UPRI_TiePntRear(1)],[CHAS_TiePntRear(2); UPRI_TiePntRear(2)],[CHAS_TiePntRear(3); UPRI_TiePntRear(3)],'Color','blue','LineWidth',TieRodWidth);
        plot3(AxesObject,[CHAS_TiePntRear(1); UPRI_TiePntRear(1)],[-CHAS_TiePntRear(2); -UPRI_TiePntRear(2)],[CHAS_TiePntRear(3); UPRI_TiePntRear(3)],'Color','blue','LineWidth',TieRodWidth);

        %scatter3(AxesObject,CHAS_TiePnt(1),CHAS_TiePnt(2),CHAS_TiePnt(3)); % CHAS_TiePnt
    end
    
    if (app.DrawRearUprightCheckBox.Value)
        % Plot Rear uprights
        fill3(AxesObject,[UPRI_LowPntRear(1), UPRI_UppPntRear(1), UPRI_TiePntRear(1)],[UPRI_LowPntRear(2), UPRI_UppPntRear(2), UPRI_TiePntRear(2)],[UPRI_LowPntRear(3), UPRI_UppPntRear(3), UPRI_TiePntRear(3)],[0.9290 0.6940 0.1250]);
        fill3(AxesObject,[UPRI_LowPntRear(1), UPRI_UppPntRear(1), UPRI_TiePntRear(1)],[-UPRI_LowPntRear(2), -UPRI_UppPntRear(2), -UPRI_TiePntRear(2)],[UPRI_LowPntRear(3), UPRI_UppPntRear(3), UPRI_TiePntRear(3)],[0.9290 0.6940 0.1250]);
    end    

    %% Calculate Instant Centers
    [X1_f,Y1_f,Z1_f,X2_f,Y2_f,Z2_f,x_IC_f,y_IC_f,z_IC_f] = calculateIC(CHAS_LowForFront, CHAS_LowAftFront, UPRI_LowPntFront, CHAS_UppForFront, CHAS_UppAftFront, UPRI_UppPntFront);
    [X1_r,Y1_r,Z1_r,X2_r,Y2_r,Z2_r,x_IC_r,y_IC_r,z_IC_r] = calculateIC(CHAS_LowForRear, CHAS_LowAftRear, UPRI_LowPntRear, CHAS_UppForRear, CHAS_UppAftRear, UPRI_UppPntRear);

    % Change Instant Center direction if it gets negative
    if (y_IC_f(1) < 0)
        y_IC_f(1) = y_IC_f(1) * -1;
    end

    if (y_IC_r(1) < 0)
        y_IC_r(1) = y_IC_r(1) * -1;
    end
    
    Camber_f = app.CamberfrontEditField.Value;
    Camber_r = app.CamberrearEditField.Value;

    trackCamber_f = calculateStaticCamberOffset(Camber_f, track_f, tire_radius);
    trackCamber_r = calculateStaticCamberOffset(Camber_r, track_r, tire_radius);

    X_IC1_TireGround_f = [x_fa x_fa];
    Y_IC1_TireGround_f = [-y_IC_f trackCamber_f/2];
    Z_IC1_TireGround_f = [z_IC_f 0];
    X_IC2_TireGround_f = [x_fa x_fa];
    Y_IC2_TireGround_f = [y_IC_f -trackCamber_f/2];
    Z_IC2_TireGround_f = [z_IC_f 0];

    X_IC1_TireGround_r = [x_ra x_ra];
    Y_IC1_TireGround_r = [-y_IC_r trackCamber_r /2];
    Z_IC1_TireGround_r = [z_IC_r 0];
    X_IC2_TireGround_r = [x_ra x_ra];
    Y_IC2_TireGround_r = [y_IC_r -trackCamber_r /2];
    Z_IC2_TireGround_r = [z_IC_r 0];

    if (app.DrawInstantCentersFrontCheckBox.Value)
        % Plot Lines between Instantaneous Center points and Tire-Ground Contact points 
        plot3(AxesObject,X_IC1_TireGround_f,Y_IC1_TireGround_f,Z_IC1_TireGround_f,'b-.')
        plot3(AxesObject,X_IC2_TireGround_f,Y_IC2_TireGround_f,Z_IC2_TireGround_f,'b-.')
        
        % Plot Instant Centers
        scatter3(AxesObject,x_fa,y_IC_f,z_IC_f,'filled','SizeData', 40, 'MarkerFaceColor', 'b');
        scatter3(AxesObject,x_fa,-y_IC_f,z_IC_f,'filled','SizeData', 40, 'MarkerFaceColor', 'b');
    end

    if (app.DrawInstantCentersRearCheckBox.Value)
        % Plot Lines between Instantaneous Center points and Tire-Ground Contact points 
        plot3(AxesObject,X_IC1_TireGround_r,Y_IC1_TireGround_r,Z_IC1_TireGround_r,'b-.')
        plot3(AxesObject,X_IC2_TireGround_r,Y_IC2_TireGround_r,Z_IC2_TireGround_r,'b-.')
        
        % Plot Instant Centers
        scatter3(AxesObject,x_ra,y_IC_r,z_IC_r,'filled','SizeData', 40, 'MarkerFaceColor', 'b');
        scatter3(AxesObject,x_ra,-y_IC_r,z_IC_r,'filled','SizeData', 40, 'MarkerFaceColor', 'b');
    end


    %% Roll Centers

    % Find the intersection, Roll Center point, of the two lines
    [x_RC_Front,y_RC_Front,z_RC_Front] = calculateRC(X_IC1_TireGround_f-x_fa,Y_IC1_TireGround_f,Z_IC1_TireGround_f,X_IC2_TireGround_f-x_fa,Y_IC2_TireGround_f,Z_IC2_TireGround_f);
    [x_RC_Rear,y_RC_Rear,z_RC_Rear] = calculateRC(X_IC1_TireGround_r-x_ra,Y_IC1_TireGround_r,Z_IC1_TireGround_r,X_IC2_TireGround_r-x_fa,Y_IC2_TireGround_r,Z_IC2_TireGround_r);

    if (app.DrawFrontRollCenterCheckBox.Value)
        % Roll Center Front
        %scatter3(AxesObject,x_fa, 0, h_rc_f);
        scatter3(AxesObject,x_fa,round(y_RC_Front,5),z_RC_Front,'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
    end

    if (app.DrawRearRollCenterCheckBox.Value)
        % Roll Center Rear
        %scatter3(AxesObject,x_fa+wheelbase, 0, h_rc_r);
        scatter3(AxesObject,x_ra,round(y_RC_Rear,5),z_RC_Rear,'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
    end

    %% Roll axis
    if (app.DrawRollAxisCheckBox.Value)
        % Draw Roll axis
        plot3(AxesObject,[x_fa; x_fa+wheelbase],[0; 0],[z_RC_Front; z_RC_Rear], 'Color', '#0288d1')
    end

    %% Calculate Pitxch Center
    
    %% Calculate KPI (Kingpin Inclination)
    KPI_f = calculateKPI(UPRI_LowPntFront, UPRI_UppPntFront);
    KPI_r = calculateKPI(UPRI_LowPntRear, UPRI_UppPntRear);

    %% Calculate Caster
    Caster_f = calculateCaster(UPRI_LowPntFront, UPRI_UppPntFront);
    Caster_r = calculateCaster(UPRI_LowPntRear, UPRI_UppPntRear);

    %% Output
    app.RollCenterheightfrontLabel.Text = "Roll Center height front " + num2str(z_RC_Front) + " [mm]";
    app.RollCenterheightrearLabel.Text = "Roll Center height rear " + num2str(z_RC_Rear) + " [mm]";
    %app.RollCenterMomentRadiusLabel.Text = "Roll Center Moment Radius " + num2str(h_rc_mr) + " [mm]";
    app.KPIfrontLabel.Text = "KPI front " + num2str(KPI_f) + " [deg]";
    app.KPIrearLabel.Text = "KPI rear " + num2str(KPI_r) + " [deg]";
    app.CasterfrontLabel.Text = "Caster front " + num2str(Caster_f) + " [deg]";
    app.CasterrearLabel.Text = "Caster rear " + num2str(Caster_r) + " [deg]";
   
end

function [X1,Y1,Z1,X2,Y2,Z2,x_intersect,y_intersect,z_intersect] = calculateIC(CHAS_LowFor, CHAS_LowAft, UPRI_LowPnt, CHAS_UppFor, CHAS_UppAft, UPRI_UppPnt)
    % Combine CHAS_LowFor and CHAS_LowAft
    X1 = [UPRI_LowPnt(1),(CHAS_LowFor(1)+CHAS_LowAft(1))/2];
    Y1 = [UPRI_LowPnt(2),(CHAS_LowFor(2)+CHAS_LowAft(2))/2];
    Z1 = [UPRI_LowPnt(3),(CHAS_LowFor(3)+CHAS_LowAft(3))/2];

    X2 = [UPRI_UppPnt(1),(CHAS_UppFor(1)+CHAS_UppAft(1))/2];
    Y2 = [UPRI_UppPnt(2),(CHAS_UppFor(2)+CHAS_UppAft(2))/2];
    Z2 = [UPRI_UppPnt(3),(CHAS_UppFor(3)+CHAS_UppAft(3))/2];

    % Find Instant Center of two linwe
    p1 = polyfit(Y1,Z1,1);
    p2 = polyfit(Y2,Z2,1);
    y_intersect = fzero(@(x) polyval(p1-p2,x),1);
    z_intersect = polyval(p1,y_intersect);
    x_intersect = (UPRI_LowPnt(1) + UPRI_UppPnt(1))/2;
end

function [x_intersect,y_intersect,z_intersect] = calculateRC(X1,Y1,Z1,X2,Y2,Z2)
    
    if (Y1(1) > 0)
        Y1(1) = Y1(1) * -1;
        Y2(1) = Y2(1) * -1;
    end

    p1 = polyfit(Y1,Z1,1);
    p2 = polyfit(Y2,Z2,1);
    
    y_intersect = fzero(@(x) polyval(p1-p2,x),1);
    z_intersect = polyval(p1,y_intersect);
    x_intersect = (X1(1)+X2(1))/2;
end

function [x_intersect,y_intersect,z_intersect] = calculatePC(X1,Y1,Z1,X2,Y2,Z2)
    
end

function KPI = calculateKPI(UPRI_LowPnt, UPRI_UppPnt)
    A = [0,0,0];
    B = [1,0,0];
    C = [0,0,1];

    % normal vector to plane ABC
    N = cross(B-A, C-A);
    
    % angle between plane and line, equals pi/2 - angle between UPRI_LowPnt-UPRI_UppPnt and N
    alpha = abs( pi/2 - acos( dot(UPRI_UppPnt-UPRI_LowPnt, N)/norm(N)/norm(UPRI_UppPnt-UPRI_LowPnt) ) );
    
    % conversion rad to degree
    KPI = rad2deg(alpha);
end

function trackCamber = calculateStaticCamberOffset(Camber, track, tireRadius)
    trackCamber = track + (tireRadius * sin(Camber*pi/180)) / sin((90-Camber)*pi/180) * -2;
end

function Caster = calculateCaster(UPRI_LowPnt, UPRI_UppPnt)
    A = [0,0,0];
    B = [0,1,0];
    C = [0,0,1];

    % normal vector to plane ABC
    N = cross(B-A, C-A);
    
    % angle between plane and line, equals pi/2 - angle between UPRI_LowPnt-UPRI_UppPnt and N
    alpha = abs( pi/2 - acos( dot(UPRI_UppPnt-UPRI_LowPnt, N)/norm(N)/norm(UPRI_UppPnt-UPRI_LowPnt) ) );
    
    %// you probably want it in degrees: 
    Caster = rad2deg(alpha);
end