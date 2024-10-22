%% plotCar.m
% Plots the car suspension and calculates kinematic data such as Roll Centers, Instant Centers and so on.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function plotCar(app, CallingApp, AxesObject, params)

    % Define where to plot the car
    if nargin == 2  % If only two input arguments are given use the standard object for plotting.
        AxesObject = app.UIAxes;
    end

    % cla old plot
    cla(AxesObject);

    %% Set viewing angle
    view(AxesObject, -37.5, 30);

    if nargin == 4
        drawCoG = params.drawCoG;

        drawFrontUpright = params.drawFrontUpright;
        drawRearUpright = params.drawRearUpright;

        drawFrontTires = params.drawFrontTires;
        drawRearTires = params.drawRearTires;

        drawFrontWheelCenterAxis = params.drawFrontWheelCenterAxis;
        drawRearWheelCenterAxis = params.drawRearWheelCenterAxis;
        drawFrontRollCenter = params.drawFrontRollCenter;
        drawRearRollCenter = params.drawRearRollCenter;

        drawFrontLowerAarm = params.drawFrontLowerAarm;
        drawFrontUpperAarm = params.drawFrontUpperAarm;

        drawRearLowerAarm = params.drawRearLowerAarm;
        drawRearUpperAarm = params.drawRearUpperAarm;

        drawFrontTieRod = params.drawFrontTieRod;
        drawRearTieRod = params.drawRearTieRod;

        drawFrontInstantCenters = params.drawFrontInstantCenters;
        drawRearInstantCenters = params.drawRearInstantCenters;

        drawPitchCenter = params.drawPitchCenter;

        drawRollAxis = params.drawRollAxis;

        drawFrontRocker = params.drawFrontRocker;
        drawRearRocker = params.drawRearRocker;
        drawFrontPushPullRod = params.drawFrontPushPullRod;
        drawRearPushPullRod = params.drawRearPushPullRod;
        drawFrontCoilover = params.drawFrontCoilover;
        drawRearCoilover = params.drawRearCoilover;

        drawFrontARB = params.drawFrontARB;
        drawRearARB = params.drawRearARB;

        % Set view back to old view point
        view(AxesObject, params.viewPoint);
    else
        drawCoG = app.DrawCoGCheckBox.Value;

        drawFrontUpright = app.DrawFrontUprightCheckBox.Value;
        drawRearUpright = app.DrawRearUprightCheckBox.Value;

        drawFrontTires = app.DrawFrontTiresCheckBox.Value;
        drawRearTires = app.DrawRearTiresCheckBox.Value;

        drawFrontWheelCenterAxis = app.DrawWheelaxisCheckBox.Value;
        drawRearWheelCenterAxis = app.DrawWheelaxisCheckBox.Value;
        drawFrontRollCenter = app.DrawFrontRollCenterCheckBox.Value;
        drawRearRollCenter = app.DrawRearRollCenterCheckBox.Value;

        drawFrontLowerAarm = app.DrawFrontSuspensionCheckBox.Value;
        drawFrontUpperAarm = app.DrawFrontSuspensionCheckBox.Value;
        drawFrontTieRod = app.DrawFrontSuspensionCheckBox.Value;

        drawRearLowerAarm = app.DrawRearSuspensionCheckBox.Value;
        drawRearUpperAarm = app.DrawRearSuspensionCheckBox.Value;
        drawRearTieRod = app.DrawRearSuspensionCheckBox.Value;

        drawFrontInstantCenters = app.DrawInstantCentersFrontCheckBox.Value;
        drawRearInstantCenters = app.DrawInstantCentersRearCheckBox.Value;
        
        drawPitchCenter = app.DrawPitchCenterCheckBox.Value;

        drawRollAxis = app.DrawRollAxisCheckBox.Value;

        drawFrontRocker = app.DrawFrontSuspensionCheckBox.Value;
        drawRearRocker = app.DrawRearSuspensionCheckBox.Value;
        drawFrontPushPullRod = app.DrawFrontSuspensionCheckBox.Value;
        drawRearPushPullRod = app.DrawRearSuspensionCheckBox.Value;
        drawFrontCoilover = app.DrawFrontSuspensionCheckBox.Value;
        drawRearCoilover = app.DrawRearSuspensionCheckBox.Value;

        drawFrontARB = app.DrawFrontSuspensionCheckBox.Value;
        drawRearARB = app.DrawRearSuspensionCheckBox.Value;
    end

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

    camber_f = app.CamberfrontEditField.Value;
    camber_r = app.CamberrearEditField.Value;

    toe_f = app.ToefrontEditField.Value;
    toe_r = app.ToerearEditField.Value;

    tire_radius = app.TireradiusmmEditField.Value;
    tire_width = app.TirewidthmmEditField.Value;    
    rim_diameter = app.RimdiametermmEditField.Value;
   
    % hold UIAxes to plot other points and lines
    hold(AxesObject,'on');
    
    %% CoG
    if (drawCoG)       
        % Draw CoG
        %scatter3(AxesObject,x_cog, y_cog, z_cog);
        plot3(AxesObject,[x_cog; x_cog - 100], [0; 0], [z_cog; z_cog], 'Color', 'r');
        plot3(AxesObject,[x_cog; x_cog], [0; 100], [z_cog; z_cog], 'Color', 'g');
        plot3(AxesObject,[x_cog; x_cog], [0; 0], [z_cog; z_cog + 100], 'Color', 'b');
    end 
    
    
    %% Wheel Center axis
    if drawFrontWheelCenterAxis
        plot3(AxesObject,[x_fa; x_fa],[-track_f/2; track_f/2],[tire_radius tire_radius])
    end

    if drawRearWheelCenterAxis
        plot3(AxesObject,[x_fa+wheelbase; x_fa+wheelbase],[-track_r/2; track_r/2],[tire_radius tire_radius])
    end
    
    axis(AxesObject, "equal");
    axis(AxesObject, 'off');

    %% Plot A-Arms
    %% Front Suspension
    CHAS_LowForFront = app.CHAS_LowForFront + [0, 0, zOffset];
    CHAS_LowAftFront = app.CHAS_LowAftFront + [0, 0, zOffset];
    UPRI_LowPntFront = app.UPRI_LowPntFront + [0, 0, zOffset];

    CHAS_UppForFront = app.CHAS_UppForFront + [0, 0, zOffset];
    CHAS_UppAftFront = app.CHAS_UppAftFront + [0, 0, zOffset];
    UPRI_UppPntFront = app.UPRI_UppPntFront + [0, 0, zOffset];

    CHAS_TiePntFront = app.CHAS_TiePntFront + [0, 0, zOffset];
    UPRI_TiePntFront = app.UPRI_TiePntFront + [0, 0, zOffset];

    NSMA_PPAttPnt_L_Front = app.NSMA_PPAttPnt_L_Front + [0, 0, zOffset];
    CHAS_AttPnt_L_Front = app.CHAS_AttPnt_L_Front + [0, 0, zOffset]; % Coilover Chassis Point
    CHAS_RocAxi_L_Front = app.CHAS_RocAxi_L_Front + [0, 0, zOffset];
    CHAS_RocPiv_L_Front = app.CHAS_RocPiv_L_Front + [0, 0, zOffset];
    ROCK_RodPnt_L_Front = app.ROCK_RodPnt_L_Front + [0, 0, zOffset];
    ROCK_CoiPnt_L_Front = app.ROCK_CoiPnt_L_Front + [0, 0, zOffset];

    % ARB
    NSMA_UBarAttPnt_L_Front = app.NSMA_UBarAttPnt_L_Front + [0, 0, zOffset];
    UBAR_AttPnt_L_Front = app.UBAR_AttPnt_L_Front + [0, 0, zOffset];
    CHAS_PivPnt_L_Front = app.CHAS_PivPnt_L_Front + [0, 0, zOffset];

    %% Rear Suspension
    CHAS_LowForRear = app.CHAS_LowForRear + [0, 0, zOffset];
    CHAS_LowAftRear = app.CHAS_LowAftRear + [0, 0, zOffset];
    UPRI_LowPntRear = app.UPRI_LowPntRear + [0, 0, zOffset];

    CHAS_UppForRear = app.CHAS_UppForRear + [0, 0, zOffset];
    CHAS_UppAftRear = app.CHAS_UppAftRear + [0, 0, zOffset];
    UPRI_UppPntRear = app.UPRI_UppPntRear + [0, 0, zOffset];

    CHAS_TiePntRear = app.CHAS_TiePntRear + [0, 0, zOffset];
    UPRI_TiePntRear = app.UPRI_TiePntRear + [0, 0, zOffset];

    NSMA_PPAttPnt_L_Rear = app.NSMA_PPAttPnt_L_Rear + [0, 0, zOffset];
    CHAS_AttPnt_L_Rear = app.CHAS_AttPnt_L_Rear + [0, 0, zOffset]; % Coilover Chassis Point
    CHAS_RocAxi_L_Rear = app.CHAS_RocAxi_L_Rear + [0, 0, zOffset];
    CHAS_RocPiv_L_Rear = app.CHAS_RocPiv_L_Rear + [0, 0, zOffset];
    ROCK_RodPnt_L_Rear = app.ROCK_RodPnt_L_Rear + [0, 0, zOffset];
    ROCK_CoiPnt_L_Rear = app.ROCK_CoiPnt_L_Rear + [0, 0, zOffset];

    % ARB
    NSMA_UBarAttPnt_L_Rear = app.NSMA_UBarAttPnt_L_Rear + [0, 0, zOffset];
    UBAR_AttPnt_L_Rear = app.UBAR_AttPnt_L_Rear + [0, 0, zOffset];
    CHAS_PivPnt_L_Rear = app.CHAS_PivPnt_L_Rear + [0, 0, zOffset];

    %% Calculate Heave
    Heave = 0;

    % Heave Front
    [CHAS_LowForFront, CHAS_LowAftFront, CHAS_UppForFront, CHAS_UppAftFront, CHAS_TiePntFront] = setHeave(Heave, CHAS_LowForFront, CHAS_LowAftFront, CHAS_UppForFront, CHAS_UppAftFront, CHAS_TiePntFront, UPRI_TiePntFront, UPRI_LowPntFront, UPRI_UppPntFront);

    % Heave Rear
    [CHAS_LowForRear, CHAS_LowAftRear, CHAS_UppForRear, CHAS_UppAftRear, CHAS_TiePntRear] = setHeave(Heave, CHAS_LowForRear, CHAS_LowAftRear, CHAS_UppForRear, CHAS_UppAftRear, CHAS_TiePntRear, UPRI_TiePntRear, UPRI_LowPntRear, UPRI_UppPntRear);

    %% Plot Front Suspension
    if drawFrontLowerAarm
        % Plot fl lower A-arm
        plot3(AxesObject,[CHAS_LowForFront(1); UPRI_LowPntFront(1); CHAS_LowAftFront(1)],[CHAS_LowForFront(2); UPRI_LowPntFront(2); CHAS_LowAftFront(2)],[CHAS_LowForFront(3); UPRI_LowPntFront(3); CHAS_LowAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_LowForFront(1); UPRI_LowPntFront(1); CHAS_LowAftFront(1)],[-CHAS_LowForFront(2); -UPRI_LowPntFront(2); -CHAS_LowAftFront(2)],[CHAS_LowForFront(3); UPRI_LowPntFront(3); CHAS_LowAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        %scatter3(AxesObject,CHAS_LowFor(1),CHAS_LowFor(2),CHAS_LowFor(3)); % CHAS_LowFor
            
        %scatter3(AxesObject,CHAS_TiePnt(1),CHAS_TiePnt(2),CHAS_TiePnt(3)); % CHAS_TiePnt
    end

    if drawFrontUpperAarm
        % Plot fl upper A-arm
        plot3(AxesObject,[CHAS_UppForFront(1); UPRI_UppPntFront(1); CHAS_UppAftFront(1)],[CHAS_UppForFront(2); UPRI_UppPntFront(2); CHAS_UppAftFront(2)],[CHAS_UppForFront(3); UPRI_UppPntFront(3); CHAS_UppAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_UppForFront(1); UPRI_UppPntFront(1); CHAS_UppAftFront(1)],[-CHAS_UppForFront(2); -UPRI_UppPntFront(2); -CHAS_UppAftFront(2)],[CHAS_UppForFront(3); UPRI_UppPntFront(3); CHAS_UppAftFront(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    end

    if drawFrontTieRod
        % Plot Tie rod
        plot3(AxesObject,[CHAS_TiePntFront(1); UPRI_TiePntFront(1)],[CHAS_TiePntFront(2); UPRI_TiePntFront(2)],[CHAS_TiePntFront(3); UPRI_TiePntFront(3)],'Color','blue','LineWidth',TieRodWidth);
        plot3(AxesObject,[CHAS_TiePntFront(1); UPRI_TiePntFront(1)],[-CHAS_TiePntFront(2); -UPRI_TiePntFront(2)],[CHAS_TiePntFront(3); UPRI_TiePntFront(3)],'Color','blue','LineWidth',TieRodWidth);   
    end

    if drawFrontUpright
        % Plot front uprights
        fill3(AxesObject,[UPRI_LowPntFront(1), UPRI_UppPntFront(1), UPRI_TiePntFront(1)],[UPRI_LowPntFront(2), UPRI_UppPntFront(2), UPRI_TiePntFront(2)],[UPRI_LowPntFront(3), UPRI_UppPntFront(3), UPRI_TiePntFront(3)],[0.9290 0.6940 0.1250]);
        fill3(AxesObject,[UPRI_LowPntFront(1), UPRI_UppPntFront(1), UPRI_TiePntFront(1)],[-UPRI_LowPntFront(2), -UPRI_UppPntFront(2), -UPRI_TiePntFront(2)],[UPRI_LowPntFront(3), UPRI_UppPntFront(3), UPRI_TiePntFront(3)],[0.9290 0.6940 0.1250]);
    end

    SpringDiameter = 5;
    CoiloverColor = [0.17 0.17 0.17];

    if drawFrontPushPullRod
        % Right Side
        plot3(AxesObject,[NSMA_PPAttPnt_L_Front(1),ROCK_RodPnt_L_Front(1)],[NSMA_PPAttPnt_L_Front(2),ROCK_RodPnt_L_Front(2)],[NSMA_PPAttPnt_L_Front(3),ROCK_RodPnt_L_Front(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);

        % Left Side
        plot3(AxesObject,[NSMA_PPAttPnt_L_Front(1),ROCK_RodPnt_L_Front(1)],[-NSMA_PPAttPnt_L_Front(2),-ROCK_RodPnt_L_Front(2)],[NSMA_PPAttPnt_L_Front(3),ROCK_RodPnt_L_Front(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    end

    % Plot Rockers
    if drawFrontRocker
        % Right Rocker
        fill3(AxesObject,[CHAS_RocPiv_L_Front(1),ROCK_RodPnt_L_Front(1),ROCK_CoiPnt_L_Front(1)],[CHAS_RocPiv_L_Front(2),ROCK_RodPnt_L_Front(2),ROCK_CoiPnt_L_Front(2)],[CHAS_RocPiv_L_Front(3),ROCK_RodPnt_L_Front(3),ROCK_CoiPnt_L_Front(3)],[0.9290 0.6940 0.1250]);

        % Left Rocker
        fill3(AxesObject,[CHAS_RocPiv_L_Front(1),ROCK_RodPnt_L_Front(1),ROCK_CoiPnt_L_Front(1)],[-CHAS_RocPiv_L_Front(2),-ROCK_RodPnt_L_Front(2),-ROCK_CoiPnt_L_Front(2)],[CHAS_RocPiv_L_Front(3),ROCK_RodPnt_L_Front(3),ROCK_CoiPnt_L_Front(3)],[0.9290 0.6940 0.1250]);
    end

    % Plot Coilover
    if drawFrontCoilover
        % Right Coilover
        plot3(AxesObject,[CHAS_AttPnt_L_Front(1),ROCK_CoiPnt_L_Front(1)],[CHAS_AttPnt_L_Front(2),ROCK_CoiPnt_L_Front(2)],[CHAS_AttPnt_L_Front(3),ROCK_CoiPnt_L_Front(3)],'Color',CoiloverColor,'LineWidth',SpringDiameter);

        % Left Coilover
        plot3(AxesObject,[CHAS_AttPnt_L_Front(1),ROCK_CoiPnt_L_Front(1)],[-CHAS_AttPnt_L_Front(2),-ROCK_CoiPnt_L_Front(2)],[CHAS_AttPnt_L_Front(3),ROCK_CoiPnt_L_Front(3)],'Color',CoiloverColor,'LineWidth',SpringDiameter);
    end

    % Plot ARB (Antirollbar)
    if drawFrontARB
        plot3(AxesObject, [NSMA_UBarAttPnt_L_Front(1), UBAR_AttPnt_L_Front(1), CHAS_PivPnt_L_Front(1), CHAS_PivPnt_L_Front(1), UBAR_AttPnt_L_Front(1), NSMA_UBarAttPnt_L_Front(1)],[NSMA_UBarAttPnt_L_Front(2), UBAR_AttPnt_L_Front(2), CHAS_PivPnt_L_Front(2), -CHAS_PivPnt_L_Front(2), -UBAR_AttPnt_L_Front(2), -NSMA_UBarAttPnt_L_Front(2)],[NSMA_UBarAttPnt_L_Front(3), UBAR_AttPnt_L_Front(3), CHAS_PivPnt_L_Front(3), CHAS_PivPnt_L_Front(3), UBAR_AttPnt_L_Front(3), NSMA_UBarAttPnt_L_Front(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    end

    if drawRearLowerAarm
        % Plot rl lower A-arm
        plot3(AxesObject,[CHAS_LowForRear(1); UPRI_LowPntRear(1); CHAS_LowAftRear(1)],[CHAS_LowForRear(2); UPRI_LowPntRear(2); CHAS_LowAftRear(2)],[CHAS_LowForRear(3); UPRI_LowPntRear(3); CHAS_LowAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_LowForRear(1); UPRI_LowPntRear(1); CHAS_LowAftRear(1)],[-CHAS_LowForRear(2); -UPRI_LowPntRear(2); -CHAS_LowAftRear(2)],[CHAS_LowForRear(3); UPRI_LowPntRear(3); CHAS_LowAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        %scatter3(AxesObject,CHAS_LowFor(1),CHAS_LowFor(2),CHAS_LowFor(3)); % CHAS_LowFor
        %scatter3(AxesObject,CHAS_TiePnt(1),CHAS_TiePnt(2),CHAS_TiePnt(3)); % CHAS_TiePnt
    end

    if drawRearUpperAarm
        % Plot rl upper A-arm
        plot3(AxesObject,[CHAS_UppForRear(1); UPRI_UppPntRear(1); CHAS_UppAftRear(1)],[CHAS_UppForRear(2); UPRI_UppPntRear(2); CHAS_UppAftRear(2)],[CHAS_UppForRear(3); UPRI_UppPntRear(3); CHAS_UppAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
        plot3(AxesObject,[CHAS_UppForRear(1); UPRI_UppPntRear(1); CHAS_UppAftRear(1)],[-CHAS_UppForRear(2); -UPRI_UppPntRear(2); -CHAS_UppAftRear(2)],[CHAS_UppForRear(3); UPRI_UppPntRear(3); CHAS_UppAftRear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);  
    end

    if drawRearTieRod
        % Plot Tie rod
        plot3(AxesObject,[CHAS_TiePntRear(1); UPRI_TiePntRear(1)],[CHAS_TiePntRear(2); UPRI_TiePntRear(2)],[CHAS_TiePntRear(3); UPRI_TiePntRear(3)],'Color','blue','LineWidth',TieRodWidth);
        plot3(AxesObject,[CHAS_TiePntRear(1); UPRI_TiePntRear(1)],[-CHAS_TiePntRear(2); -UPRI_TiePntRear(2)],[CHAS_TiePntRear(3); UPRI_TiePntRear(3)],'Color','blue','LineWidth',TieRodWidth);
    end
    
    if drawRearUpright
        % Plot Rear uprights
        fill3(AxesObject,[UPRI_LowPntRear(1), UPRI_UppPntRear(1), UPRI_TiePntRear(1)],[UPRI_LowPntRear(2), UPRI_UppPntRear(2), UPRI_TiePntRear(2)],[UPRI_LowPntRear(3), UPRI_UppPntRear(3), UPRI_TiePntRear(3)],[0.9290 0.6940 0.1250]);
        fill3(AxesObject,[UPRI_LowPntRear(1), UPRI_UppPntRear(1), UPRI_TiePntRear(1)],[-UPRI_LowPntRear(2), -UPRI_UppPntRear(2), -UPRI_TiePntRear(2)],[UPRI_LowPntRear(3), UPRI_UppPntRear(3), UPRI_TiePntRear(3)],[0.9290 0.6940 0.1250]);
    end    

    SpringDiameter = 5;

    if drawRearPushPullRod
        % Right Side
        plot3(AxesObject,[NSMA_PPAttPnt_L_Rear(1),ROCK_RodPnt_L_Rear(1)],[NSMA_PPAttPnt_L_Rear(2),ROCK_RodPnt_L_Rear(2)],[NSMA_PPAttPnt_L_Rear(3),ROCK_RodPnt_L_Rear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);

        % Left Side
        plot3(AxesObject,[NSMA_PPAttPnt_L_Rear(1),ROCK_RodPnt_L_Rear(1)],[-NSMA_PPAttPnt_L_Rear(2),-ROCK_RodPnt_L_Rear(2)],[NSMA_PPAttPnt_L_Rear(3),ROCK_RodPnt_L_Rear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    end

    % Plot Rockers
    if drawRearRocker
        % Right Rocker
        fill3(AxesObject,[CHAS_RocPiv_L_Rear(1),ROCK_RodPnt_L_Rear(1),ROCK_CoiPnt_L_Rear(1)],[CHAS_RocPiv_L_Rear(2),ROCK_RodPnt_L_Rear(2),ROCK_CoiPnt_L_Rear(2)],[CHAS_RocPiv_L_Rear(3),ROCK_RodPnt_L_Rear(3),ROCK_CoiPnt_L_Rear(3)],[0.9290 0.6940 0.1250]);

        % Left Rocker
        fill3(AxesObject,[CHAS_RocPiv_L_Rear(1),ROCK_RodPnt_L_Rear(1),ROCK_CoiPnt_L_Rear(1)],[-CHAS_RocPiv_L_Rear(2),-ROCK_RodPnt_L_Rear(2),-ROCK_CoiPnt_L_Rear(2)],[CHAS_RocPiv_L_Rear(3),ROCK_RodPnt_L_Rear(3),ROCK_CoiPnt_L_Rear(3)],[0.9290 0.6940 0.1250]);
    end

    % Plot Coilover
    if drawRearCoilover
        % Right Coilover
        plot3(AxesObject,[CHAS_AttPnt_L_Rear(1),ROCK_CoiPnt_L_Rear(1)],[CHAS_AttPnt_L_Rear(2),ROCK_CoiPnt_L_Rear(2)],[CHAS_AttPnt_L_Rear(3),ROCK_CoiPnt_L_Rear(3)],'Color',CoiloverColor,'LineWidth',SpringDiameter);

        % Left Coilover
        plot3(AxesObject,[CHAS_AttPnt_L_Rear(1),ROCK_CoiPnt_L_Rear(1)],[-CHAS_AttPnt_L_Rear(2),-ROCK_CoiPnt_L_Rear(2)],[CHAS_AttPnt_L_Rear(3),ROCK_CoiPnt_L_Rear(3)],'Color',CoiloverColor,'LineWidth',SpringDiameter);
    end

    % Plot ARB (Antirollbar)
    if drawRearARB
        plot3(AxesObject, [NSMA_UBarAttPnt_L_Rear(1), UBAR_AttPnt_L_Rear(1), CHAS_PivPnt_L_Rear(1), CHAS_PivPnt_L_Rear(1), UBAR_AttPnt_L_Rear(1), NSMA_UBarAttPnt_L_Rear(1)],[NSMA_UBarAttPnt_L_Rear(2), UBAR_AttPnt_L_Rear(2), CHAS_PivPnt_L_Rear(2), -CHAS_PivPnt_L_Rear(2), -UBAR_AttPnt_L_Rear(2), -NSMA_UBarAttPnt_L_Rear(2)],[NSMA_UBarAttPnt_L_Rear(3), UBAR_AttPnt_L_Rear(3), CHAS_PivPnt_L_Rear(3), CHAS_PivPnt_L_Rear(3), UBAR_AttPnt_L_Rear(3), NSMA_UBarAttPnt_L_Rear(3)],'Color',[0.9290 0.6940 0.1250],'LineWidth',A_ArmWidth);
    end

    %% Calculate Instant Centers Front view
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

    if drawFrontInstantCenters
        % Plot Lines between Instantaneous Center points and Tire-Ground Contact points 
        plot3(AxesObject,X_IC1_TireGround_f,Y_IC1_TireGround_f,Z_IC1_TireGround_f,'b-.')
        plot3(AxesObject,X_IC2_TireGround_f,Y_IC2_TireGround_f,Z_IC2_TireGround_f,'b-.')
        
        % Plot Instant Centers
        scatter3(AxesObject,x_fa,y_IC_f,z_IC_f,'filled','SizeData', 40, 'MarkerFaceColor', 'b');
        scatter3(AxesObject,x_fa,-y_IC_f,z_IC_f,'filled','SizeData', 40, 'MarkerFaceColor', 'b');
    end

    if drawRearInstantCenters
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

    %% Calculate Instant Centers side view
    

    if drawFrontRollCenter
        % Roll Center Front
        %scatter3(AxesObject,x_fa, 0, h_rc_f);
        scatter3(AxesObject,x_fa,round(y_RC_Front,5),z_RC_Front,'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
    end

    if drawRearRollCenter
        % Roll Center Rear
        %scatter3(AxesObject,x_fa+wheelbase, 0, h_rc_r);
        scatter3(AxesObject,x_ra,round(y_RC_Rear,5),z_RC_Rear,'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
    end

    %% Roll axis
    if drawRollAxis
        % Draw Roll axis
        plot3(AxesObject,[x_fa; x_fa+wheelbase],[0; 0],[z_RC_Front; z_RC_Rear], 'Color', '#0288d1');
    end

    %% Calculate Pitch Center (1 == x, 2 == y, 3 == z)
    sideViewIC_f = calculate2DIC(CHAS_LowForFront, CHAS_LowAftFront, CHAS_UppForFront, CHAS_UppAftFront);
    sideViewIC_r = calculate2DIC(CHAS_LowForRear, CHAS_LowAftRear, CHAS_UppForRear, CHAS_UppAftRear);

    if drawPitchCenter %% ( !! Change to front and rear !! )
        %% Right Side Pitch Center ( !! Adjust Calculations when calculating with roll (left / right independent) !! )
        scatter3(AxesObject,sideViewIC_f(1),track_f/2,sideViewIC_f(2),'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
        scatter3(AxesObject,sideViewIC_r(1),track_r/2,sideViewIC_r(2),'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
    
        plot3(AxesObject,[x_fa; sideViewIC_f(1)],[track_f/2; track_f/2],[0; sideViewIC_f(2)],'b-.');
        plot3(AxesObject,[x_ra; sideViewIC_r(1)],[track_r/2; track_r/2],[0; sideViewIC_r(2)],'b-.');
        
        %% Left Side Pitch Center
        scatter3(AxesObject,sideViewIC_f(1),-track_f/2,sideViewIC_f(2),'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
        scatter3(AxesObject,sideViewIC_r(1),-track_r/2,sideViewIC_r(2),'filled','SizeData', 40, 'MarkerFaceColor', '#0288d1');
    
        plot3(AxesObject,[x_fa; sideViewIC_f(1)],[-track_f/2; -track_f/2],[0; sideViewIC_f(2)],'b-.');
        plot3(AxesObject,[x_ra; sideViewIC_r(1)],[-track_r/2; -track_r/2],[0; sideViewIC_r(2)],'b-.');
    end

    %% Calculate KPI (Kingpin Inclination)
    KPI_f = calculateKPI(UPRI_LowPntFront, UPRI_UppPntFront);
    KPI_r = calculateKPI(UPRI_LowPntRear, UPRI_UppPntRear);

    %% Calculate Caster
    Caster_f = calculateCaster(UPRI_LowPntFront, UPRI_UppPntFront);
    Caster_r = calculateCaster(UPRI_LowPntRear, UPRI_UppPntRear);

    %% Calculate Mechanical Trail
    mechanicalTrail_f = calculateMechanicalTrail(tire_radius, Caster_f);
    mechanicalTrail_r = calculateMechanicalTrail(tire_radius, Caster_r);

    %% Calculate Wheel Center Y
    WheelCenterY_f = trackCamber_f/2;
    WheelCenterY_r = trackCamber_r/2;

    %% DEBUG
    steeringAngle_f = 0;
    steeringAngle_r = 0;

    %% Plot front tires
    if drawFrontTires
        plotTires(AxesObject, tire_radius, tire_width, rim_diameter, track_f, x_fa, camber_f, toe_f, KPI_f, steeringAngle_f)
        plotTires(AxesObject, tire_radius, tire_width, rim_diameter, -track_f, x_fa, camber_f, toe_f, KPI_f, steeringAngle_f);
    end
    
    %% Plot rear tires
    if drawRearTires
        plotTires(AxesObject, tire_radius, tire_width, rim_diameter, track_r, x_fa+wheelbase, camber_r, toe_r, KPI_r, steeringAngle_r)
        plotTires(AxesObject, tire_radius, tire_width, rim_diameter, -track_r, x_fa+wheelbase, camber_r, toe_r, KPI_r, steeringAngle_r)
    end

    %% Output
    app.RollCenterheightfrontLabel.Text = "Roll Center height front " + num2str(z_RC_Front) + " [mm]";
    app.RollCenterheightrearLabel.Text = "Roll Center height rear " + num2str(z_RC_Rear) + " [mm]";
    %app.RollCenterMomentRadiusLabel.Text = "Roll Center Moment Radius " + num2str(h_rc_mr) + " [mm]";
    app.KPIfrontLabel.Text = "KPI front " + num2str(KPI_f) + " [deg]";
    app.KPIrearLabel.Text = "KPI rear " + num2str(KPI_r) + " [deg]";
    app.CasterfrontLabel.Text = "Caster front " + num2str(Caster_f) + " [deg]";
    app.CasterrearLabel.Text = "Caster rear " + num2str(Caster_r) + " [deg]";
    app.MechanicalTrailfrontLabel.Text = "Mechanical Trail front " + num2str(mechanicalTrail_f) + " [mm]";
    app.MechanicalTrailrearLabel.Text = "Mechanical Trail rear " + num2str(mechanicalTrail_r) + " [mm]";
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
    % Calculate the Pitch Center of the car
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
    trackCamber = track - (tireRadius * sin(Camber*pi/180)) / sin((90-Camber)*pi/180) * 2;
    b = abs(trackCamber - track);
    a = abs((b*sin(deg2rad(90)))/(sin(deg2rad(Camber))));
    rideHeightCamber = sqrt(a^2-2*a*b*cos(deg2rad(90-Camber))+b^2)/2;
end

function Caster = calculateCaster(UPRI_LowPnt, UPRI_UppPnt)
    A = [0,0,0];
    B = [0,1,0];
    C = [0,0,1];

    % normal vector to plane ABC
    N = cross(B-A, C-A);
    
    % angle between plane and line, equals pi/2 - angle between UPRI_LowPnt-UPRI_UppPnt and N
    alpha = abs( pi/2 - acos( dot(UPRI_UppPnt-UPRI_LowPnt, N)/norm(N)/norm(UPRI_UppPnt-UPRI_LowPnt) ) );
    
    % conversion rad to degree
    Caster = rad2deg(alpha);
end

function mechanicalTrail = calculateMechanicalTrail(tire_radius, Caster)
    % Calculate the mechanical trail based on the tire radius and the
    % caster angle of the wheel
    mechanicalTrail = tan(deg2rad(Caster))*tire_radius;
end

function plotTires(AxesObject, tire_radius, tire_width, rim_diameter, track, xOffset, camber, toe, KPI, steeringAngle)

    % Calculate R and r for given tire and rim diameter
    r = ( tire_radius*2 - rim_diameter ) / 400;  
    R = tire_radius*2 / 200 - r; 

    u = linspace(0,2*pi,100);
    v = linspace(0,2*pi,100);
    
    % Mesh the grid for the tires
    [u,v]=meshgrid(u,v);
    
    %% Calculate tire coordinates
    x = (R+r.*cos(v))*sin(u) + xOffset;
    y = tire_width / 2.*sin(v) + track/2;
    z = (R+r.*cos(v))*cos(u) + tire_radius;
       
    %% With camber
    if track > 0  % Adjusts camber for left and right side tire
        camber = camber * -1;
    end

    tire = mesh(AxesObject,x,y,z,'FaceColor','#36353a','LineStyle','none','FaceLighting','gouraud');
    rotate(tire, [1 0 0], camber, [xOffset track/2 tire_radius]);

    %% With Toe
    if track < 0  % Adjusts toe for left and right side tire
        toe = toe * -1;
    end

    if KPI == 0
        toez = 0;
    else
        toez = sin(deg2rad(90-KPI))/sin(deg2rad(KPI));
    end

    KPIvector = [1 0 toez]; % Vector for Kingpin Inclination

    rotate(tire, KPIvector, toe, [xOffset track/2 tire_radius]);

    %% With steering angle ( !! move Tie Rod and calculate angle from that !!  steering angle for left and right wheel)
    rotate(tire, KPIvector, steeringAngle, [xOffset track/2 tire_radius]);
end

function adjustTieRod(steeringAngle, UPRI_TiePnt, CHAS_TiePnt, UPRI_LowPnt)
    

    A = [0,0,0];
    B = [0,1,0];
    C = [0,0,1];

    % normal vector to plane ABC
    N = cross(B-A, C-A);
    
    % angle between plane and line, equals pi/2 - angle between UPRI_LowPnt-UPRI_UppPnt and N
    uprightAngle = abs( pi/2 - acos( dot(UPRI_TiePnt-UPRI_LowPnt, N)/norm(N)/norm(UPRI_TiePnt-UPRI_LowPnt) ) );

    uprightAngle = uprightAngle + steeringAngle;

    %steeringArmLength = 

    CHAS_TiePnt = CHAS_TiePnt + [0 10 0];
    UPRI_TiePnt = UPRI_TiePnt + [];    
end

function [CHAS_LowFor, CHAS_LowAft, CHAS_UppFor, CHAS_UppAft, CHAS_TiePnt] = setHeave(Heave, CHAS_LowFor, CHAS_LowAft, CHAS_UppFor, CHAS_UppAft, CHAS_TiePnt, UPRI_TiePnt, UPRI_LowPnt, UPRI_UppPnt)
    
    %% Chassis points movement
    CHAS_LowFor = CHAS_LowFor + [0 0 Heave];
    CHAS_LowAft = CHAS_LowAft + [0 0 Heave];
    CHAS_UppFor = CHAS_UppFor + [0 0 Heave];
    CHAS_UppAft = CHAS_UppAft + [0 0 Heave];
    CHAS_TiePnt = CHAS_TiePnt + [0 0 Heave];
    
    %% Upright point movement

%     UPRI_TiePnt = s
%     UPRI_LowPnt =
% 
%     lengthUpperAarm = CHAS_UppFor(2) - UPRI_UppPnt(2);
%     
%     UPRI_UppPnt = UPRI_UppPnt *
end

function sideViewIC = calculate2DIC(CHAS_LowFor, CHAS_LowAft, CHAS_UppFor, CHAS_UppAft) 
    x1 = CHAS_LowFor(1);
    y1 = CHAS_LowFor(3);
    x2 = CHAS_LowAft(1);
    y2 = CHAS_LowAft(3);
    x3 = CHAS_UppFor(1);
    y3 = CHAS_UppFor(3);
    x4 = CHAS_UppAft(1);
    y4 = CHAS_UppAft(3);

    p1 = polyfit([x1 x2],[y1 y2],1);
    p2 = polyfit([x3 x4],[y3 y4],1);
    
    x_intersect = fzero(@(x) polyval(p1-p2,x),1);
    z_intersect = polyval(p1,x_intersect);

    sideViewIC = [x_intersect, z_intersect];
end

function calculate3DIC(CHAS_LowFor, CHAS_LowAft, CHAS_UppFor, CHAS_UppAft, UPRI_LowPnt, UPRI_UppPnt)
    patch(CHAS_LowFor, CHAS_LowAft, UPRI_LowPnt)
    patch(CHAS_UppFor, CHAS_UppAft, UPRI_UppPnt)
end