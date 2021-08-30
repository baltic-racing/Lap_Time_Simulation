%% plotAeroData.m
% Plots Aero Data from the aerodyanmic setup screen.
%
% By Eric Dornieden, Baltic Racing
% Copyright (C) 2021, Baltic Racing, all rights reserved.

function plotAeroData(A_S,c_w, c_l, downforce_multiplier, aero_pv, t_L, p_L, R_L)
%plotAeroData plots the aero data over velocity
%   Detailed explanation goes here

clf

% Speed from 0 to 30m/s in 0.5 steps
vV = (1:60);

rho_L = p_L/(R_L*(t_L+273.15));

for i = 1:60
    vV(i) = i/2;    % Velocity is (i*0.5) = Speed from 0 to 30m/s in 0.5 steps
    Faero(i) = downforce_multiplier * c_l * A_S * rho_L * vV(i)^2 / 2; 
    FL(i) = c_w * A_S * rho_L * vV(i)^2 / 2;
    disp(num2str(i));
    disp(num2str(Faero(i)));
end

plot(vV, Faero,'b', vV, FL,'r', vV, Faero*aero_pv);
title('Aero Data','FontSize',12) 
xlabel('Velocity [m/s]','FontSize',10) 
ylabel('Aero Force [N]','FontSize',10)
legend('Downforce [N]','Drag [N]','Downforce Front [N]')

annotation('textbox', [0.2, 0.8, 0.1, 0.1],...
                        'String', {['L/D = ' num2str(Faero(60)/FL(60)),],...
                        ['Downforce (@15 m/s) = ' num2str(Faero(30))],...
                        ['Drag (@15 m/s) = ' num2str(FL(30))],...
                        ['Rho = ' num2str(rho_L)]})
end

