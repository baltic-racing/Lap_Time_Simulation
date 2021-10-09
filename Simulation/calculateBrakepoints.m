function [BrakeIndexes, NonBrakeApexes, vRev] = calculateBrakepoints(FB_in, Track, ApexIndexes, vAPEXmax, m_tot, downforce_multiplier, c_l, c_d, A_S, rho_L, ConstantDownforce, c_l_DRS, DRS_status, aero_ph, aero_pv, vV, k_R, FG, h_COG, wheelbase, track, lr, lf, GAMMA, TIRparam, FWZ_fl_stat, FWZ_fr_stat, FWZ_rl_stat, FWZ_rr_stat, R, s, brakeBias_setup, brakeFunction)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

    BrakeIndexes = [];              
    NonBrakeApexes = [];

    vRev = zeros(1,length(Track)-1);                        % Initialise vRev with zeros

    k = length(ApexIndexes);

    % Reset FB
    FB = FB_in.*ones(1,length(Track));            

    if brakeFunction == 0
        while k >= 1

            Counter = 0;
            j = ApexIndexes(k);
            vRev(j-1) = vAPEXmax(k);

            if vV(j-1)>vAPEXmax(k)
                while vRev(j-1) < vV(j-1) 

%                     if vRev(j) == 0
%                         vRev(j) = vV(j);
%                     end

                    if R(ApexIndexes(k)) > 0
                        f = -1;
                    else
                        f = 1;
                    end

                    %[~, FL(j-1), Fdr(j-1), FVY(j-1), ~, ~] = calculateVehicleResistancesForces(k_R, FG, rho_L, vRev(j), c_w, A_S, m_tot, R(j), FVX(j), FVX_f(j), c_d_DRS, DRS_status(j-1), rpmpointer, n_Mmax, FB_fl(i), FB_fr(i), FB_rl(i), FB_rr(i));

                    Faero(j-1) = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vRev(j), ConstantDownforce, c_l_DRS, DRS_status(j-1));   % [N] Aerodynamic force
                    FR(j-1) = k_R*(FG+Faero(j-1));                                              % [N] Rolling Resistance
                    FL(j-1) = rho_L*vRev(j)^2/2*c_d*A_S;                                      % [N] Air resistance
                    Fdr(j-1) = FR(j-1)+FL(j-1);                                                 % [N] Overall Resistance 

                    % Wheel load transfer due to drag forces (Radlastverlagerung in Folge von Aerokräften) 
                    [dFWZrl_aero(j-1), dFWZrr_aero(j-1), dFWZfl_aero(j-1), dFWZfr_aero(j-1)] = calculateAeroforceOnWheels(Faero(j-1), aero_ph, aero_pv);

                    % First aproxmiated acceleration for further
                    % calculations
                    aVX(j-1) = (-Fdr(j-1)-FB(j-1))/m_tot;

                    % Dynamic wheel load displacement in longitudinal direction (Dynamische Radlastverlagerung in Längsrichtung = 0 angenommen)
                    [dFWZfl_x(j-1), dFWZfr_x(j-1), dFWZrl_x(j-1), dFWZrr_x(j-1)] = calculateWheelloadLongDisp(h_COG, 0, aVX(j-1), wheelbase); % Loads = 0 assumed

                    FVY(j-1) = m_tot*vRev(j)^2/R(j-1);    % [N] Centrifugal force (Zentrifugalkraft)

                    % Lenkwinkel, Schwimmwinkel, Gierrate, Schräglaufwinkel  
                    delta(j-1) = f*atan((wheelbase/1000)/sqrt(R(j-1)^2-(lr/1000)^2));   % [rad] Lenkwinkel
                    beta(j-1) = f*atan((lr/1000)/sqrt(R(j-1)^2-(lr/1000)^2));  % [rad] Schwimmwinkel
                    psi1(j-1) = vRev(j)/R(j-1);                               % [rad/s] Gierrate
                    alpha_f(j-1) = 180/pi*(delta(j-1)-(lf/1000)/vRev(j)*psi1(j-1)-beta(j-1));  % [°] Schräglaufwinkel vorne
                    alpha_r(j-1) = 180/pi*((lr/1000)/vRev(j)*psi1(j-1)-beta(j-1));           % [°] Schräglaufwinkel hinten

                    % Dynamic wheel load displacement in lateral direction (Dynamische Radlastverlagerung in Querrichtung)
                    [dFWZfl_y(j-1), dFWZfr_y(j-1), dFWZrl_y(j-1), dFWZrr_y(j-1)] = calculateWheelloadLatDisp(h_COG, track, lr, lf, wheelbase, FVY(j-1));

                    % Wheel loads (Radlasten)
                    FWZ_fl(j-1) = FWZ_fl_stat + dFWZfl_aero(j-1) + dFWZfl_x(j-1) + dFWZfl_y(j-1); % [N] Front left wheel load (Radlast vorne links)
                    FWZ_fr(j-1) = FWZ_fr_stat + dFWZfr_aero(j-1) + dFWZfr_x(j-1) + dFWZfr_y(j-1); % [N] Front right wheel load (Radlast vorne rechts)
                    FWZ_rl(j-1) = FWZ_rl_stat + dFWZrl_aero(j-1) + dFWZrl_x(j-1) + dFWZrl_y(j-1); % [N] Rear left wheel load (Radlast hinten links)
                    FWZ_rr(j-1) = FWZ_rr_stat + dFWZrr_aero(j-1) + dFWZrr_x(j-1) + dFWZrr_y(j-1); % [N] Rear right wheel load (Radlast hinten rechts)   

                    % Maximum transmissible tire forces in lateral direction (Maximal übertragbare Reifenkräfte in Querrichtung)    
                    %[FWYmax_f(j-1), FWYmax_r(j-1)] = calculateLatTireforces(FWZ_fl(j-1), FWZ_fr(j-1),FWZ_rl(j-1), FWZ_rr(j-1), GAMMA, TIRparam, alpha_f(j-1), alpha_r(j-1));

                    % Maximum transmissible tire forces in longitudinal direction (Maximal übertragbare Reifenkräfte in Längsrichtung)
                    [FWXmax_fl(j-1), FWXmax_fr(j-1), FWXmax_rl(j-1), FWXmax_rr(j-1), FWXmax_f(j-1), FWXmax_r(j-1)] = calculateLongiTireforces(FWZ_fl(j-1), FWZ_fr(j-1),FWZ_rl(j-1), FWZ_rr(j-1), GAMMA, TIRparam, alpha_f(j-1), alpha_r(j-1));

                    [~, ~, ~, ~, aRev(j-1), ~, ~] = calculateDeceleration(FB(j-1), m_tot, Fdr(j-1), FWXmax_fl(j-1), FWXmax_fr(j-1), FWXmax_rl(j-1), FWXmax_rr(j-1), brakeBias_setup);

%                     if vRev(j-1) == 0
%                         vRev(j-1) = vRev(j);
%                     end

                    vRev(j-2) = sqrt(vRev(j-1)^2-2*aRev(j-1)*(s(j-1)-s(j-2)));

                    j = j - 1;
                    Counter = Counter + 1;
                end

            end

            if Counter > 0
                BrakeIndexes = [BrakeIndexes j:ApexIndexes(k)-1];
            end

            if k > 1 && j < ApexIndexes(k-1)
                NonBrakeApexes = [NonBrakeApexes k-1];
                k = k - 1;
            end

            k = k - 1;             
        end
    else
        BrakeIndizes = [];

%         for k = 1:length(ApexIndexes)
% 
%             Counter = 0;
%             j = ApexIndizes(k);
%             vRev(j-1) = vAPEXmax(k);
% 
%             while vRev(j-1) < vV(j-1)
%                 Faero(j-1) = interp1(v,FA,vRev(j-1)*3.6,'linear','extrap'); % [N] Abtriebskraft interpoliert
%                 FR(j-1) = k_R*(FG+Faero(j-1));
%                 FL(j-1) = rho_L*vRev(j-1)^2/2*c_w*A_S;
%                 Fdr(j-1) = FR(j-1)+FL(j-1);
%                 aRev(j-1) = (-Fdr(j-1)-FB(j-1))/m_ges;
%                 vRev(j-2) = sqrt(vRev(j-1)^2-2*aRev(j-1)*(s(j-1)-s(j-2)));
%                 j = j - 1;
%                 Counter = Counter + 1;
%             end
% 
%             if Counter > 0
%                 BrakeIndizes = [BrakeIndizes j-1:ApexIndizes(k)-1];
%             end
% 
%         end
        %for k = 1:length(ApexIndexes)
        while k >= 1
            
                Counter = 0;
                j = ApexIndexes(k);
                %vRev(j-1) = vAPEXmax(k);

                while vRev(j-1) < vV(j-1)
                    
                    %[~, FL(j-1), Fdr(j-1), FVY(j-1), ~, ~] = vehicle_resistances_forces(k_R, FG, rho_L, vRev(j-1), c_w, A_S, m_tot, R(j-1), FVX(j-1), FVX_f(j-1), FB(j-1), c_d_DRS, DRS_status(j-1));
                    
                    Faero(j-1) = calculateAeroforce(downforce_multiplier, c_l, A_S, rho_L, vRev(j-1), ConstantDownforce, c_l_DRS, DRS_status(j-1));   % [N] Aerodynamic force
                    FR(j-1) = k_R*(FG+Faero(j-1));                                              % [N] Rolling Resistance
                    FL(j-1) = rho_L*vRev(j-1)^2/2*c_d*A_S;                                      % [N] Air resistance
                    Fdr(j-1) = FR(j-1)+FL(j-1);                                                 % [N] Overall Resistance 
                    
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
    end
end

