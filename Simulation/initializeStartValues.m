function sim = initializeStartValues(sim, FB, Track, ApexIndexes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    load('RandomizedCellData.mat');

    trackLength = length(Track);

    sim.aRev = zeros(1,trackLength-1);
    sim.aVX = zeros(1,trackLength-1);
    sim.aVY = zeros(1,trackLength-1);
    sim.BPPsignal = zeros(1,trackLength-1);
    sim.cZ_fl = zeros(1,trackLength);
    sim.cZ_fr = zeros(1,trackLength);
    sim.cZ_rl = zeros(1,trackLength);
    sim.cZ_rr = zeros(1,trackLength);
    sim.dFWZrl_aero = zeros(1,trackLength-1);
    sim.dFWZrr_aero = zeros(1,trackLength-1);
    sim.dFWZfl_aero = zeros(1,trackLength-1);
    sim.dFWZfr_aero = zeros(1,trackLength-1);
    sim.dFWZrl_x = zeros(1,trackLength-1);
    sim.dFWZrr_x = zeros(1,trackLength-1);
    sim.dFWZfl_x = zeros(1,trackLength-1);
    sim.dFWZfr_x = zeros(1,trackLength-1);
    sim.dFWZrl_y = zeros(1,trackLength-1);
    sim.dFWZrr_y = zeros(1,trackLength-1);
    sim.dFWZfl_y = zeros(1,trackLength-1);
    sim.dFWZfr_y = zeros(1,trackLength-1);

    sim.E_Accu = zeros(1,trackLength);
    sim.E_Accu_Recu = zeros(1,trackLength);
    sim.E_heat = zeros(1,trackLength);

    sim.Faero = zeros(1,trackLength-1);
    sim.FB = FB.*ones(1,trackLength);
    sim.FVX = zeros(1,trackLength-1);
    sim.FVX_f = zeros(1,trackLength-1);
    sim.FVX_fr = zeros(1,trackLength-1);
    sim.FVX_fl = zeros(1,trackLength-1);
    sim.FVX_rl = zeros(1,trackLength-1);
    sim.FVX_rr = zeros(1,trackLength-1);
    sim.FVXre = zeros(1,trackLength-1);
    sim.FVXid = zeros(1,trackLength-1);
    sim.FR = zeros(1,trackLength-1);
    sim.FL = zeros(1,trackLength-1);
    sim.Fdr = zeros(1,trackLength-1);
    sim.FVY = zeros(1,trackLength-1);
    sim.FWXmax_f = zeros(1,trackLength);
    sim.FWXmax_r = zeros(1,trackLength);
    sim.FWXmax_fl = zeros(1,trackLength);
    sim.FWXmax_fr = zeros(1,trackLength);
    sim.FWXmax_rl = zeros(1,trackLength);
    sim.FWXmax_rr = zeros(1,trackLength);
    sim.FWYf = zeros(1,trackLength-1);
    sim.FWYr = zeros(1,trackLength-1);
    sim.FWYmax_f = zeros(1,trackLength);
    sim.FWYmax_r = zeros(1,trackLength);
    sim.FWYmax_fl = zeros(1,trackLength);
    sim.FWYmax_fr = zeros(1,trackLength);
    sim.FWYmax_rl = zeros(1,trackLength);
    sim.FWYmax_rr = zeros(1,trackLength);
    sim.FWZtot = zeros(1,trackLength);
    sim.FWZ_rl = zeros(1,trackLength);
    sim.FWZ_rr = zeros(1,trackLength);
    sim.FWZ_fl = zeros(1,trackLength);
    sim.FWZ_fr = zeros(1,trackLength);
    sim.FWZr = zeros(1,trackLength-1);
    sim.FWZf = zeros(1,trackLength-1);
    sim.Mi = zeros(1,trackLength-1);
    sim.M_tractive = zeros(1,trackLength-1);
    sim.ni = zeros(1,trackLength-1);
    sim.P_M = zeros(1,trackLength-1);
    sim.P_el = zeros(1,trackLength-1);
    sim.P_tractive = zeros(1,trackLength-1);
    sim.Rdyn_fl = zeros(1,trackLength);
    sim.Rdyn_fr = zeros(1,trackLength);
    sim.Rdyn_rl = zeros(1,trackLength);
    sim.Rdyn_rr = zeros(1,trackLength);
    sim.slipY_f = zeros(1,trackLength);
    sim.slipY_r = zeros(1,trackLength);
    sim.t = zeros(1,trackLength-1);

    sim.alpha_fr = zeros(1,length(Track)-1);
    sim.alpha_rl = zeros(1,length(Track)-1);
    sim.alpha_rr = zeros(1,length(Track)-1);
    sim.alpha_rl = zeros(1,length(Track)-1);

    sim.l_contact_patch_fl = zeros(1,trackLength);
    sim.l_contact_patch_fr = zeros(1,trackLength);
    sim.l_contact_patch_rl = zeros(1,trackLength);
    sim.l_contact_patch_rr = zeros(1,trackLength);

    sim.kappa_rl = zeros(1,trackLength);
    sim.kappa_rr = zeros(1,trackLength);
    sim.kappa_fl = zeros(1,trackLength);
    sim.kappa_fr = zeros(1,trackLength);

    sim.delta = zeros(1,trackLength);
    sim.beta = zeros(1,trackLength);
    sim.psi1 = zeros(1,trackLength);
    sim.alpha_f = zeros(1,trackLength);
    sim.alpha_r = zeros(1,trackLength);

    sim.TC = zeros(1,trackLength);
    sim.TC_front = zeros(1,trackLength);
    sim.ABS = zeros(1,trackLength);
    %DRS_status = zeros(1,trackLength-1);
    sim.gearSelection = zeros(1,trackLength-1);

    sim.Tirelimit = zeros(1,trackLength);
    sim.vAPEXmax = zeros(1,length(ApexIndexes));
    sim.vV = zeros(1,trackLength-1);   
    sim.vVYmax = zeros(1,trackLength-1);
    sim.A_accu_cell = zeros(1,trackLength);
    sim.P_Mloss = zeros(1, trackLength);
    sim.P_Bh = zeros(1, trackLength);
    %     M_eff_inter = zeros(l, trackLength);
    sim.motor_eff = zeros(1, trackLength);
    sim.Capacity_Cellpack = Parrallelcellgroups(1:131,1);
    sim.SOC_Cellpack(1:131,1) = 1;
    sim.Voltage_Cellpack(1:131,1) = 4.2;
    sim.V_i = zeros(1,trackLength);
    sim.VirtualCurrent_Cellpack = zeros(1,trackLength);
    sim.Current_Cellpack_Pointer = zeros(1,trackLength);
    sim.Energy_Cellpack = zeros(1,trackLength);
    sim.Energy_Cellpack_Total = zeros(1,trackLength);
    sim.SOC_Pointer = zeros(131,trackLength);
    sim.Current_Cellpack_Pointer_Voltage = zeros(1,trackLength);
end

