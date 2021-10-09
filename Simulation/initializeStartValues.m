function [aRev,aVX,aVY,BPPsignal,cZ_fl,cZ_fr,cZ_rl,cZ_rr,dFWZrl_aero,dFWZrr_aero,dFWZfl_aero,dFWZfr_aero,dFWZrl_x,dFWZrr_x,dFWZfl_x,dFWZfr_x...
    ,dFWZrl_y,dFWZrr_y,dFWZfl_y,dFWZfr_y,E_Accu,E_Accu_Recu,E_heat,Faero,FB,FVX,FVX_f,FVX_fr,FVX_fl,FVX_rl,FVX_rr,FVXre,FVXid,FR,FL,Fdr,FVY,...
    FWXmax_f,FWXmax_r,FWXmax_fl,FWXmax_fr,FWXmax_rl,FWXmax_rr,FWYf,FWYr,FWYmax_f,FWYmax_r,FWYmax_fl,FWYmax_fr,FWYmax_rl,FWYmax_rr,FWZtot,FWZ_rl,...
    FWZ_rr,FWZ_fl,FWZ_fr,FWZr,FWZf,Mi,M_tractive,ni,P_M,P_el,P_tractive,Rdyn_fl,Rdyn_fr,Rdyn_rl,Rdyn_rr,slipY_f,slipY_r,t,l_contact_patch_fl,...
    l_contact_patch_fr,l_contact_patch_rl,l_contact_patch_rr,kappa_rl,kappa_rr,kappa_fl,kappa_fr,delta,beta,psi1,alpha_f,alpha_r,TC,TC_front,ABS,...
    gearSelection,Tirelimit,vAPEXmax,vV,vVYmax,A_accu_cell,P_Mloss,P_Bh,motor_eff,Capacity_Cellpack,SOC_Cellpack,Voltage_Cellpack,V_i,VirtualCurrent_Cellpack,...
    Current_Cellpack_Pointer,Energy_Cellpack,Energy_Cellpack_Total,SOC_Pointer,Current_Cellpack_Pointer_Voltage] = initializeStartValues(FB, Track, ApexIndexes)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    load('RandomizedCellData.mat');

    trackLength = length(Track);

    aRev = zeros(1,trackLength-1);
    aVX = zeros(1,trackLength-1);
    aVY = zeros(1,trackLength-1);
    BPPsignal = zeros(1,trackLength-1);
    cZ_fl = zeros(1,trackLength);
    cZ_fr = zeros(1,trackLength);
    cZ_rl = zeros(1,trackLength);
    cZ_rr = zeros(1,trackLength);
    dFWZrl_aero = zeros(1,trackLength-1);
    dFWZrr_aero = zeros(1,trackLength-1);
    dFWZfl_aero = zeros(1,trackLength-1);
    dFWZfr_aero = zeros(1,trackLength-1);
    dFWZrl_x = zeros(1,trackLength-1);
    dFWZrr_x = zeros(1,trackLength-1);
    dFWZfl_x = zeros(1,trackLength-1);
    dFWZfr_x = zeros(1,trackLength-1);
    dFWZrl_y = zeros(1,trackLength-1);
    dFWZrr_y = zeros(1,trackLength-1);
    dFWZfl_y = zeros(1,trackLength-1);
    dFWZfr_y = zeros(1,trackLength-1);

    E_Accu = zeros(1,trackLength);
    E_Accu_Recu = zeros(1,trackLength);
    E_heat = zeros(1,trackLength);

    Faero = zeros(1,trackLength-1);
    FB = FB.*ones(1,trackLength);
    FVX = zeros(1,trackLength-1);
    FVX_f = zeros(1,trackLength-1);
    FVX_fr = zeros(1,trackLength-1);
    FVX_fl = zeros(1,trackLength-1);
    FVX_rl = zeros(1,trackLength-1);
    FVX_rr = zeros(1,trackLength-1);
    FVXre = zeros(1,trackLength-1);
    FVXid = zeros(1,trackLength-1);
    FR = zeros(1,trackLength-1);
    FL = zeros(1,trackLength-1);
    Fdr = zeros(1,trackLength-1);
    FVY = zeros(1,trackLength-1);
    FWXmax_f = zeros(1,trackLength);
    FWXmax_r = zeros(1,trackLength);
    FWXmax_fl = zeros(1,trackLength);
    FWXmax_fr = zeros(1,trackLength);
    FWXmax_rl = zeros(1,trackLength);
    FWXmax_rr = zeros(1,trackLength);
    FWYf = zeros(1,trackLength-1);
    FWYr = zeros(1,trackLength-1);
    FWYmax_f = zeros(1,trackLength);
    FWYmax_r = zeros(1,trackLength);
    FWYmax_fl = zeros(1,trackLength);
    FWYmax_fr = zeros(1,trackLength);
    FWYmax_rl = zeros(1,trackLength);
    FWYmax_rr = zeros(1,trackLength);
    FWZtot = zeros(1,trackLength);
    FWZ_rl = zeros(1,trackLength);
    FWZ_rr = zeros(1,trackLength);
    FWZ_fl = zeros(1,trackLength);
    FWZ_fr = zeros(1,trackLength);
    FWZr = zeros(1,trackLength-1);
    FWZf = zeros(1,trackLength-1);
    Mi = zeros(1,trackLength-1);
    M_tractive = zeros(1,trackLength-1);
    ni = zeros(1,trackLength-1);
    P_M = zeros(1,trackLength-1);
    P_el = zeros(1,trackLength-1);
    P_tractive = zeros(1,trackLength-1);
    Rdyn_fl = zeros(1,trackLength);
    Rdyn_fr = zeros(1,trackLength);
    Rdyn_rl = zeros(1,trackLength);
    Rdyn_rr = zeros(1,trackLength);
    slipY_f = zeros(1,trackLength);
    slipY_r = zeros(1,trackLength);
    t = zeros(1,trackLength-1);

    l_contact_patch_fl = zeros(1,trackLength);
    l_contact_patch_fr = zeros(1,trackLength);
    l_contact_patch_rl = zeros(1,trackLength);
    l_contact_patch_rr = zeros(1,trackLength);

    kappa_rl = zeros(1,trackLength);
    kappa_rr = zeros(1,trackLength);
    kappa_fl = zeros(1,trackLength);
    kappa_fr = zeros(1,trackLength);

    delta = zeros(1,trackLength);
    beta = zeros(1,trackLength);
    psi1 = zeros(1,trackLength);
    alpha_f = zeros(1,trackLength);
    alpha_r = zeros(1,trackLength);

    TC = zeros(1,trackLength);
    TC_front = zeros(1,trackLength);
    ABS = zeros(1,trackLength);
    %DRS_status = zeros(1,trackLength-1);
    gearSelection = zeros(1,trackLength-1);

    Tirelimit = zeros(1,trackLength);
    vAPEXmax = zeros(1,length(ApexIndexes));
    vV = zeros(1,trackLength-1);   
    vVYmax = zeros(1,trackLength-1);
    A_accu_cell = zeros(1,trackLength);
    P_Mloss = zeros(1, trackLength);
    P_Bh = zeros(1, trackLength);
    %     M_eff_inter = zeros(l, trackLength);
    motor_eff = zeros(1, trackLength);
    Capacity_Cellpack = Parrallelcellgroups(1:131,1);
    SOC_Cellpack(1:131,1) = 1;
    Voltage_Cellpack(1:131,1) = 4.2;
    V_i = zeros(1,trackLength);
    VirtualCurrent_Cellpack = zeros(1,trackLength);
    Current_Cellpack_Pointer = zeros(1,trackLength);
    Energy_Cellpack = zeros(1,trackLength);
    Energy_Cellpack_Total = zeros(1,trackLength);
    SOC_Pointer = zeros(131,trackLength);
    Current_Cellpack_Pointer_Voltage = zeros(1,trackLength);
end

