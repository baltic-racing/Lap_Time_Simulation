function FX = MF52_Fx_ps(ALPHA,Fz_r,GAMMA,KAPPA,TIRparam)
   %% Inputs (Parameters)

    % VERTICAL
    FZ0   = TIRparam.FNOMIN;            % [N] Nominal wheel load
    R0    = TIRparam.UNLOADED_RADIUS;   % [m] Unloaded tire radius
    
    % Vertical Force Range
    FZMIN = TIRparam.FZMIN;             % [N] Minimum wheel load
    FZMAX = TIRparam.FZMAX;             % [N] Maximum wheel load
    
    % SCALING_COEFFICIENTS
    LFZO  = TIRparam.LFZO;          % Scale factor of nominal (rated) load
    LCX   = TIRparam.LCX;           % Scale factor of Fx shape factor
	LMUX  = TIRparam.LMUX;          % Scale factor of Fx peak friction coefficient  	
    LEX   = TIRparam.LEX;           % Scale factor of Fx curvature factor
    LKX   = TIRparam.LKX;           % Scale factor of Fx slip stiffness
    LHX   = TIRparam.LHX;           % Scale factor of Fx horizontal shift
    LVX   = TIRparam.LVX;           % Scale factor of Fx vertical shift
    LGAX  = TIRparam.LGAX;          % Scale factor of camber for Fx
    LCY   = TIRparam.LCY;			% Scale factor of Fy shape factor	
    LMUY  = TIRparam.LMUY;          % Scale factor of Fy peak friction coefficient
    LEY   = TIRparam.LEY;           % Scale factor of Fy curvature factor
    LKY   = TIRparam.LKY;           % Scale factor of Fy cornering stiffness
    LHY   = TIRparam.LHY;           % Scale factor of Fy horizontal shift
    LVY   = TIRparam.LVY;           % Scale factor of Fy vertical shift
    LGAY  = TIRparam.LGAY;          % Scale factor of camber for Fy
    LTR   = TIRparam.LTR;           % Scale factor of Peak of pneumatic trail
    LRES  = TIRparam.LRES;          % Scale factor for offset of residual torque
    LGAZ  = TIRparam.LGAZ;          % Scale factor of camber for Mz
    LXAL  = TIRparam.LXAL;          % Scale factor of alpha influence on Fx
    LYKA  = TIRparam.LYKA;			% Scale factor of alpha influence on Fy
    LVYKA = TIRparam.LYKA;			% Scale factor of kappa induced Fy
    LS    = TIRparam.LS;			% Scale factor of Moment arm of Fx	
    LSGKP = TIRparam.LSGKP;			
    LSGAL = TIRparam.LSGAL;			
    LGYR  = TIRparam.LGYR;			
    LMX   = TIRparam.LMX;			% Scale factor of overturning couple	
    LVMX  = TIRparam.LVMX;          % Scale factor of Mx vertical shift
    LMY   = TIRparam.LMY;           % Scale factor of rolling resistance torque

    ALPHA  =  ALPHA.*pi/180;    % [rad] Slip angle
    Fz     =  Fz_r;             % [N]   Wheel load
    GAMMA  =  GAMMA.*pi/180;    % [rad] Camber angle
    
GAMMAY = GAMMA .* LGAY; %31 (%48 lgay=lg
GAMMAZ = GAMMA .* LGAZ; %47 (%63 lgaz = lg
FZ0PR  = FZ0  .*  LFZO; %15,  NEED LFZO NOT LFZ0 TO MATCH TIRE PROP FILE
DFZ    = (Fz-FZ0PR) ./ FZ0PR; %14,  (%30)

%% Longitudinal Force (Pure longitudinal Slip)

% Coefficients - pure longitudinal slip       
PCX1 = TIRparam.PCX1;           % Shape faxtor Cfx for longitudinal force
PDX1 = TIRparam.PDX1;           % Longitudinal friction Mux at Fznom
PDX2 = TIRparam.PDX2;           % Variation of friction Mux with load
PDX3 = TIRparam.PDX3;           % Variation of friction Mux with camber
PEX1 = TIRparam.PEX1;           % Longitudinal curvature Efx at Fznom
PEX2 = TIRparam.PEX2;           % Variation of curvature Efx with load
PEX3 = TIRparam.PEX3;           % Variation of curvature Efx with load squared
PEX4 = TIRparam.PEX4;           % Factor in curvature Efx while driving
PKX1 = TIRparam.PKX1;           % Longitudinal slip stiffness Kfx/Fz at Fznom
PKX2 = TIRparam.PKX2;           % Variation of slip stiffness Kfx/Fz with load
PKX3 = TIRparam.PKX3;           % Exponent in slip stiffness Kfx/Fz with load
PHX1 = TIRparam.PHX1;           % Horizontal shift Shx at Fznom
PHX2 = TIRparam.PHX2;           % Variation of shift Shx with load
PVX1 = TIRparam.PVX1;           % Vertical shift Svx/Fz at Fznom
PVX2 = TIRparam.PVX2;           % Variation of shift Svx/Fz with load

% Calculation - Longitudinal Force (Pure longitudinal Slip)
GAMMAX = GAMMA .* LGAX;
SHX = (PHX1 + PHX2 .* DFZ) .* LHX;
KAPPAX = KAPPA + SHX;
CX = PCX1 .* LCX;
MUX = (PDX1 + PDX2 .* DFZ) .* (1 - PDX3 .* GAMMAX.^2) .* LMUX;
DX = MUX .* Fz;
EX = (PEX1 + PEX2 .* DFZ + PEX3 .* DFZ.^2) .* (1 - PEX4 .* sign(KAPPAX)) .* LEX;
KX = Fz .* (PKX1 + PKX2 .* DFZ) .* exp(PKX3 .* DFZ) .* LKX;
BX = KX ./ (CX .* DX);
SVX = Fz .* (PVX1 + PVX2 .* DFZ) .* LVX .* LMUX;

FX0 = DX .* sin(CX .* atan(BX .* KAPPAX - EX .* (BX .* KAPPAX - atan(BX .* KAPPA)))) + SVX;

FX = FX0;
