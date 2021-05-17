function FY = MF52_Fy_ps(ALPHA,Fz_r,GAMMA,KAPPA,TIRparam)
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
DFZ    = (Fz - FZ0PR) ./ FZ0PR; %14,  (%30)

%% Lateral Force (Pure lateral Slip)

% Coefficients - pure lateral slip       
PCY1 = TIRparam.PCY1;           % Shape faxtor Cfx for longitudinal force
PDY1 = TIRparam.PDY1;           % Longitudinal friction Mux at Fznom
PDY2 = TIRparam.PDY2;           % Variation of friction Mux with load
PDY3 = TIRparam.PDY3;           % Variation of friction Mux with camber
PEY1 = TIRparam.PEY1;           % Longitudinal curvature Efx at Fznom
PEY2 = TIRparam.PEY2;           % Variation of curvature Efx with load
PEY3 = TIRparam.PEY3;           % Variation of curvature Efx with load squared
PEY4 = TIRparam.PEY4;           % Factor in curvature Efx while driving
PKY1 = TIRparam.PKY1;           % Longitudinal slip stiffness Kfx/Fz at Fznom
PKY2 = TIRparam.PKY2;           % Variation of slip stiffness Kfx/Fz with load
PKY3 = TIRparam.PKY3;           % Exponent in slip stiffness Kfx/Fz with load
PHY1 = TIRparam.PHY1;           % Horizontal shift Shx at Fznom
PHY2 = TIRparam.PHY2;           % Variation of shift Shx with load
PHY3 = TIRparam.PHY3;           % Variation of shift Shx with load
PVY1 = TIRparam.PVY1;           % Vertical shift Svx/Fz at Fznom
PVY2 = TIRparam.PVY2;           % Variation of shift Svx/Fz with load
PVY3 = TIRparam.PVY3;           % Vertical shift Svx/Fz at Fznom
PVY4 = TIRparam.PVY4;           % Variation of shift Svx/Fz with load

% Calculation - Lateral Force (Pure lateral Slip)
GAMMAY = GAMMA .* LGAY; 
SHY = (PHY1 + PHY2 .* DFZ) .* LHY + PHY3 .* GAMMAY; 
ALPHAY = ALPHA + SHY; 
CY = PCY1 .* LCY; 
MUY = (PDY1 + PDY2 .* DFZ) .* (1 - PDY3 .* GAMMAY.^2) .* LMUY; 
DY = MUY .* Fz; 
EY = (PEY1 + PEY2 .* DFZ) .* (1 - (PEY3 + PEY4 .* GAMMAY) .* sign(ALPHAY)) .* LEY; 
KY = PKY1 .* FZ0 .* sin(2 * atan(Fz ./ (PKY2 .* FZ0 .* LFZO))) .* (1 - PKY3 .* abs(GAMMAY)) .* LFZO .* LKY; 
BY = KY ./ (CY .* DY); 
SVY = Fz .* ((PVY1 + PVY2 .* DFZ) .* LVY + (PVY3 + PVY4 .* DFZ) .* GAMMAY) .* LMUY; 

FY0 = DY .* sin(CY .* atan(BY .* ALPHAY - EY .* (BY .* ALPHAY - atan(BY .* ALPHAY)))) + SVY; 

FY = FY0;
