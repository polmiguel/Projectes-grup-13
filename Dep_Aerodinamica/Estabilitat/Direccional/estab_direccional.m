clc; clear all; close all;
format long;

% Flight conditions
cruiseAltitude = 3000; % m
v = 138; % m/s
[T,P,rho] = airConditions(cruiseAltitude);

% TAIL PARAMETERS
cr_t = 3 ; %Root chord
ct_t = 1.5 ; %Tip chord
b_t = 3; % wingspan
Sec_t = b_t*(cr_t+ct_t) / 2; % Wing section
AR_t = b_t^2/Sec_t ;   % aspect ratio
TR_t = ct_t/cr_t  ;   % taper ratio
DE25_t = 15 ; % sweep angle at c/4 (deg)
ETIP_t = 0; % tip twist (deg, negative for washout)
A0p_t = [ -1.15 -1.1 ]; % root and tip section zero-lift angles (deg)
CM0p_t = [ -0.005 -0.005 ]; % root and tip section free moments
CDP_t = [ 0.005676252644454  -0.005889143467220   0.007192733912257  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.006329508442611  -0.008927584478817   0.009711871881867 ] ;  % tip section CD0, k1 and k2

%FUSELAGE PARAMETERS
l_v = 15; %braç de palanca de estabilitzador
l_f = 15; %longitud del fuselatge
d = 2.5;          % Diámetro del fuselaje (m)
eta_v = 0.9; %eficiència de l'estabilitzador vertical
S_f = l_f * d;    % Área proyectada lateral del fuselaje (m^2)
k_f = 0.1;        % Factor empírico para el fuselaje

% WING PARAMETERS
cr = 2 ; %Root chord
ct = 1.3 ; %Tip chord
b = 22; %wingspan
Sec = b*(cr+ct) / 2; %Wing section
AR = b^2/Sec ;   % aspect ratio
TR = ct/cr  ;   % taper ratio
DE25 = 0 ; % sweep angle at c/4 (deg)
ETIP = -5; % tip twist (deg, negative for washout)
A0p = [ -1.22 -1.15 ]; % root and tip section zero-lift angles (deg)
CM0p = [ -0.0086 -0.0086 ]; % root and tip section free moments
CDP = [ 0.003656878539270  -0.002865324621379   0.006216731672035  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.004606879413827  -0.003416310943263   0.006571601776373 ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)
YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)

% Reynolds determination
mu = 1.8034e-5;
mu2 = calculateViscosity(T)/rho;
Re_root = calculateReynolds(v,cr,mu);
Re_tip = calculateReynolds(v,ct,mu);

% Simulation data (by the time being only longitudinal analysis)
N = 100 ; % number of panels along the span
ALPHA = [ -2. 0. 2. 3.4 4. 8. ] ; % angles of attack for analysis (deg) 

%% Wing 
[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

a_wb = polyfit(ALPHA, force_coeff(7,:), 1);
CL_ala = force_coeff(7,4);
cL = force_coeff(7,:);
cM = force_coeff(5,:); % Coeficiente de momento de cabeceo en el borde de ataque
coeff_cM_curve = polyfit(cL, cM, 1); % Ajuste lineal
cM0_wb = coeff_cM_curve(2); % Coeficiente de momento de cabeceo libre
mac_wb = mac;

%% Tail 
[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR_t,TR_t,N,DE25_t,ETIP_t,A0p_t,CM0p_t,CDP_t,YF_pos,CF_ratio,DE_flap,FlapCorr); 

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

a_t = polyfit(ALPHA, force_coeff(7,:), 1);
CL_cua = force_coeff(7,4);


%% Cn_beta

% Cálculo del coeficiente de volumen del estabilizador vertical
V_v = (Sec_t * l_v) / (Sec * b);

% Cálculo de la pendiente de la curva de sustentación del estabilizador
CL_alpha_cua = (2 * pi) / (1 + (2 / AR_t)); % En radianes^-1

% Contribución del estabilizador vertical
Cn_beta_vertical = eta_v * V_v * CL_alpha_cua;

% Contribución del fuselaje
Cn_beta_fuselaje = -k_f * (S_f * l_v) / (Sec * b);

% Cálculo total de Cn_beta
Cn_beta_total = Cn_beta_vertical + Cn_beta_fuselaje

