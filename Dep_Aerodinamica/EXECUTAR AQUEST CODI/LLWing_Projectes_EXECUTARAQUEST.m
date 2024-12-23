% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodin√ mica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear all; close all;
format long;

% -------------------------------------------------------------------------
% INPUT DATA
% -------------------------------------------------------------------------

%Flight conditions
cruiseAltitude = 3000; % m
v = 138.8; % m/s
[T,P,rho] = airConditions(cruiseAltitude);

%Plane conditions
cl_design = 0.2;
xcg = 0;

% Wing planform (assumes planar wing)
cr = 3 ; %Root chord
ct = 1.95 ; %Tip chord
b = 22; %wingspan
Sec = b*(cr+ct) / 2; %Wing section
Span=22;

AR = b^2/Sec ;   % aspect ratio
TR = ct/cr  ;   % taper ratio
DE25 = 0 ; % sweep angle at c/4 (deg)

ETIP = -5; % tip twist (deg, negative for washout) %CANVIAR AQUEST TWIST PER AEROELASTICITAT

%Reynolds determination
mu = 1.8034e-5;
mu2 = calculateViscosity(T)/rho;
Re_root = calculateReynolds(v,cr,mu);
Re_tip = calculateReynolds(v,ct,mu);

% Sections data (uses linear interpolation between root and tip)

A0p = [ -1.22 -1.15 ]; % root and tip section zero-lift angles (deg)
CM0p = [ -0.0086 -0.0086 ]; % root and tip section free moments
%CDP CALCULATED FROM HW2_nacaPlot archive function (ROOT 2m AND TIP 1m aproximation Re)
CDP = [ 0.003656878539270  -0.002865324621379   0.006216731672035  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.004606879413827  -0.003416310943263   0.006571601776373 ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)

YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)

N = 75 ; % number of panels along the span

ALPHA = [ -2. 0. 2. 3. 3.4 4.0 8.0 ] ; % angles of attack for analysis (deg) 


LiftTotal = 0.5*rho*Sec*v^2*cl_design
% -------------------------------------------------------------------------
% LIFTING LINE SOLUTION
% -------------------------------------------------------------------------

% Wing discretization (lenghts are dimensionless with the wing span)

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using the Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

% -------------------------------------------------------------------------
% POSTPROCESS ( your work starts here...)
% -------------------------------------------------------------------------

%% PART TREBALL PROJECTES
figure;
hold on;

    plot(ALPHA,force_coeff(7,:),'DisplayName',num2str(TR))
    polyfit(force_coeff(7,:),ALPHA,1)


legend (Interpreter="latex");
legend show;
title('lift VS ALPHA curve, CRUISE');
grid on;
hold off;


figure;
hold on;

    plot(force_coeff(7,:),force_coeff(11,:),'DisplayName',num2str(TR))
    polyfit(force_coeff(11,:),force_coeff(7,:),2)


title('Total Drag (induced + profile) VS lift coefficients curve, CRUISE');
legend (Interpreter="latex");
legend show;
grid on;
hold off;

figure; hold on;

    cl1 = cl_local(:, 1);
    CL1 = force_coeff(7, 1);
    cl2 = cl_local(:, 2);
    CL2 = force_coeff(7, 2);
    
    Cla = (cl1 - cl2) ./ (CL1 - CL2);
    Clb = cl1 - Cla * CL1;
    
    spanx = linspace(-0.5*Span,0.5*Span,N);
    
    Clmax = 0.45; %tip a reynold 7M
    Cl_span = (Clmax - Clb)./Cla;
    Cl_stall = min(Cl_span);
    
    
    plot(spanx,Clb+Cl_stall.*Cla,'DisplayName',num2str(TR),'LineWidth',2);
    
    


legend (Interpreter="latex");
legend show;
yline(Clmax,LineWidth=2);
ylim([0 1]);
xlim([-11 11]);
title('Lift distribution CRUISE');
grid on;
hold off;
%% TAKE-OFF
%Flight conditions
cruiseAltitude = 10; % m
v = 55; % m/s
[T,P,rho] = airConditions(cruiseAltitude);

%Plane conditions
cl_design_to = 1.2;
xcg = 0;

% Wing planform (assumes planar wing)
cr = 3 ; %Root chord
ct = 1.95 ; %Tip chord
b = 22; %wingspan
Sec = b*(cr+ct) / 2; %Wing section
Span=22;

AR = b^2/Sec ;   % aspect ratio
TR = ct/cr  ;   % taper ratio
DE25 = 0 ; % sweep angle at c/4 (deg)

ETIP = -5; % tip twist (deg, negative for washout)

%Reynolds determination
mu = 1.8034e-5;
mu2 = calculateViscosity(T)/rho;
Re_root = calculateReynolds(v,cr,mu);
Re_tip = calculateReynolds(v,ct,mu);

% Sections data (uses linear interpolation between root and tip)

A0p = [ -1.355 -1.36 ]; % root and tip section zero-lift angles (deg)
CM0p = [ -0.0025 -0.0025 ]; % root and tip section free moments
CDP = [ 0.00476039394846042 -0.00369115090706098 0.00632578060880848  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.00435255288044264 -0.00315060448108891 0.00658385993719104 ] ;  % tip section CD0, k1 and k2

YF_pos = [ 0.13 0.64 ]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.3 ;  % flap_chord/chord ratio
DE_flap = 10; % flap deflection (deg, positive:down)
FlapCorr = 0.9 ; % flap effectiviness (<=1)

N = 100 ; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. 12 14] ; % angles of attack for analysis (deg)

LiftTotal_TakeOff = 0.5*rho*Sec*v^2*cl_design_to

% Wing discretization (lenghts are dimensionless with the wing span)

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

% Assembly of the influence coefficients matrix (needed once)

[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

% Solve circulations for different angles of attack

[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

% Loads calculation using the Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

[cl_local_takeoff,force_coeff_takeoff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

figure;
hold on;

    plot(ALPHA,force_coeff_takeoff(7,:),'DisplayName',num2str(TR))
    polyfit(force_coeff_takeoff(7,:),ALPHA,1)


legend (Interpreter="latex");
legend show;
title('lift VS ALPHA curve, TAKEOFF');
grid on;
hold off;

figure; hold on;

    cl1_takeoff = cl_local_takeoff(:, 2);
    CL1_takeoff = force_coeff_takeoff(7, 2);
    cl2_takeoff = cl_local_takeoff(:, 3);
    CL2_takeoff = force_coeff_takeoff(7, 3);
    
    Cla_takeoff = (cl1_takeoff - cl2_takeoff) ./ (CL1_takeoff - CL2_takeoff);
    Clb_takeoff = cl1_takeoff - Cla_takeoff * CL1_takeoff;
    
    spanx = linspace(-0.5*Span,0.5*Span,N);
    
    Clmax = 1.3; %tip a reynold 7M
    Cl_span_takeoff = (Clmax - Clb_takeoff)./Cla_takeoff;
    Cl_stall_takeoff = min(Cl_span_takeoff);
    
    
    plot(spanx,Clb_takeoff+Cl_stall_takeoff.*Cla_takeoff,'LineWidth',2);
    
    


legend (Interpreter="latex");
legend show;
yline(Clmax,LineWidth=2);
axis([-0.5*b 0.5*b 0 1.5]); 
title('Lift distribution TAKE-OFF');
grid on;
hold off;
%% DISTRIBUCI” LIFT PER ESTRUS
distribucioLiftPerEstrus_CRUISE = Clb+Cl_stall.*Cla 

distribucioLiftperEstrus_TAKEOFF = Clb_takeoff+Cl_stall_takeoff.*Cla_takeoff
