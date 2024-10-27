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
cruiseAltitude = 5000;
[T,P,rho] = airConditions(cruiseAltitude);

%Plane conditions
cl_design = 0.55;
xcg = 1.4;

% Wing planform (assumes planar wing)
cr = 2 ; %Root chord
ct = 1.8 ; %Tip chord
b = 20; %wingspan
Sec = b*(cr+ct) / 2; %Wing section
Span=20;

AR = b^2/Sec ;   % aspect ratio
TR = ct/cr  ;   % taper ratio
DE25 = 0 ; % sweep angle at c/4 (deg)

ETIP = 0; % tip twist (deg, negative for washout)

% Sections data (uses linear interpolation between root and tip)

A0p = [ -1.1 -0.9 ]; % root and tip section zero-lift angles (deg)
CM0p = [ 0.0033 0.0056 ]; % root and tip section free moments
CDP = [ 0.0063 -0.0015 0.0062  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.0107 -0.0068 0.0089 ] ;  % tip section CD0, k1 and k2

% Flap/aileron (symmetrical deflection)

YF_pos = [ 0.0 0.0]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.2 ;  % flap_chord/chord ratio
DE_flap = 0.0; % flap deflection (deg, positive:down)
FlapCorr = 1.0 ; % flap effectiviness (<=1)

% Simulation data (by the time being only longitudinal analysis)

N = 100 ; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. ] ; % angles of attack for analysis (deg) 

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
for (valor=0.1:0.1:2)
    cr = 2 ; %Root chord
    ct = valor ; %Tip chord
    b = 20; %wingspan
    Sec = b*(cr+ct) / 2; %Wing section
    Span=20;
    
    AR = b^2/Sec ;   % aspect ratio
    TR = ct/cr  ;   % taper ratio
    DE25 = 0 ; % sweep angle at c/4 (deg)
    

    % Wing discretization (lenghts are dimensionless with the wing span)

    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

    % Assembly of the influence coefficients matrix (needed once)

    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

    % Solve circulations for different angles of attack

    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

    % Loads calculation using the Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

    plot(ALPHA,force_coeff(7,:),'DisplayName',num2str(TR))
    polyfit(force_coeff(7,:),ALPHA,1)

end

legend (Interpreter="latex");
legend show;
title('lift VS ALPHA curve, different TR values');
grid on;
hold off;


figure;
hold on;
for (valor=0.1:0.1:2)
    
    cr = 2 ; %Root chord
    ct = valor ; %Tip chord
    b = 20; %wingspan
    Sec = b*(cr+ct) / 2; %Wing section
    Span=20;
    
    AR = b^2/Sec ;   % aspect ratio
    TR = ct/cr  ;   % taper ratio
    DE25 = 0 ; % sweep angle at c/4 (deg)

    % Wing discretization (lenghts are dimensionless with the wing span)

    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

    % Assembly of the influence coefficients matrix (needed once)

    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

    % Solve circulations for different angles of attack

    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

    % Loads calculation using the Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

    plot(force_coeff(7,:),force_coeff(11,:),'DisplayName',num2str(TR))
    polyfit(force_coeff(11,:),force_coeff(7,:),2)

end
title('Total Drag (induced + profile) VS lift coefficients curve, different TR values');
legend (Interpreter="latex");
legend show;
grid on;
hold off;

figure; hold on;
for (valor=0:0.5:2)
    cr = 2 ; %Root chord
    ct = valor ; %Tip chord
    b = 20; %wingspan
    Sec = b*(cr+ct) / 2; %Wing section
    Span=20;
    
    AR = b^2/Sec ;   % aspect ratio
    TR = ct/cr  ;   % taper ratio
    DE25 = 0 ; % sweep angle at c/4 (deg)
    

    % Wing discretization (lenghts are dimensionless with the wing span)

    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

    % Assembly of the influence coefficients matrix (needed once)

    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

    % Solve circulations for different angles of attack

    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

    % Loads calculation using the Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

    cl1 = cl_local(:, 1);
    CL1 = force_coeff(7, 1);
    cl2 = cl_local(:, 2);
    CL2 = force_coeff(7, 2);
    
    Cla = (cl1 - cl2) ./ (CL1 - CL2);
    Clb = cl1 - Cla * CL1;
    
    spanx = linspace(-0.5*Span,0.5*Span,N);
    
    Clmax = 1.6267; %angle de 16.5 deg BUSCAR !!!!!!!!!!!!!!!
    Cl_span = (Clmax - Clb)./Cla;
    Cl_stall = min(Cl_span);
    
    
    plot(spanx,Clb+Cl_stall.*Cla,'DisplayName',num2str(TR));
    

end
legend (Interpreter="latex");
legend show;
yline(1.6267);
title('Lift distribution revisar');
grid on;
hold off;


figure;
hold on;
for (valor=0:0.5:5)
    
    cr = 2 ; %Root chord
    ct = 0.8 ; %Tip chord
    b = 20; %wingspan
    Sec = b*(cr+ct) / 2; %Wing section
    Span=20;
    
    AR = b^2/Sec ;   % aspect ratio
    TR = ct/cr  ;   % taper ratio
    DE25 = 0 ; % sweep angle at c/4 (deg)

    ETIP = valor; % tip twist (deg, negative for washout)

    % Wing discretization (lenghts are dimensionless with the wing span)

    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

    % Assembly of the influence coefficients matrix (needed once)

    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

    % Solve circulations for different angles of attack

    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

    % Loads calculation using the Kutta-Joukowsky theorem (costlier, but general... and illustrative!)

    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;
    
    % plot(force_coeff(7,:),force_coeff(11,:),'DisplayName',num2str(ETIP))
    % polyfit(force_coeff(11,:),force_coeff(7,:),2)

     cl1 = cl_local(:, 1);
    CL1 = force_coeff(7, 1);
    cl2 = cl_local(:, 2);
    CL2 = force_coeff(7, 2);

    Cla = (cl1 - cl2) ./ (CL1 - CL2);
    Clb = cl1 - Cla * CL1;

    spanx = linspace(-0.5*Span,0.5*Span,N);

    Clmax = 1.6267; %angle de 16.5 deg BUSCAR !!!!!!!!!!!!!!!
    Cl_span = (Clmax - Clb)./Cla;
    Cl_stall = min(Cl_span);


    plot(spanx,Clb+Cl_stall.*Cla,'DisplayName',num2str(ETIP));

end
title('Total Drag (induced + profile) VS lift coefficients curve, different ETIP values');
legend (Interpreter="latex");
legend show;
grid on;
hold off;
