clc; clear all; close all;
format long;

%% INPUTS TAKE OFF
cruiseAltitude = 0; 
Weight = 11654*9.81; 

[T,P,rho] = airConditions(cruiseAltitude);
A0p = [ -1.355 -1.36 ]; 
CM0p = [ -0.0025 -0.0025 ]; 
CDP = [ 0.00476039394846042 -0.00369115090706098 0.00632578060880848  ;   
        0.00435255288044264 -0.00315060448108891 0.00658385993719104 ] ;  
YF_pos = [ 0.13 0.64 ]; 
CF_ratio = 0.3 ;  
DE_flap = 10; 
FlapCorr = 0.9 ;
cr = 3 ;
ct = 1.9 ;
b = 22 ;
Sec = b*(cr+ct) / 2 ;   
AR = b^2/Sec ;   
TR = ct/cr  ;   
DE25 = 0 ; 
ETIP = -5; 
N = 100 ; 
ALPHA = [ -10. -8.0 -4.0 0. 3.4 4.0 8.0 9.5 10. 12 14] ; 

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;    
[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

V_min_takeoff = sqrt(Weight*2/(force_coeff(7,8)*Sec*rho))

Clmax = 1.45; 
CL_alpha = polyfit(ALPHA,force_coeff(7,:),1);
cl1 = cl_local(:, 5);
CL1 = force_coeff(7, 5);
cl2 = cl_local(:, 6);
CL2 = force_coeff(7, 6);
Cla = (cl1 - cl2) ./ (CL1 - CL2);
Clb = cl1 - Cla * CL1;
CL_span = (Clmax - Clb)./Cla;
CL_stall= min(CL_span);

ALPHA_stall_takeoff = (CL_stall - CL_alpha(2))/CL_alpha(1)



%% INPUTS CRUISE
Power_av = 3000000*0.8;
cruiseAltitude = 3000; 
Weight = 11500*9.81;

[T,P,rho] = airConditions(cruiseAltitude);
cr = 3 ; 
ct = 1.9 ; 
b = 22; 
Sec = b*(cr+ct) / 2; 
AR = b^2/Sec ;   
TR = ct/cr  ;   
DE25 = 0 ; 
ETIP = -5; 
A0p = [ -1.22 -1.15 ]; 
CM0p = [ -0.0086 -0.0086 ]; 
CDP = [ 0.003656878539270  -0.002865324621379   0.006216731672035  ;   
        0.004606879413827  -0.003416310943263   0.006571601776373 ] ;  
YF_pos = [ 0.0 0.0]; 
CF_ratio = 0.2 ;  
DE_flap = 0.0; 
FlapCorr = 1.0 ; 

[c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 
[inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;
[GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;    
[cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

V_min_cruise = sqrt(Weight*2/(force_coeff(7,5)*Sec*rho))

Clmax = 1.6; % angle stall reynolds 6.5*10^6 aprox
CL_alpha = polyfit(ALPHA,force_coeff(7,:),1);
cl1 = cl_local(:, 5);
CL1 = force_coeff(7, 5);
cl2 = cl_local(:, 6);
CL2 = force_coeff(7, 6);
Cla = (cl1 - cl2) ./ (CL1 - CL2);
Clb = cl1 - Cla * CL1;
CL_span = (Clmax - Clb)./Cla;
CL_stall= min(CL_span);

ALPHA_stall_cruise = (CL_stall - CL_alpha(2))/CL_alpha(1)

Drag_cruise = 0.5*rho*V_min_cruise^2*0.0137*(pi*1.28^2);

RoC = (Power_av - Drag_cruise*V_min_cruise)/Weight;