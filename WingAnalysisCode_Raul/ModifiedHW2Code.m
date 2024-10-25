% HOMEWORK 2 %
clc; clear; close all;

% Import data:

formatSpec = "%f";
polarsFileID = fopen("AirfoilData\polarData23018AscentVFinal.txt","r");
airfoilPolars = (fscanf(polarsFileID,formatSpec,[7 inf]))';

% Extract columns data:

Polars.alfa_CL = airfoilPolars(:,[1 2]);
Polars.CL_CD = airfoilPolars(:,[2 3]);
Polars.alfa_CM = airfoilPolars(:,[1 5]);

linearLimit = Polars.alfa_CM(:,1) >= -5 & Polars.alfa_CM(:,1) <= 10;

alfaCMLinreg = polyfit(Polars.alfa_CM(linearLimit,1),Polars.alfa_CM(linearLimit,2),1);
Cm0 = alfaCMLinreg(2);

% Plots

figure()
plot(Polars.alfa_CL(:,1), Polars.alfa_CL(:,2))

figure()
plot(Polars.CL_CD(:,1), Polars.CL_CD(:,2))

figure()
hold on
plot(Polars.alfa_CM(:,1), Polars.alfa_CM(:,2))
plot(Polars.alfa_CM(linearLimit,1),polyval(alfaCMLinreg,Polars.alfa_CM(linearLimit,1)))


Cm0_tip = Cm0;
Cm0_root = Cm0;

%% Wing Data

b=20; % wingspan [m]
ct=2.0; % tip chord [m]
cr=2.0; % root chord [m]
c_av=(ct+cr)*0.5; % average chord [m]
XCG=0.6; % distance from LE to center of gravity [m]
TR=ct/cr; % Taper Ratio [-]
S=(b/2)*(cr+ct); % wing surface [m2]
AR=b^2/S; % Aspect Ratio [-]
DE25=0; % sweep angle at c/4 [deg]
ETIP=0; % tip twist [deg]

% Flap/aileron (symmetrical deflection)

YF_pos = [0.0 0.0];     % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.0;         % flap_chord/chord ratio
DE_flap = 0.0;          % flap deflection (deg, positive:down)
FlapCorr = 0.0;         % flap effectiviness (<=1)

ALPHA=linspace(-12,16,29); % angles for numerical ansys

%% Preliminary Wing Calculations

% 2.1. C_L and alfa_L0 calculation
[ force_coeff, cl_local ,mac, ~ ] = LLWing_function(AR,TR,DE25,ETIP,YF_pos ,CF_ratio ,DE_flap ,FlapCorr,ALPHA,Cm0_tip,Cm0_root);
reg_CL=polyfit(ALPHA,force_coeff(7,:),1); %in force_coeff(7,:) we find each CL(alpha)
alpha_L0=-reg_CL(2)/reg_CL(1);
CL_alpha=reg_CL(1); % derivative of CL(alpha)

figure;
hold on
plot(ALPHA, force_coeff(7,:));
plot(ALPHA, polyval(reg_CL, ALPHA), '-' ,'LineWidth',1) ;
xlabel('alpha [deg]');
ylabel('CL');
title('Lift coefficient curve');
grid on;

% 2.2. dCM_le/dCL and CM0 and Xac
reg_CM=polyfit(force_coeff(7,:),force_coeff(5,:),1); %CM (row 5) is a function of CL (row 7)
dCM_dCL=reg_CM(1);
CM0=-reg_CM(2)/reg_CM(1);
XAC=-dCM_dCL*mac*b; %mean aerdynamic chord (mac) is dimentionless in respect to wingspan

figure;
plot(force_coeff(5,:), polyval(reg_CM, force_coeff(5,:)), '-' ,'LineWidth',1) ;
xlabel('CL');
ylabel('CM_LE');
title('Moment over leading edge coefficient (CL)');
grid on;
%xlim([-0.4 0.22])

% 2.3. Basic and additional lift coefs
Cl_1=cl_local(:,1);
Cl_2=cl_local(:,2);
CL_1=force_coeff(7,1);
CL_2=force_coeff(7,2);
Cl_a=(Cl_1-Cl_2)/(CL_1-CL_2);
Cl_b=Cl_1-Cl_a*CL_1;

figure; hold on;
plot(Cl_b,'r','LineWidth', 1);
plot(Cl_a,'b', 'LineWidth',1);
xlabel('Wingspan (m)');
ylabel('Lift coefficient Cl');
title('Lift distributions');
legend('Cl basic', 'Cl additional', 'Location','east');
grid on

% 2.4. Parabolic drag curve
reg_CD=polyfit(force_coeff(7,:),force_coeff(11,:),2); %CL row 7 and CD row 11

figure; hold on;
plot(force_coeff(7,:),force_coeff(11,:),'o');
plot(force_coeff(7,:), polyval(reg_CD, force_coeff(7,:)), '-' ,'LineWidth',1) ;
xlabel('Wing Lift Coefficient CL');
ylabel('Wing Drag Coefficient CD');
title('Quadratic wing drag coefficient');
grid on;
%ylim([0,0.03]);

% 2.5. CM_cg vs CL curve of the wing
CM_cg=-CM0-force_coeff(7,:)*(XCG-XAC)/(mac*b);

figure;
plot(force_coeff(7,:),CM_cg,'-','LineWidth', 1); %CM vs CL
xlabel('Wing Lift Coefficient CL');
ylabel('Moment about center of gravity coefficient (CM_cg)');
title('Moment coefficient about CG in front of CL');
grid on


%% Adaptation Analysis

% 3.1. Adjust of wing sweep to guarantee stability margin

DE25_vector = [15 15.5 16 16.5 17 17.5 18 18.5 19 19.5 20];  % Sweep angles to test
RES=zeros(length(DE25_vector),2);

for i = 1:length(DE25_vector)
    DE25 = DE25_vector(i);
    [ force_coeff, ~ ,mac, ~ ] = LLWing_function(AR,TR,DE25,ETIP,YF_pos ,CF_ratio ,DE_flap ,FlapCorr,ALPHA,Cm0_tip,Cm0_root);
    reg_CM=polyfit(force_coeff(7,:),force_coeff(5,:),1);
    dCM_dCL=reg_CM(1);
    XAC=-dCM_dCL*mac*b;
    Distance_percent = ((XAC -XCG)/(mac*b)) * 100;    % Distance between Xac and x_cg in percent of mac. Should be between 10% and 20%
    RES(i,1) = DE25;                                  % Results matrix with sweep angles in col. 1 and the respective Xac-x_cg distances in col.2
    RES(i,2) = Distance_percent;
end

figure
plot(RES(:,1),RES(:,2),'b',LineWidth=2);
yline(20,'--k')
yline(10,'--k','Acceptable Stability Margin Range', 'Interpreter','latex')
grid on
title(["Evoluci\'o del marge"; " d'estabilitat amb l'angle de fletxa"], 'Interpreter', 'latex',FontSize=15)
xlabel("\textit{Sweep Angle} $\Lambda$ [deg]",'Interpreter','latex',FontSize=15)
ylabel("\textit{Stability Margin} [\%]", 'Interpreter','latex',FontSize=15)

% Xac with original sweep is slightly too far back ((Xac -x_cg)/(mac*b)) *
% 100 = 20.7966 % Sweep angle has to be lower to bring Xac forward
% Analizing RES it is possible to say that, for example, 16.5 deg or 16 deg
% would be more appropiate angles.

DE25 = 16.5;

% 3.2. Washout needed to trim the wing (CM_cg=0)

W_S = 20;                                                   % Wing Load (W/S) [kgf/m^2]
CL_design = (W_S*9.81)/(1.225*((80*(1000/3600))^2)/2);      % Design CL // (20 kgf/m^2 * 9,81 N / 1 kgf) = L/S --> L = CL * q_inf * S  --> L/S = CL * q_inf  --> CL = (L/S)/q_inf
Alpha_design = (CL_design/CL_alpha)+alpha_L0;               % Design Alpha, use the CL-alpha curve obtained earlier to get alpha for the design CL

ETIP_vector = -30:1:30;                            % Washout/Washin angles to test
CM_cg_etip = zeros(length(ETIP_vector),1);           % CM_cg vector for alpha=alpha_design for each ETIP angle tested 
ALPHA = [-7 -6 -5 -4 -3 -2 -1 0 1 2 3 4 Alpha_design 6 7];

for i = 1:length(ETIP_vector)
    ETIP = ETIP_vector(i);
    [ force_coeff, ~ ,mac, ~ ] = LLWing_function(AR,TR,DE25,ETIP,YF_pos ,CF_ratio ,DE_flap ,FlapCorr,ALPHA,Cm0_tip,Cm0_root);
    reg_CM=polyfit(force_coeff(7,:),force_coeff(5,:),1);
    dCM_dCL=reg_CM(1);
    CM0=-reg_CM(2)/reg_CM(1);
    XAC=-dCM_dCL*mac*b;
    CM_cg_etip(i) = -CM0-force_coeff(7,13)*(XCG-XAC)/(mac*b);
end

regCMcgETIP = polyfit(CM_cg_etip,ETIP_vector,2);                        % 2nd order regression to fit the data obtained
ETIP_CMcg0 = regCMcgETIP(3);                    

figure
hold on
plot(-0.5:0.001:0.5,polyval(regCMcgETIP,-0.5:0.001:0.5),'r')
plot(CM_cg_etip,ETIP_vector,'bo')
yline(ETIP_CMcg0,'--k','$\epsilon$ at $C_{M,cg}=0$','Interpreter','latex')
grid on
title(["Evoluci\'o del moment al voltant"; " del centre de gravetat amb la torsi\'o"], 'Interpreter', 'latex','FontSize',15)
ylabel("\textit{Twist angle} $\epsilon$ [deg]",'Interpreter','latex',FontSize=15)
xlabel("$C_{M,cg}$", 'Interpreter','latex',FontSize=15)
legend("Regressi\'o","Dades calculades",'Interpreter','latex')

ETIP = ETIP_CMcg0;

% Comprovar que CM_cg = 0 i que la distància xac-xcg esta entre 10% i 20%
% per ETIP i DE25 trobats

[ force_coeff, ~ ,mac, ~ ] = LLWing_function(AR,TR,DE25,ETIP,YF_pos ,CF_ratio ,DE_flap ,FlapCorr,ALPHA,Cm0_tip,Cm0_root);
reg_CM=polyfit(force_coeff(7,:),force_coeff(5,:),1); 
dCM_dCL=reg_CM(1);
CM0=-reg_CM(2)/reg_CM(1);
XAC=-dCM_dCL*mac*b; 
CM_cg=-CM0-force_coeff(7,13)*(XCG-XAC)/(mac*b)
Distance_percent = ((XAC -XCG)/(mac*b)) * 100
ETIP
DE25


% 3.3. Wings CL_max and stall speed

[ force_coeff, cl_local ,mac, N ] = LLWing_function(AR,TR,DE25,ETIP,YF_pos ,CF_ratio ,DE_flap ,FlapCorr,ALPHA,Cm0_tip,Cm0_root);

Cl_1=cl_local(:,1);
Cl_2=cl_local(:,2);
CL_1=force_coeff(7,1);
CL_2=force_coeff(7,2);
Cl_a=(Cl_1-Cl_2)/(CL_1-CL_2);
Cl_b=Cl_1-Cl_a*CL_1;

figure; hold on;
plot(ALPHA,cl_local([1 25 50 75 100],:),'LineWidth', 1);
xlabel('Alpha (deg)');
ylabel('cl_y');
title('Cl local at some specific wing stations');
grid on

Cl_max = max(root.Cl);
Cl_max = Cl_max.*ones(N,1);

CL_y = (Cl_max - Cl_b)./Cl_a;
CL_max = min(CL_y)*ones(N,1);

Cl_y_CLmax = Cl_b + Cl_a .* CL_max;

v_max = sqrt((2/(1.225*CL_max))*W_S*9.81);
V_max = v_max * 3600/1000;

figure; 
subplot(2,1,1)
hold on;
ylim([-0.4 2.2]);
plot(linspace(-100,100,100),Cl_y_CLmax,'b','LineWidth', 1);
plot(linspace(-100,100,100),Cl_a,'k','LineWidth', 1);
plot(linspace(-100,100,100),Cl_b,'g','LineWidth', 1);
yline(Cl_max(1),'--b',"$C_{l,max}$",'Interpreter','latex','FontSize', 15,'LineWidth', 1);
xlabel('Wing span [percent of span]','Interpreter','latex','FontSize',15);
ylabel('Lift coefficient $C_l$','Interpreter','latex','FontSize',15);
legend('$C_l(y) \quad at \quad C_{L,max}$','$C_l,a$','$C_l,b$', 'Interpreter','latex')
grid on
subplot(2,1,2)
ylim([-0.4 2.2]);
hold on;
plot(linspace(-100,100,100),CL_y,'r','LineWidth', 1);
yline(CL_max(1),'--r',"$C_{L,max}$",'Interpreter','latex','FontSize', 15,'LineWidth', 1);
xlabel('Wing span [percent of span]','Interpreter','latex','FontSize',15);
ylabel('Lift coefficient $C_L$','Interpreter','latex','FontSize',15);
legend('$C_L$ for wich $C_{l,max}$ is reached at each station','$C_{L,max}$', 'Interpreter','latex')
grid on

% 3.4. Parabolic drag model curve and corresponding L/D curve:

% CD/CL:

CL_Curve = force_coeff(7, :);
CD_Curve = force_coeff(11, :);
RCD = polyfit(CL_Curve.^2, CD_Curve, 1);
CL_Plot = linspace(-0.8, 0.8, 100);
CD_Plot = RCD(1)*CL_Plot.^2 + RCD(2);

figure; hold on;
plot(CL_Curve, CD_Curve, 'or');
plot(CL_Plot, CD_Plot, 'b');
xlim([-0.8 0.8]);
xlabel('CL');
ylabel('CD');
title('Corba polar de resistència aerodinàmica de la ala');
legend('Valors obtinguts mitjançant el programa', 'Aproximació parabòlica de la corba polar', 'location', 'northeast');
grid on; hold off;

% L/D VS. Freestream velocity:

V_FS = linspace(10, 50, 100);
CL = 2*W_S*9.81./(1.225*V_FS.^2);       % Expressió derivada d'L = 1/2*S*CL*Rho*V^2, on L/S = W_S en [N/M^2]
CD = RCD(1)*CL.^2 + RCD(2);
EFF = CL./CD;

figure; hold on;
plot(V_FS, EFF, 'b');
xlabel('Velocitat de vol (M/S)');
ylabel('Eficiència aerodinàmica');
title('Corba L/D en funció de la velocitat de vol');
grid on; hold off;

% V -> (L/D)_MAX:

EFF_M = max(EFF);
V_EFF = V_FS(EFF==max(EFF));

% 3.5. CM_CG/CL for a series of flap deflections:

YFpos = [0.0 1.0];      % Full-span trailing edge flap
CFratio = 0.2;          % Indicated flap chord to chord ratio
FlapCorr = 1.0;         % Maximum flap effectiviness

XCG = 1.4;              % Indicated center of mass
ALPHA = -10:1:10;       % Studied attack angles
DEF = -20:5:20;         % Studied deflexions

CL_RANGE = sqrt(RCD(2)/RCD(1));             % CL at maximum-range flight conditions
CL_AUTONOMY = sqrt((3*RCD(2))/RCD(1));      % CL at minimum-descent-speed flight conditions

CL_Plot = -1.5:0.1:1.5;

figure; hold on;

for i = 1:length(DEF)

    [force_coeff, ~, mac, ~] = LLWing_function(AR, TR, DE25, ETIP, YF_pos, CF_ratio, DEF(i), FlapCorr, ALPHA, Cm0_tip, Cm0_root);

    CL_Curve = force_coeff(7, :);
    CM_Curve = force_coeff(5, :);
    RCM = polyfit(CL_Curve, CM_Curve, 1);

    S_W = mac*b;

    XAC = - RCM(1)*S_W;
    CM0 = - RCM(2)/RCM(1);
    CM_Plot = - CL_Plot*(XAC - XCG)/S_W + CM0;

    plot(CL_Plot, CM_Plot);

end

xline(CL_AUTONOMY, 'k--', 'C_{L} Màxima autonomia', 'LineWidth', 0.5);
xline(CL_RANGE, 'k--', 'C_{L} Màxim abast', 'LineWidth', 0.5);
title('Coeficient de moment (CG) en funció del coeficient de sustentació per a diferents \delta_e');
ylabel('CM_CG');
xlabel('CL');
legend([cellstr(num2str(DEF','\\delta_e = %3d'))], 'location', 'southwest');
grid on;

%% Functions

function [force_coeff, cl_local , mac, N] = LLWing_function(AR, TR, DE25, ETIP, YF_pos, CF_ratio, DE_flap, FlapCorr, ALPHA, Cm0_tip, Cm0_root)

    % Sections data (uses linear interpolation between root and tip)

    A0p = [ -0.9028 -0.9048 ]; % root and tip section zero-lift angles (deg)
    CM0p = [ Cm0_root Cm0_tip ]; % root and tip section free moments
    CDP = [ 0.007 -0.0035 0.0073  ;   % root section CD0, k1 and k2  (airfoil CD curve)
        0.0092 -0.0048 0.0088 ] ;  % tip section CD0, k1 and k2
    
    % Simulation data (by the time being only longitudinal analysis)

    N = 100 ; % number of panels along the span

    % LIFTING LINE SOLUTION %

    % Wing discretization (lenghts are dimensionless with the wing span)

    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr);
   
    % Assembly of the influence coefficients matrix (needed once)

    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h);

    % Solve circulations for different angles of attack

    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

    % Loads calculation using the Kutta-Joukowsky theorem (costlier, but
    % general... and illustrative!)

    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

end