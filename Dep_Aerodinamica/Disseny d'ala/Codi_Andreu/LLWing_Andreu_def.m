% -------------------------------------------------------------------------     
% PROGRAM LLWING: Weissinger lifting-line method for trapezoidal wings
% 220024 Aerodin√ mica ESEIAAT-UPC
% e.ortega@upc.edu
% -------------------------------------------------------------------------     

clc; clear all; close all;
format long;

% -------------------------------------------------------------------------
% POSTPROCESS ( your work starts here...)
% -------------------------------------------------------------------------

%% INPUTS
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
cr = 3 ; %Root chord
ct = 1.95 ; %Tip chord
b = 22; %wingspan
Sec = b*(cr+ct) / 2; %Wing section
Span=22;

AR = b^2/Sec ;   % aspect ratio
TR = ct/cr  ;   % taper ratio
xcg = 1;
rho = 1.225;
DE25 = 0 ; % sweep angle at c/4 (deg)
ETIP = -5; % tip twist (deg, negative for washout)

%% EFICIENCIA b, ct a la vegada

i = 1;
mat_eff = zeros(16,13);
for valor = 0.5:0.1:2
    
    ct = valor;
    j = 1;
    for valor2 = 15:1:27

        b = valor2;
        Sec = b*(cr+ct) / 2;   
        AR = b^2/Sec ;   
        TR = ct/cr  ;   

        [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

        [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

        [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;
    
        [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;

        mat_eff(i,j) = force_coeff (7 ,7)/force_coeff(11,7); % CANVIAR EL "6" SEGONS ANGLE D'ATAC DE LA CONDICI” DE VOL 

        j = j + 1;
        
    end    
   
    i = i + 1;
end

figure(1); 
hold on;
grid on;
for k = 1:size(mat_eff, 2)
    plot(0.5:0.1:2, mat_eff(:,k), 'DisplayName', ['b = ', num2str(k + 14)]);
end
title('CL/CD at 8 angle of attack', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
xlabel('ct', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
ylabel('CL/CD', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
legend('Interpreter', 'tex');
hold off;

%% INPUTS
A0p = [ -1.355 -1.36 ]; % root and tip section zero-lift angles (deg)
CM0p = [ -0.0025 -0.0025 ]; % root and tip section free moments
CDP = [ 0.00476039394846042 -0.00369115090706098 0.00632578060880848  ;  
        0.00435255288044264 -0.00315060448108891 0.00658385993719104 ] ;  % tip section CD0, k1 and k2

YF_pos = [ 0.13 0.64 ]; % 2y/b initial and final position of the flap/aileron in the half-wing
CF_ratio = 0.3 ;  % flap_chord/chord ratio
DE_flap = 10; % flap deflection (deg, positive:down)
FlapCorr = 0.9 ; % flap effectiviness (<=1)

N = 100 ; % number of panels along the span

ALPHA = [ -10. -8.0 -4.0 0. 4.0 8.0 10. 12 14] ; % angles of attack for analysis (deg)

cr = 2 ;
ct = 1.3 ;
b = 22 ;
Sec = b*(cr+ct) / 2 ;    % wing area
AR = b^2/Sec ;   % aspect ratio
TR = ct/cr  ;   % taper ratio
xcg = 1;
rho = 1.225;
DE25 = 0 ; % sweep angle at c/4 (deg)
ETIP = -5; % tip twist (deg, negative for washout)

%% TWIST

Clmax = 1.45; % angle stall reynolds 6.5*10^6 aprox
spanx = linspace(-0.5*b,0.5*b,N);

CL_stall = zeros(1, 7);
ALPHA_stall = zeros(2, 7);

cont = 1; 
for i = -7:2:5
    ETIP = i;
    ALPHA_stall(1,cont) = ETIP;

    [c4nods,c75nods,chord,s_pan,h,Cm0_y,normals,mac,S] = geo(AR,TR,N,DE25,ETIP,A0p,CM0p,CDP,YF_pos,CF_ratio,DE_flap,FlapCorr); 

    [inv_A,wake_len] = infcoeff(N,c4nods,c75nods,normals,h) ;

    [GAMMA,Ui,ncases] = getcirc(N,ALPHA,inv_A,normals) ;

    [cl_local,force_coeff] = KuttaJoukowsky(N,c4nods,h,GAMMA,Ui,s_pan,Cm0_y,chord,CDP,ncases,wake_len,S,mac,ALPHA) ;
    
    CL_alpha = polyfit(ALPHA,force_coeff(7,:),1);

    cl1 = cl_local(:, 5);
    CL1 = force_coeff(7, 5);
    cl2 = cl_local(:, 6);
    CL2 = force_coeff(7, 6);
    Cla = (cl1 - cl2) ./ (CL1 - CL2);
    Clb = cl1 - Cla * CL1;
    CL_span = (Clmax - Clb)./Cla;
    CL_stall(1,cont) = min(CL_span);

    ALPHA_stall(2,cont) = (CL_stall(1,cont) - CL_alpha(2))/CL_alpha(1);
    
    % Cl(y, epsilon)
    figure(2);
    grid on;
    hold on;
    plot(spanx, Clb + CL_stall(1,cont).*Cla, 'DisplayName', ['\epsilon = ', num2str(i), '^\circ']);
    legend_entries1{cont} = ['\epsilon = ', num2str(i), '^\circ'];
    title('C_l(\epsilon) distribution', 'FontSize',14, 'FontWeight', 'bold', 'Interpreter', 'tex');
    xlabel('y', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
    ylabel('C_l', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
    axis([-0.5*b 0.5*b 0 1.5]); 

    % CL(alpha, epsilon)
    figure(3);
    grid on;
    hold on;
    plot(ALPHA, force_coeff(7,:), 'DisplayName', ['\epsilon = ',num2str(i), '^\circ']);
    legend_entries2{cont} = ['\epsilon = ', num2str(i), '^\circ'];
    title('C_L(\alpha)', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
    xlabel('\alpha', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
    ylabel('C_L', 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex');
    
    cont = cont + 1;
end

figure(2);
yline(Clmax, '--k', 'LineWidth', 0.7, 'DisplayName', '');
legend(legend_entries1, 'Interpreter', 'tex');

figure(3);
hold on;
legend(legend_entries2, 'Interpreter', 'tex');
plot(ALPHA_stall(2,:), CL_stall(1,:), 'k', 'LineWidth',1.5, 'DisplayName', 'Stall');


figure(4);
plot(ALPHA_stall(1,:),ALPHA_stall(2,:));
title('\alpha_{stall} (\epsilon)');
xlabel('\epsilon', 'FontSize',14, 'FontWeight', 'bold');
ylabel('\alpha_{stall}', 'FontSize',14, 'FontWeight', 'bold');