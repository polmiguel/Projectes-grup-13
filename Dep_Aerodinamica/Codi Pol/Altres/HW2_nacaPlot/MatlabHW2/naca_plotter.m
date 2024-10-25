clear;
clf;
o = 3; %nombre de perfils que hi ha
cont = 0;
%assignin('base','naca_noviscositat',readmatrix(strcat("24112_SENSE_VISCOSITAT.csv")));
%assignin('base','naca_ambviscositat',readmatrix(strcat("24112_AMB_VISCOSITAT.csv")));
assignin('base','naca_re2300Mn4',readmatrix(strcat("Resultats_2_300M.csv")));
%assignin('base','naca_re3Mn9',readmatrix(strcat("24112_re3M_n9.csv")));
    
    for k = 3:3
       assignin('base','naca_r05',readmatrix(strcat("HW2_nacaPlot/MatlabHW2/Re_500/xf-naca2",int2str(4),"112-jf-500000.csv")));
       %assignin('base','naca_r10',readmatrix(strcat("Re_1000/xf-naca2",int2str(4),"112-jf-1000000.csv")));
       
       hold on;
       
       %cont = cont +7;
       txt = strcat ('NACA 2',int2str(4),'112 - 500.000') ;
       %plot(naca_r05(:,2),naca_r05(:,3),DisplayName=txt)%plot del angle i de cl
       
       
       
   end
   %plot(naca_noviscositat(:,1),naca_noviscositat(:,2),DisplayName='No viscositat')
   %plot(naca_ambviscositat(:,1),naca_ambviscositat(:,2),DisplayName='NACA 24112 - 3.000.000')
   %plot(naca_re2Mn4(:,1),naca_re2Mn4(:,2),DisplayName='re2Mn4')
   %plot(naca_re3Mn9(:,1),naca_re3Mn9(:,2),DisplayName='re3Mn9')
   legend show
   hold off;

   %figure;
   grid on;
   hold on;
   plot(naca_re2300Mn4(:,2),naca_re2300Mn4(:,5), DisplayName='2M')
   plot(naca_r05(:,2),naca_r05(:,5),DisplayName='500m')
   hold off;
   Cd_Root = polyfit(naca_re2300Mn4(:,2),naca_re2300Mn4(:,3),2)
   Cd_tip = polyfit(naca_r05(:,2),naca_r05(:,3),2)

   %si voleu plotejar mes coses us recomano que obriu l'arxiu csv
   % alla trobareu l'ordre del vector naca (naca(:,1) es l'angle, naca(:,2) es el cl,etc)
   % quan volgueu representar l'altre reynolds simplement poseu al plot
   % naca_r10 i aleshores s'us fara la representacio de les condicions en
   % Re = 1.000.000.
   % qualsevol cosa em podeu trucar,np 
   % P.D ignoreu lo de la matriu, es una merda i no funciona, feu nomes cas
   % als naca. 