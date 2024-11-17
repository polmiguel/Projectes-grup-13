clear;
clc;

%assignin('base','naca_noviscositat',readmatrix(strcat("24112_SENSE_VISCOSITAT.csv")));
%assignin('base','naca_ambviscositat',readmatrix(strcat("24112_AMB_VISCOSITAT.csv")));
%assignin('base','naca_re2300Mn4',readmatrix(strcat("Resultats_2_300M.csv")));
%assignin('base','naca_re3Mn9',readmatrix(strcat("24112_re3M_n9.csv")));
assignin('base','ReRoot',readmatrix(strcat("T1_Re15.304_M0.36_N9.0.csv")));
assignin('base','ReTip',readmatrix(strcat("T1_Re7.650_M0.36_N9.0.csv")));
           
       
   %plot(naca_noviscositat(:,1),naca_noviscositat(:,2),DisplayName='No viscositat')
   %plot(naca_ambviscositat(:,1),naca_ambviscositat(:,2),DisplayName='NACA 24112 - 3.000.000')
   %plot(naca_re2Mn4(:,1),naca_re2Mn4(:,2),DisplayName='re2Mn4')
   %plot(naca_re3Mn9(:,1),naca_re3Mn9(:,2),DisplayName='re3Mn9')
   

   figure;
   grid on;
   hold on;
   plot(ReRoot(:,2),ReRoot(:,5), DisplayName='Root')
   plot(ReTip(:,2),ReTip(:,5),DisplayName='Tip')
   hold off;
   Cd_Root = polyfit(ReRoot(:,2),ReRoot(:,3),2)
   Cd_tip = polyfit(ReTip(:,2),ReTip(:,3),2)

