close all
clear
clc
load('IJADE_CURVE.mat');
total2(:,:,1)=SAVE_PSO_mean_cg_curve;
load('SHADE_CURVE.mat')
total2(:,:,2)=SAVE_PSO_mean_cg_curve;
load('JADE_CURVE.mat')
total2(:,:,3)=SAVE_PSO_mean_cg_curve;

a=linspace(1,300000,600);
b=linspace(1,300000,600);
c=linspace(1,300000,600);
t=1;
for i=[1:30]%[3 8 17 26 ]
    if i==1
        figure(1);
        t=1;
    elseif i==16
        figure(2);
        t=1;
    end
    subplot(3,5,t)
    tit=['F',num2str(i)];%Convergence curve
   % colororder(newcolors);
    semilogy(a,total2(i,:,1),'r','LineWidth',2);
    hold on
    semilogy(b,total2(i,:,2),...
        c,total2(i,:,3),...
        'LineWidth',1.5);
%     title(tit,'Fontname', 'Times New Roman','FontSize',10)
    xlabel('Function Evaluations','Fontname', 'Times New Roman','FontSize',10,'Color',[0 0 0]);
    ylabel([tit,' Mean error f(x)-f(x*)'],'Fontname', 'Times New Roman','FontSize',10,'Color',[0 0 0]);%Times New RomanMicrosft YaHei UI
     axis auto
     t=t+1;
     if mod(t,6)==0
         legend('IJADE','JADE','SHADE','Location','southoutside','NumColumns',7,'Box','off','Fontname', 'Times New Roman','FontSize',10);
     end
end
