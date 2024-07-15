clc,close all,clear all

load('res_w_P_L2.mat')
load('res_w_P_L10.mat')
load('res_w_P_L50.mat')
load('res_w_P_L100.mat')
load('res_wo_P.mat')


fS = 20; lW = 3; 

x2_lim = 3*ones(1,size(tVec_P,2));

dx_lim = [-3;3]*ones(1,size(tVec_P,2));

u_lim = [-0.1;0.1]*ones(1,size(tVec_P,2));


x1 = [DATA_w_Param.x(1,:);DATA2_w_Param.x(1,:);DATA3_w_Param.x(1,:);DATA4_w_Param.x(1,:);DATA_wo_Param.x(1,:)];
x2 = [DATA_w_Param.x(2,:);DATA2_w_Param.x(2,:);DATA3_w_Param.x(2,:);DATA4_w_Param.x(2,:);DATA_wo_Param.x(2,:)];
x3 = [DATA_w_Param.x(3,:);DATA2_w_Param.x(3,:);DATA3_w_Param.x(3,:);DATA4_w_Param.x(3,:);DATA_wo_Param.x(3,:)];
x4 = [DATA_w_Param.x(4,:);DATA2_w_Param.x(4,:);DATA3_w_Param.x(4,:);DATA4_w_Param.x(4,:);DATA_wo_Param.x(4,:)];
x5 = [DATA_w_Param.x(5,:);DATA2_w_Param.x(5,:);DATA3_w_Param.x(5,:);DATA4_w_Param.x(5,:);DATA_wo_Param.x(5,:)];
x6 = [DATA_w_Param.x(6,:);DATA2_w_Param.x(6,:);DATA3_w_Param.x(6,:);DATA4_w_Param.x(6,:);DATA_wo_Param.x(6,:)];


u1 = [DATA_w_Param.u(1,:);DATA2_w_Param.u(1,:);DATA3_w_Param.u(1,:);DATA4_w_Param.u(1,:);DATA_wo_Param.u(1,:)];
u2 = [DATA_w_Param.u(2,:);DATA2_w_Param.u(2,:);DATA3_w_Param.u(2,:);DATA4_w_Param.u(2,:);DATA_wo_Param.u(2,:)];
u3 = [DATA_w_Param.u(3,:);DATA2_w_Param.u(3,:);DATA3_w_Param.u(3,:);DATA4_w_Param.u(3,:);DATA_wo_Param.u(3,:)];

tComp = [DATA_w_Param.tComp(1,:);DATA2_w_Param.tComp(1,:);DATA3_w_Param.tComp(1,:);DATA4_w_Param.tComp(1,:);DATA_wo_Param.tComp(1,:)];



figure(1)
subplot(3,1,1);plot(tVec_P,x1,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
set(gca,'fontweight','bold','fontsize', fS)
title('state position [m]')
ylabel('radial',"FontSize",fS);legend('off')
subplot(3,1,2);plot(tVec_P,x2,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,x2_lim,'--k','LineWidth',lW);
grid on;legend('off')
set(gca,'fontweight','bold','fontsize', fS)
ylabel('along track')
subplot(3,1,3);plot(tVec_P,x3,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
xlabel('time [s]');ylabel('cross track')
set(gca,'fontweight','bold','fontsize', fS)
legend('off')
figure(2)
subplot(3,1,1);plot(tVec_P,x4,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,dx_lim,'--k','LineWidth',lW);legend('off')
set(gca,'fontweight','bold','fontsize', fS)
title('state velocity [m/s]')
ylabel('radial')
subplot(3,1,2);plot(tVec_P,x5,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,dx_lim,'--k','LineWidth',lW);legend('off')
set(gca,'fontweight','bold','fontsize', fS)
ylabel('along track')
subplot(3,1,3);plot(tVec_P,x6,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,dx_lim,'--k','LineWidth',lW);legend('off')
set(gca,'fontweight','bold','fontsize', fS)
xlabel('time [s]');ylabel('cross track')

figure(3)
subplot(3,1,1);plot(tVec_P,u1,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,u_lim,'--k','LineWidth',lW);legend('off')
title('input [N/kg]')
set(gca,'fontweight','bold','fontsize', fS)
%title('input, red Line = without parametrization, blue Line = without parametrization  ')
ylabel('u1')
subplot(3,1,2);plot(tVec_P,u2,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,u_lim,'--k','LineWidth',lW);legend('off')
set(gca,'fontweight','bold','fontsize', fS)
ylabel('u2')
subplot(3,1,3);plot(tVec_P,u3,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
hold on; plot(tVec_P,u_lim,'--k','LineWidth',lW);
set(gca,'fontweight','bold','fontsize', fS)
xlabel('time [s]');ylabel('u3')
legend('off')


% figure(4)
% %title('Cost Function, red Line = without parametrization, blue Line = without parametrization  ')
% plot(tVec_P, J_w_P,tVec2_P, J2_w_P,'linewidth',lW);grid on;
% set(gca,'fontweight','bold','fontsize', fS)
% xlabel('time [s]');ylabel('J')


figure(5)
%title('computational time')
plot(tVec_P*2, tComp,'linewidth',lW);grid on; legend('pRSHMPC_1','pRSHMPC_2','pRSHMPC_3','pRSHMPC_4','RSHMPC','Orientation','horizontal','Location','best')
set(gca,'fontweight','bold','fontsize', fS)
xlabel('time-step k');ylabel('comp time[s]'); legend('off')


for i = 1:size(tComp,1)
    avg_tComp(i) = mean(tComp(i,:));
    sum_tComp(i) = sum(tComp(i,:));
end


avg_tComp
sum_tComp

cost_fcn_res

conePrint




