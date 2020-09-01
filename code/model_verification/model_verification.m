clc
close all
clear variables
format shortG


%% KO - expoenetial fn vs. my model
% fitted on 2020-08-28 by SGBA
to_param = [-3.10902175284165,-11.4945930320253,69.4491868862909,37.5847380960043,29.9558997795475,0.170744990408810];
kslow1_param = [-35.4541893813011,5.38117058719360,104.890059050153,2.51204645753509,1174.60499208924,0.0324191306499981];
kslow2_param = [41.5408310196983,11.4792079724439,42.0528583047911,8.66747473211152,11263.9824825858,0.0265283219824804];

holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

[t1, ~, A1] = Ito(to_param, holding_p, holding_t, P1, P1_t, Ek);
[t2, ~, A2] = IKslow(kslow1_param, holding_p, holding_t, P1, P1_t, Ek);
[t3, ~, A3] = IKslow(kslow2_param, holding_p, holding_t, P1, P1_t, Ek);

% Ito = A1(:,5);
% IKslow1 = A2(:,5);
% IKslow2 = A3(:,5);

% figure(1)
% plot(t1, Ito)
% hold on
% plot(t2, IKslow1)
% plot(t3, IKslow2)
% hold off
% axis tight

amps = [17.6, 3.1, 3.6];
taus = [111.2, 1115.1, 11266.1];

t_holding = 0:0.1:holding_t;
t_P1 = 0:0.1:(P1_t-holding_t);
t_P1_shift = t_P1 + holding_t + 0.1;
t = [t_holding, t_P1_shift];
holding_exp = zeros(1, length(t_holding)); 

Ito_exp = exp_fn(t_P1, amps(1), taus(1));
Ito_exp = [holding_exp, Ito_exp];

IKslow1_exp = exp_fn(t_P1, amps(2), taus(2));
IKslow1_exp = [holding_exp, IKslow1_exp];

IKslow2_exp = exp_fn(t_P1, amps(3), taus(3));
IKslow2_exp = [holding_exp, IKslow2_exp];

% figure(2)
% plot(t, Ito_exp)
% hold on
% plot(t, IKslow1_exp)
% plot(t, IKslow2_exp)
% hold off
% axis tight

figure(3)
subplot(1,3,1)
plot(t1, A1(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(1,3,2)
plot(t2, A2(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow1_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(1,3,3)
plot(t3, A3(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow2_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')


%% KO - prediction


%% WT - expoenetial fn vs. my model
% fitted on 2020-08-28 by SGBA
to_param_wt = [36.0323917855634,41.0671435682954,54.3892960635912,53.4110530337884,34.0167375236982,0.224449754332619];
kslow1_param_wt = [17.1786204850059,3.55360447475750,32.6483194442754,1.25738530988108,1117.69914528303,0.123107346569954];
kslow2_param_wt = [-23.7191643973443,5.28511873412509,31.7611334578887,14.1815420635500,7231.18383776139,0.0511580416798235];
param_wt = [to_param_wt, kslow1_param_wt, kslow2_param_wt];

[t, S, A, ~] = AP_param2(param_wt, 70);
plot(t, S(:,1))
axis tight

plot(t, A(:,76))
% 
% [tt, SS, AA, ~] = Rasmusson_AP(70);
% plot(tt, AA(:,8))

holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

[t1_wt, ~, A1_wt] = Ito(to_param_wt, holding_p, holding_t, P1, P1_t, Ek);
[t2_wt, ~, A2_wt] = IKslow(kslow1_param_wt, holding_p, holding_t, P1, P1_t, Ek);
[t3_wt, ~, A3_wt] = IKslow(kslow2_param_wt, holding_p, holding_t, P1, P1_t, Ek);

Ito_wt = A1_wt(:,5);
IKslow1_wt = A2_wt(:,5);
IKslow2_wt = A3_wt(:,5);

amps_wt = [24.8, 17.1, 7.3];
taus_wt = [105.2, 1119.6, 7271.9];

t_holding = 0:0.1:holding_t;
t_P1 = 0:0.1:(P1_t-holding_t);
t_P1_shift = t_P1 + holding_t + 0.1;
t = [t_holding, t_P1_shift];
holding_exp = zeros(1, length(t_holding)); 

Ito_exp_wt = exp_fn(t_P1, amps_wt(1), taus_wt(1));
Ito_exp_wt = [holding_exp, Ito_exp_wt];

IKslow1_exp_wt = exp_fn(t_P1, amps_wt(2), taus_wt(2));
IKslow1_exp_wt = [holding_exp, IKslow1_exp_wt];

IKslow2_exp_wt = exp_fn(t_P1, amps_wt(3), taus_wt(3));
IKslow2_exp_wt = [holding_exp, IKslow2_exp_wt];

figure
title('WT')
subplot(1,3,1)
plot(t1_wt, Ito_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
title('I_{to} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(1,3,2)
plot(t2_wt, IKslow1_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow1_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow1} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(1,3,3)
plot(t3_wt, IKslow2_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow2_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow2} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')


%%
figure
subplot(2,3,1)
plot(t1, A1(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
title('I_{to} - KO')
legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,3,2)
plot(t2, A2(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow1_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow1} - KO')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,3,3)
plot(t3, A3(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow2_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow2} - KO')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,3,4)
plot(t1_wt, Ito_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
title('I_{to} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,3,5)
plot(t2_wt, IKslow1_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow1_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow1} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

subplot(2,3,6)
plot(t3_wt, IKslow2_wt, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow2_exp_wt, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
% legend('Simulation Model','Exponential Function')
title('I_{Kslow2} - WT')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')