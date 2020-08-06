clc
close all
clear variables
format shortG


%% KO - expoenetial fn vs. my model
to_param = [22.8776463721699,6.84699170187748,5.33698464484186,71.3850538656544,53.2735458706752,0.0668657571418254];
kslow1_param = [-67.7036577810262,10.1469012323349,182.879818956047,11.8760185824273,1415.63825792747,0.333475860087409];
kslow2_param = [-77.5193722061745,20.3845417338815,18.2043800464169,2.27717832845248,11595.3463118204,0.126565848429251];

holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

[t1, ~, A1] = Ito(to_param, holding_p, holding_t, P1, P1_t, Ek);
[t2, ~, A2] = IKslow(kslow1_param, holding_p, holding_t, P1, P1_t, Ek);
[t3, ~, A3] = IKslow(kslow2_param, holding_p, holding_t, P1, P1_t, Ek);

Ito = A1(:,5);
IKslow1 = A2(:,5);
IKslow2 = A3(:,5);

figure(1)
plot(t1, Ito)
hold on
plot(t2, IKslow1)
plot(t3, IKslow2)
hold off
axis tight

amps = [9.03951859907668, 5.02749966708095, 3.67297141335227];
taus = [144.007385253906, 1366.2353515625, 11548.158203125];

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

figure(2)
plot(t, Ito_exp)
hold on
plot(t, IKslow1_exp)
plot(t, IKslow2_exp)
hold off
axis tight

figure(3)
plot(t1, Ito, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')

figure(4)
plot(t2, IKslow1, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow1_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')

figure(5)
plot(t3, IKslow2, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow2_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')


%% WT - expoenetial fn vs. my model
to_param = [29.8745828107163,432.363185871035,131.363384854194,2.33532115101248,14.2778720091296,0.0970723758881461];
kslow1_param = [0.753290839758920,24.9516633577433,73.3435975654207,2.59579933915273,1057.01873506607,0.114673726146839];
kslow2_param = [-56.5728114280857,16.8944598117128,65.5830552797651,9.54739616634568,5534.23439976122,0.128513753910693];

holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

[t1, ~, A1] = Ito(to_param, holding_p, holding_t, P1, P1_t, Ek);
[t2, ~, A2] = IKslow(kslow1_param, holding_p, holding_t, P1, P1_t, Ek);
[t3, ~, A3] = IKslow(kslow2_param, holding_p, holding_t, P1, P1_t, Ek);

Ito = A1(:,5);
IKslow1 = A2(:,5);
IKslow2 = A3(:,5);

amps = [14.1040588378906, 9.57069625854487, 7.0088455200195];
taus = [91.0120239257812, 989.449768066406, 5482.32666015625];

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

figure(1)
plot(t1, Ito, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')

figure(2)
plot(t2, IKslow1, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow1_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')

figure(3)
plot(t3, IKslow2, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKslow2_exp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
