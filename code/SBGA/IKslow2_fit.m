clc
close all
clear variables
format shortG


%% main
% holding_p = -70; %mV
% holding_t = 450; %ms
% P1 = 50; %mV
% P1_t = 25*1000; % ms
% Ek = -91.1;
% X0 = [22.5, 7.7, 45.2, 5.7, 1200.0, 0.16];
% 
% [t, ~, A] = IKslow(X0, holding_p, holding_t, P1, P1_t, Ek);
% plot(t, A(:,5))

y = [3.6, 11266.1];
tol = [0.2, 10.0];
[best_amps, best_taus, best_gens, best_chroms] = IKslow_SBGA(y, tol, 30, 6, 4);

param = best_chroms(end,:);
holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;
[t, ~, A] = IKslow(param, holding_p, holding_t, P1, P1_t, Ek);

t_holding = 0:0.1:holding_t;
t_P1 = 0:0.1:(P1_t-holding_t);
t_P1_shift = t_P1 + holding_t + 0.1;
tt = [t_holding, t_P1_shift];
holding_exp = zeros(1, length(t_holding)); 

Iexp = exp_fn(t_P1, y(1), y(2));
Iexp = [holding_exp, Iexp];

plot(t, A(:,5), 'LineWidth',2, 'Color','red')
hold on
plot(tt, Iexp, '--', 'LineWidth',2, 'Color','black')
hold off
axis tight
legend('Simulation Model','Exponential Function')
ylabel('Current (pA/pF)')
xlabel('Time (ms)')

% save('IKslow2_fit')
