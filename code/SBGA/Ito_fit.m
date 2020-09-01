clc
close all
clear variables
format shortG


%% Main
% holding_p = -70; %mV
% holding_t = 450; %ms
% P1 = 50; %mV
% P1_t = 25*1000; % ms
% Ek = -91.1;
% X0 = [30.0, 30.0, 13.5, 33.5, 7.0, 0.4067];
% [t, ~, A] = Ito(best_chroms(end,:), holding_p, holding_t, P1, P1_t, Ek);
% plot(t, A(:,5))

y = [17.6, 111.2];
tol = [0.2, 1.0];
[best_amps, best_taus, best_gens, best_chroms] = Ito_SBGA(y, tol, 30, 6, 4);

param = best_chroms(end,:);
holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;
[t, ~, A] = Ito(param, holding_p, holding_t, P1, P1_t, Ek);

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

save('Ito_KO')
