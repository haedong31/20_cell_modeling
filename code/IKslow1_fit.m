clc
close all
clear variables
format long


holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

% X0 = [22.5, 7.7, 45.2, 5.7, 1200.0, 0.16];
% [t, ~, A] = IKslow(X0, holding_p, holding_t, P1, P1_t, Ek);
% plot(t, A(:,5))

y = [3.1, 1115.1];
[best_amps, best_taus, best_gens, best_chroms] = IKslow1_AGA(6, y, 30, 6, 4);

[t1, ~, A1] = IKslow(param_ko(1,:), holding_p, holding_t, P1, P1_t, Ek);
IKslow1 = A1(:,5);
[t2, ~, A2] = IKslow(param_ko(2,:), holding_p, holding_t, P1, P1_t, Ek);
IKslow2 = A2(:,5);

plot(t1, IKslow1)
hold on
plot(t2, IKslow2)
hold off
legend('IKslow1','IKslow2')
