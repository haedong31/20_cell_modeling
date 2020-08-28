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

[best_amps, best_taus, best_gens, best_chroms] = IKslow_SBGA([], [], 30, 6, 4);
save('IKslow2_fit')
