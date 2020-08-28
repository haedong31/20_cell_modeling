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

[best_amps, best_taus, best_gens, best_chroms] = Ito_SBGA[], [], 30, 6, 4);
save('Ito_fit')
