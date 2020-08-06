clc
close all
clear variables


%% Main
% holding_p = -70; %mV
% holding_t = 450; %ms
% P1 = 50; %mV
% P1_t = 25*1000; % ms
% Ek = -91.1;
% X0 = [30.0, 30.0, 13.5, 33.5, 7.0, 0.4067];
% [t, ~, A] = Ito(best_chroms(end,:), holding_p, holding_t, P1, P1_t, Ek);
% plot(t, A(:,5))

df = readtable("./potassium-KO.xlsx");
num_obs = 35;
num_var = 6;
amp = df.A3;
tau = df.Tau3;

best_amp_container = zeros(1, num_obs);
best_tau_container = zeros(1, num_obs);
num_iters_container = zeros(1, num_obs);
best_chrom_container = zeros(num_obs, num_var);

fprintf('### Iter 1/%i \n', num_obs)
y = [amp(1), tau(1)];
[best_amps, best_taus, best_gens, best_chroms] = Ito_AGA(num_var, y, 30, 6, 4);
best_amp_container(1) = best_amps(end);
best_tau_container(1) = best_taus(end);
num_iters_container(1) = best_gens(end);
best_chrom_container(1,:) = best_chroms(end,:);

for i=2:num_obs
    fprintf('### Iter %i/%i \n', i, num_obs)

    y = [amp(i), tau(i)];
    [best_amps, best_taus, best_gens, best_chroms] = Ito_AGA_seq(num_var, y, best_chrom_container(i-1,:), 30, 6, 4);
    best_amp_container(i) = best_amps(end);
    best_tau_container(i) = best_taus(end);
    num_iters_container(i) = best_gens(end);
    best_chrom_container(i,:) = best_chroms(end,:);

    save('Ito_fit')
end
