clc
close all
clear variables
format long


%% data
cap = 254.3;
load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
% load('ds_Ktrace_wt.mat')
% ds_Ktrace = ds_Ktrace_wt;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};
ds_Ktrace.I = ds_Ktrace.I ./ cap;

[peak, peak_idx] = max(ds_Ktrace.I);
Ktrace_wo_ss = ds_Ktrace.I;
Ktrace_wo_ss(peak_idx:end) = Ktrace_wo_ss(peak_idx:end) - Iss_amp;

init_to = [-13.5655  128.4098  321.7877  127.2189   58.4796];
init_Kslow1 = [-0.0613    0.0097    0.2070    0.0128    1.1628];
init_Kslow1 = init_Kslow1*1000;
init_Kslow2 = [-0.0717    0.0123    0.0245    0.0399    8.6985];
init_Kslow2 = init_Kslow2*1000;
init_param = [init_to init_Kslow1 init_Kslow2];

% finetuning
% new_param = [-0.039528373401724   0.075202336459042   0.321175496958008   0.104292617589531   0.045989602921346    -0.187040228812866   0.007589258786537   0.217235124305737   0.011179502733295   1.096912251230504  -0.081358869811239   0.001300876833982   0.267302261070732   0.029954020941346   8.703977664513110];
% new_param = new_param*1000;

holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

[t, ~, A, ~] = IKsum(init_param, holding_p, holding_t, P1, P1_t, Ek);
IKsum_sim = A(:,5) + A(:,10) + A(:,15);

figure(1)
plot(ds_Ktrace.time, ds_Ktrace.I, 'LineWidth',2)
hold on
plot(t, IKsum_sim, 'LineWidth',2)
hold off
title('After Finetuning')
legend('Experimental','Simulated')

[peak, peak_idx] = max(IKsum_sim);
Iss_amp = 3.1;
IKsum_sim_ss = IKsum_sim;
IKsum_sim_ss((peak_idx+1):end) = IKsum_sim_ss((peak_idx+1):end) + Iss_amp;

plot(ds_Ktrace.time, ds_Ktrace.I, 'LineWidth',2)
hold on
plot(t, IKsum_sim_ss, 'LineWidth',2)
hold off
title('Whole K+ Current Trace - KO')
ylabel('pA/pF')
xlabel('Time(ms)')
legend('Experimental','Simulated')

dtw(ds_Ktrace.I, IKsum_sim_ss)

Iss = zeros(length(t), 1);
Iss(1:peak_idx-1) = 0;
Iss(peak_idx:end) = Iss_amp;
plot(t, A(:,5), 'LineWidth',2)
hold on
plot(t, A(:,10), 'LineWidth',2)
plot(t, A(:,15), 'LineWidth',2)
plot(t, Iss, 'LineWidth',2)
hold off
title('Key K+ Component Currents - KO')
ylabel('pA/pF')
xlabel('Time(ms)')
legend('Ito', 'IKslow1', 'IKslow2', 'Iss')


%% run AGA
% [best_fits, best_chrom] = IKsum_AGA(Ktrace_wo_ss, 15, 100, 30, 6, 4);
