clc
close all
clear variables
format shortG

holding_p = -70; %mV
holding_t = 50; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

set(0, 'DefaultAxesFontSize',14, 'DefaultAxesFontWeight','bold')


%% average IKsum data of KO group
cap_ko = 254.3;
raw_trc = readtable('k_trace.csv');

% remove tail after 25 sec
raw_trc = raw_trc(1:500001,:);

% pull out KO trace
raw_ko = raw_trc.KO;
raw_ko = raw_ko ./ cap_ko;
raw_ko(raw_ko < 0) = 0;
raw_time = raw_trc.time;

% chop off early phase
[~, peak_idx] = max(raw_ko);
raw_ko_rd = raw_ko(peak_idx:end);
raw_time_rd = raw_time(peak_idx:end);
raw_time_rd = raw_time_rd - raw_time_rd(1);

raw_trc_rd = table(raw_time_rd, raw_ko_rd);
raw_trc_rd.Properties.VariableNames = {'time', 'current'};
ds_trc = downsample(raw_trc_rd, 400);


%% manually fit IKsum
% X(1) = 0.18064; % alpha_a
% X(2) = 0.00095; % beta_i
% X(3) = 22.5; % ass
% X(4) = 1200; % tau_i
% X(5) = 22.5; % ass
% X(6) = 1200; % tau_i
% X(7:9) = [0.4067, 0.16, 0.16];

X(1) = 0.18064; % alpha_a
X(2) = 0.00085; % beta_i
X(3) = 22.5; % ass
X(4) = 1065.1; % tau_i
X(5) = 22.5; % ass
X(6) = 11266.1; % tau_i
X(7:9) = [0.1548, 0.0223, 0.0256];
Iss = 2.88;

% fit IKsum
[t, S, A] = IKsum(X, holding_p, holding_t, 50, P1_t, Ek);

Ito_trc = A(:,5);
IKslow1_trc = A(:,10);
IKslow2_trc = A(:,15);

sub_IKsum_trc = Ito_trc + IKslow1_trc + IKslow2_trc;

[~, peak_idx] = max(sub_IKsum_trc);
sub_IKsum_trc_rd = sub_IKsum_trc(peak_idx:end);
t_rd = t(peak_idx:end);

IKsum_trc_rd = sub_IKsum_trc_rd + Iss;

% statistics (peaks and taus)
[to_peak, to_pidx] = max(Ito_trc);
[Kslow1_peak, Kslow1_pidx] = max(IKslow1_trc);
[Kslow2_peak, Kslow2_pidx] = max(IKslow2_trc);

tto_rd = t(to_pidx:end);
tKslow1_rd = t(Kslow1_pidx:end);
tKslow2_rd = t(Kslow2_pidx:end);

[~, to_tau_idx] = min(abs(to_peak*exp(-1) - Ito_trc(to_pidx:end)));
[~, Kslow1_tau_idx] = min(abs(Kslow1_peak*exp(-1) - IKslow1_trc(Kslow1_pidx:end)));
[~, Kslow2_tau_idx] = min(abs(Kslow2_peak*exp(-1) - IKslow2_trc(Kslow2_pidx:end)));

to_tau = tto_rd(to_tau_idx);
Kslow1_tau = tKslow1_rd(Kslow1_tau_idx);
Kslow2_tau = tKslow2_rd(Kslow2_tau_idx);

fprintf('Ito - Peak: %f | Tau: %f \n', to_peak, to_tau)
fprintf('IKslow1 - Peak: %f | Tau: %f \n', Kslow1_peak, Kslow1_tau)
fprintf('IKslow2 - Peak: %f | Tau: %f \n', Kslow2_peak, Kslow2_tau)

figure(1)
plot(ds_trc.time, ds_trc.current, 'LineWidth',2, 'Color','red')
hold on
plot(t_rd, IKsum_trc_rd, '--', 'LineWidth',2, 'Color','black')
hold off
legend('Real Data','Simulated')
axis tight

figure(2)
plot(t, Ito_trc, 'LineWidth',2)
hold on
plot(t, IKslow1_trc, 'LineWidth',2)
plot(t, IKslow2_trc, 'LineWidth',2)
hold 
legend('Ito', 'IKslow1', 'IKslow2')
axis tight

