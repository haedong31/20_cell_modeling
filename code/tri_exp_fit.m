clc
close all
clear variables
format shortG


%% import data
cap_ko = 254.3;

raw_trc = readtable('k_trace.csv');
raw_trc = raw_trc(1:500001,:);

raw_ko = raw_trc.KO;
raw_ko = raw_ko ./ cap_ko;
raw_time = raw_trc.time;

% chop off early phase
[peak, peak_idx] = max(raw_ko);
raw_ko_rd = raw_ko(peak_idx:end);
raw_ko_rd(raw_ko_rd < 0) = 0;
raw_time_rd = raw_time(peak_idx:end);
raw_time_rd = raw_time_rd - raw_time_rd(1);

raw_trc_rd = table(raw_time_rd, raw_ko_rd);
raw_trc_rd.Properties.VariableNames = {'time', 'current'};
ds_trc = downsample(raw_trc_rd, 400);


%% run AGA
[best_fits, best_gens, best_chroms] = tri_exp(7, ds_trc, 30, 6, 4);


%% modeling verification
t = ds_trc.time;
param = best_chroms(end,:);
Iss = param(1);
IKslow2 = exp_fn(t, param(2), param(5));
IKslow1 = exp_fn(t, param(3), param(6));
Ito = exp_fn(t, param(4), param(7));

IKsum = Iss + IKslow2 + IKslow1 + Ito;

figure(1)
plot(t, ds_trc.current, 'LineWidth',2, 'Color','red')
hold on
plot(t, IKsum, '--', 'LineWidth',2, 'Color','black')
hold off
legend('Read Data', 'Fitted Model')
axis tight

figure(2)
x = ones(length(t), 1)*Iss;
plot(t, ds_trc.current, 'LineWidth',2, 'Color','red')
hold on
plot(t, Ito, '--', 'LineWidth',2, 'Color','black')
plot(t, IKslow1, ':', 'LineWidth',2, 'Color','black')
plot(t, IKslow2, '-.', 'LineWidth',2, 'Color','black')
plot(t, x, 'LineWidth',2, 'Color','black')
hold off
legend('Real Data', 'Ito', 'IKslow1', 'IKslow2', 'Iss')
axis tight
