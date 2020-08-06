clc
close all
clear variables
format shortG


%% import taget data
cap_ko = 254.3;

raw_trc = readtable('k_trace.csv');
raw_trc = raw_trc(1:500001,:);

raw_ko = raw_trc.KO;
raw_ko = raw_ko ./ cap_ko;
raw_time = raw_trc.time;

% chop off early phase
[~, peak_idx] = max(raw_ko);
raw_ko_rd = raw_ko(peak_idx:end);
raw_ko_rd(raw_ko_rd < 0) = 0;
raw_time_rd = raw_time(peak_idx:end);
raw_time_rd = raw_time_rd - raw_time_rd(1);

raw_trc_rd = table(raw_time_rd, raw_ko_rd);
raw_trc_rd.Properties.VariableNames = {'time', 'current'};
ds_trc = downsample(raw_trc_rd, 400);

% peak and tau (of IKsum for initialization)
[peak, peak_idx] = max(ds_trc.current);
[~, tau_idx] = min(abs(peak*exp(-1) - ds_trc.current(peak_idx:end)));
tau = ds_trc.time(peak_idx+tau_idx);

%% run AGA
[best_amps, best_taus, best_gens, best_chroms] = direct_decomp(6, 6, [peak, tau], 30, 6, 4);


%% 


