clc
close all
clear variables


%% import relevant data
% DoE matrix
doe_mx = readtable('./FrF2_dgn_matrix.csv');

% trace
load('ds_Ktrace_ko.mat')
ds_Ktrace = ds_Ktrace_ko;
ds_Ktrace.Properties.VariableNames = {'time', 'I'};

% real experiment data
K_data = readtable('./potassium-KO.xlsx');

Iss = K_data.IssFF;
Iss = nanmean(Iss);

cap = K_data.CapFF;
cap = nanmean(cap);

% normalize the trace
ds_Ktrace.I = ds_Ktrace.I ./ cap;


%% DoE
holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
P2 = -70; % mV
P2_t = P1_t; % ms

X0 = [30.0, 30.0, 13.5, 7.0, 33.5, 33.5, 7.0, ...
      22.5, 7.7, 45.2, 5.7, 2.058, 1200.0, 45.2, 5.7, ...
      22.5, 7.7, 45.2, 5.7, 2.058, 1200.0, 45.2, 5.7];
low = [0.0, 0.0, 0.0, 2.0, 0.0, 20.0, 2.0, ...
       0.0, 2.0, 0.0, 2.0, 0.0, 170.0, 0.0, 2.0, ...
       0.0, 2.0, 0.0, 2.0, 0.0, 170.0, 0.0, 2.0];
high = [70.0, 70.0, 70.0, 50.0, 70.0, 70.0, 50.0, ...
        35.0, 14.0, 80.0, 24.0, 5.0, 5000.0, 70.0, 30.0, ...
        35.0, 14.0, 80.0, 24.0, 5.0, 5000.0, 70.0, 30.0,];

num_vars = length(X0);
num_runs = size(doe_mx);
num_runs = num_runs(1);
res = zeros(num_runs, 1); 
for i=1:num_runs
    fprintf("Exp %i/%i \n", i, num_runs)
    treat = table2array(doe_mx(i, :));
    param = zeros(1, num_vars);
    for j=1:num_vars
        if treat(j) == 1
            param(j) = high(j);
        else
            param(j) = low(j);
        end
    end
    
    [t, ~, A, ~] = Kv(param, holding_p, holding_t, P1, P1_t, P2, P2_t);
    IKsum = A(:,5) + A(:,10) + A(:,15);
    
    trace_sim = dtw(IKsum, ds_Ktrace.I);
    res(i) = trace_sim;
end

save('./doe_res.mat', 'res')