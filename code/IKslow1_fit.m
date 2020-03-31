clc
close all
clear variables


%% KO
K_ko = readtable("./potassium-KO.xlsx");
IKslow1_ko = K_ko.A2FF;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% X = [0.01176, -0.0631, 0.038198, -0.04178, 0.023391, 5.0, -0.03268, 5.0]
num_vars = 8;
lower_bd = [-0.01176, -0.1893, -0.038198, -0.12534, -0.023391, -5, -0.098, -5];
upper_bd = [0.03528, 0.0631, 0.114594, 0.04178, 0.070173, 15, 0.03268, 15];
fit_fn = @(X) IKslow1_fitness(X,IKslow1_ko,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 10);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow1 = A(:,67);
    
    % save the results
    file_path = sprintf('./GA_IKslow1_ko_%i.mat', i);
    rst(i,:) = [X, fval, fitted_IKslow1(end)];
    save(file_path, 'rst');
    disp(rst)
end


%% WT
K_wt = readtable("./potassium-WT.xlsx");
IKslow1_wt = K_wt.A2;

% voltage clamp protocol parameters
holding_p = -70; % mV
holding_t = 4.5; % sec
P1 = 50; % mV
P1_t = 29.5; % sec
P2 = 50; % mV
P2_t = 29.5; % sec

% GA parameters
% [0.01176, -0.0631, 0.038198, -0.04178, 0.023391, 5.0, -0.03268, 5.0]
num_vars = 8;
lower_bd = [-0.05, -0.1, -0.05, -0.05, -0.05, -5, -0.05, -5];
upper_bd = [0.05, 0.1, 0.05, 0.05, 0.05, 15, 0.05, 15];
fit_fn = @(X) IKslow1_fitness(X,IKslow1_wt,holding_p,holding_t,P1,P1_t,P2,P2_t);

% run GA
rst = zeros(10, 10);
for i=1:10
    fprintf('### Iter %i / 10', i)
    
    [X,fval] = ga(fit_fn,num_vars,[],[],[],[],lower_bd,upper_bd);
    
    % run the Rasmusson with the fitted parameters
    [~,~,A,~] = IKslow1(holding_p,holding_t,P1,P1_t,P2,P2_t,X);
    
    % fitted value 
    fitted_IKslow1 = A(:,67);
    
    % save the results
    file_path = sprintf('./GA_IKslow1_wt_%i.mat', i);
    rst(i,:) = [X, fval, fitted_IKslow1(end)];
    save(file_path, 'rst');
    disp(rst)
end
