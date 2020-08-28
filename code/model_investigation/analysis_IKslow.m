clc
close all
clear variables

holding_p = -70; %mV
holding_t = 50; %ms
P1 = 50; %mV
P1_t = 10000; % ms
Ek = -91.1;
% X0 = [22.5, 7.7, 45.2, 5.7, ...
%      2.058, 1200.0, 170.0, 0.16];
X = [22.5, 7.7, 45.2, 5.7, ...
    0.058, 1200.0, 170.0, 0.16];
set(0, 'DefaultAxesFontSize',14, 'DefaultAxesFontWeight','bold')


%% IKslow
[t, S, A] = IKslow(X, holding_p, holding_t, 50, P1_t, Ek);
IKslow_trc = A(:,5);

% peak and tau
[peak, peak_idx] = max(IKslow_trc);
trc_rd = IKslow_trc(peak_idx:end);
tt_rd = t(peak_idx:end);
[delta, tau_idx] = min(abs(peak*exp(-1) - trc_rd));
tau = tt_rd(tau_idx);

% state variables
aur = S(:,1);
iur = S(:,2);

% key statistics
ass = A(end,1);
iss = A(end,2);
tau_a = A(end,3);
tau_i = A(end,4);

% close-up shot
plot(t(45:170), IKslow_trc(45:170), 'LineWidth',2)
% plot(t(45:170), aur(45:170), 'LineWidth',2)
% plot(t(1:3000), iur(1:3000), 'LineWidth',2)
axis tight

fprintf('Peak: %f | Tau: %f | Peak time: %f| a_peak: %f | i_peak: %f | ass: %f | iss: %f | tau_a: %f | tau_i: %f \n', ...
    peak, tau, t(peak_idx), aur(peak_idx), iur(peak_idx), ass, iss, tau_a, tau_i)

plot(t, IKslow_trc, 'LineWidth',2)
axis tight

plot(t, aur, 'LineWidth',2)
hold on
plot(t, iur, 'LineWidth',2)
hold off
