clc
close all
clear variables

holding_p = -70; %mV
holding_t = 50; %ms
P1 = 50; %mV
P1_t = 1500; % ms
Ek = -91.1;
% X = [0.18064, 0.03577, 30.0, ...
%     0.3956, -0.06237, 30.0, ...
%     0.000152, 13.5, 7.0, 0.0067083, 33.5, 7.0, ...
%     0.00095, 33.5, 7.0, 0.051335, 33.5, 7.0, 0.4067];
X = [0.018064, 0.03577, 30.0, ...
    0.3956, -0.06237, 30.0, ...
    0.000152, 13.5, 7.0, 0.0067083, 33.5, 7.0, ...
    0.0000095, 33.5, 7.0, 0.051335, 33.5, 7.0, 0.4067];
set(0, 'DefaultAxesFontSize',14, 'DefaultAxesFontWeight','bold')


%% Ito
[t, S, A] = Ito(X, holding_p, holding_t, 50, P1_t, Ek);
Ito_trc = A(:,5);

plot(t, Ito_trc, 'LineWidth',2)
axis tight

% peak and tau
[peak, peak_idx] = max(Ito_trc);
trc_rd = Ito_trc(peak_idx:end);
tt_rd = t(peak_idx:end);
[delta, tau_idx] = min(abs(peak*exp(-1) - trc_rd));
tau = tt_rd(tau_idx);

% state variables
ato = S(:,1);
ito = S(:,2);

% close-up shot
% plot(t(86:195), Ito_trc(86:195), 'LineWidth',2)
% plot(t(86:195), ato(86:195), 'LineWidth',2)
% plot(t(1:800), ito(1:800), 'LineWidth',2)
% axis tight

% transition rates
ass = A(end,1) / (A(end,1) + A(end,2));
iss = A(end,3) / (A(end,3) + A(end,4));
tau_a = 1 / (A(end,1) + A(end,2));
tau_i = 1 / (A(end,3) + A(end,4));

fprintf('Peak: %f | Tau: %f | Peak time: %f| a_peak: %f | i_peak: %f | ass: %f | iss: %f | tau_a: %f | tau_i: %f \n', ...
    peak, tau, t(peak_idx), ato(peak_idx), ito(peak_idx), ass, iss, tau_a, tau_i)
fprintf('alpha_a: %f / beta_a: %f / alpha_i: %f / beta_i: %f \n', ...
    A(end,1), A(end,2), A(end,3), A(end,4))

figure(2)
plot(t, ato, 'LineWidth',2)
hold on
plot(t, ito, 'LineWidth',2)
hold off

% Ito & ito
figure(1)
plot(t, Ito_trc, 'LineWidth',2)
title('I_{to} under 40 mV')
axis tight

figure (2)
plot(t, ito, 'LineWidth',2)
title('i_{to} under 40 mV')
axis tight

% ato & ito
figure(2)
plot(t, ato, 'LineWidth',2)
hold on
plot(t, ito, 'LineWidth',2)
hold off
legend('ato', 'ito')

[peak, peak_idx] = max(Ito_trc);
ato(peak_idx) % ato
ito(peak_idx) % ito

% close-up shot of ato
figure(3)
plot(t(85:200), ato(85:200), 'LineWidth',2)
axis tight
title('Close-up shot of a_{to}')

% alpha_a and beta_a at 50 mV
figure(4)
plot(t, A(:,3), 'LineWidth',2)
hold on
plot(t, A(:,4), 'LineWidth',2)
hold off
legend('\alpha_{i}', '\beta_{i}')
axis tight

% alpha_a and beta_a over V
V = -80:1:60;
aa = X(1).*exp(X(2).*(V + X(3)));
ba = X(4).*exp(X(5).*(V + X(6)));
plot(V, aa, 'LineWidth',2)
hold on
plot(V, ba, 'LineWidth',2)
hold off
legend('\alpha_{a}', '\beta_{a}')

ai = (X(7).*exp( - (V+X(8))./X(9)))./( X(10).*exp( - (V+X(11))./X(12))+1);
bi = (X(13).*exp((V+X(14))./X(15)))./( X(16).*exp((V+X(17))./X(18))+1);
% ai = (X(7).*exp( - (V+X(8))./X(9)))./( X(10).*exp( - (V+X(11))./X(12)));
% bi = (X(13).*exp((V+X(14))./X(15)))./( X(16).*exp((V+X(17))./X(18)));

plot(V, ai, 'LineWidth',2)
hold on
plot(V, bi, 'LineWidth',2)
hold off
legend('\alpha_{i}','\beta_{i}')
axis tight

% V = -150:1:50;
ai1 = (X(7).*exp( - (V+33.5)./X(9)))./( X(10).*exp( - (V+X(11))./X(12))+1);
ai2 = (X(7).*exp( - (V+X(8))./X(9)))./( X(10).*exp( - (V+X(11))./X(12))+1);
plot(V, ai1, 'LineWidth',2)
hold on
plot(V, ai2, 'LineWidth',2)
hold off
legend('1', '2')

% ass & iss
ato = aa./(aa+ba);
ito = ai./(ai+bi); 

plot(V, ato, 'LineWidth',2)
hold on
plot(V, ito, 'LineWidth',2)
plot(V, ato+ito, '--', 'LineWidth',2)
hold off
legend('a_{ss}', 'i_{ss}', 'Sum')
axis tight

% tau_a & tau_i
tau_a = 1 ./ (aa+ba);
tau_i = 1 ./ (ai+bi);

plot(V, tau_a, 'LineWidth',2)
hold on
plot(V, tau_i, 'LineWidth',2)
hold off
legend('\tau_{a}','\tau_{i}')

plot(V, tau_a, 'LineWidth',2)
title('\tau_{a}')
axis tight

% current densities
P1s = -70:10:50;
peaks = zeros(1, length(P1s));
for i=1:length(P1s)
%     hold on
    [tt, SS, AA] = Ito(X, holding_p, holding_t, P1s(i), P1_t, Ek);
    trc = AA(:,5);
    [peak, peak_idx] = max(trc);
    peaks(i) = peak;
    
    trc_rd = trc(peak_idx:end);
    tt_rd = tt(peak_idx:end);
    [delta, tau_idx] = min(abs(peak*exp(-1) - trc_rd));
    tau = tt_rd(tau_idx);
    
    if (delta >= 0.5) || (trc_rd(end) > peak*exp(-1))
        disp('Tau out of time')
        tau = P1_t;
    end
    
    fprintf('Iter %i | Voltage: %i | Peak: %f | Tau %f \n', i, P1s(i), peak, tau)
    
%     plot(tt, AA(:,5))
%     hold off
end
max_peak = max(peaks);

plot(P1s, peaks/max_peak, '-o', 'LineWidth',2)
axis tight

% alpha_a and beta_a at 50 mV
figure(4)
plot(t, A(:,1), 'LineWidth',2)
hold on
plot(t, A(:,2), 'LineWidth',2)
hold off
legend('\alpha_{a}', '\beta_{a}')
axis tight

% alpha_a and beta_a over V
V = -70:1:50;
aa = X(1).*exp(X(2).*(V + X(3)));
ba = X(4).*exp(X(5).*(V + X(6)));
plot(V, aa, 'LineWidth',2)
hold on
plot(V, ba, 'LineWidth',2)
hold off
legend('\alpha_{a}', '\beta_{a}')
