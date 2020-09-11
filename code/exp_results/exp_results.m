clc
close all
clear variables

holding_p = -70; %mV
holding_t = 450; %ms
P1 = 50; %mV
P1_t = 25*1000; % ms
Ek = -91.1;

to_param_ko = [-3.10902175284165,-11.4945930320253,68.44919,37.5847380960043,29.9558997795475,0.170744990408810];
kslow1_param_ko = [41.6580470603940,5.47909573395254,31.3763688930612,9.45978973170458,1111.62338230588,0.0227149556759817];
kslow2_param_ko = [41.5408310196983,11.4792079724439,42.0528583047911,8.66747473211152,11263.9824825858,0.0265283219824804];
param_ko = [to_param_ko, kslow1_param_ko, kslow2_param_ko];

% [t1_ko, ~, A1_ko] = Ito(to_param_ko, holding_p, holding_t, P1, P1_t, Ek);
% [t2_ko, ~, A2_ko] = IKslow(kslow1_param_ko, holding_p, holding_t, P1, P1_t, Ek);
% [t3_ko, ~, A3_ko] = IKslow(kslow2_param_ko, holding_p, holding_t, P1, P1_t, Ek);

% Ito_ko = A1_ko(:,5);
% IKslow1_ko = A2_ko(:,5);
% IKslow2_ko = A3_ko(:,5);

to_param_wt = [36.0323917855634,41.0671435682954,54.3892960635912,53.4110530337884,34.0167375236982,0.224449754332619];
kslow1_param_wt = [17.1786204850059,3.55360447475750,32.6483194442754,1.25738530988108,1117.69914528303,0.123107346569954];
kslow2_param_wt = [-23.7191643973443,5.28511873412509,31.7611334578887,14.1815420635500,7231.18383776139,0.0511580416798235];
param_wt = [to_param_wt, kslow1_param_wt, kslow2_param_wt];

% [t1_wt, ~, A1_wt] = Ito(to_param_wt, holding_p, holding_t, P1, P1_t, Ek);
% [t2_wt, ~, A2_wt] = IKslow(kslow1_param_wt, holding_p, holding_t, P1, P1_t, Ek);
% [t3_wt, ~, A3_wt] = IKslow(kslow2_param_wt, holding_p, holding_t, P1, P1_t, Ek);

% Ito_wt = A1_wt(:,5);
% IKslow1_wt = A2_wt(:,5);
% IKslow2_wt = A3_wt(:,5);

%% Ito
V = -70:10:50;
numv = length(V);

to_peaks_ko = zeros(1, numv);
to_ssa_ko = zeros(1, numv);
to_ssi_ko = zeros(1, numv);
to_peaks_wt = zeros(1, numv);
to_ssa_wt = zeros(1, numv);
to_ssi_wt = zeros(1, numv);
for i=1:length(V)
    [tko, Sko, Ako] = Ito(to_param_ko, holding_p, holding_t, V(i), P1_t, Ek);
    [twt, Swt, Awt] = Ito(to_param_wt, holding_p, holding_t, V(i), P1_t, Ek);

    to_peaks_ko(i) = max(Ako(:,5));
    to_peaks_wt(i) = max(Awt(:,5));

    to_ssa_ko(i) = Ako(end,1) / (Ako(end,1) + Ako(end,2));
    to_ssi_ko(i) = Ako(end,3) / (Ako(end,3) + Ako(end,4));

    to_ssa_wt(i) = Awt(end,1) / (Awt(end,1) + Awt(end,2));
    to_ssi_wt(i) = Awt(end,3) / (Awt(end,3) + Awt(end,4));

    figure(1)
    hold on
        plot(tko, Ako(:,5))
    hold off
    title('I_{to} - KO')

    figure(2)
    hold on
        plot(twt, Awt(:,5))
    hold off
    title('I_{to} - WT')
end

figure(3)
hold on
    plot(V, to_peaks_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, to_peaks_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
legend('MGAT1KO', 'WT')
title('I_{to} - Densities')

figure(4)
hold on
    plot(V, to_ssa_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, to_ssa_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
legend('MGAT1KO', 'WT')
title('I_{to} - SSA')

% figure(5)
% hold on
%     plot(V, to_ssi_ko, '-o', 'Color','red')
%     plot(V, to_ssi_wt, '-o', 'Color','blue')
% hold off
% title('WT SS')

%% IKslow1
V = -70:10:50;
numv = length(V);

kslow1_peaks_ko = zeros(1, numv);
kslow1_ssa_ko = zeros(1, numv);
kslow1_ssi_ko = zeros(1, numv);
kslow1_peaks_wt = zeros(1, numv);
kslow1_ssa_wt = zeros(1, numv);
kslow1_ssi_wt = zeros(1, numv);
for i=1:length(V)
    [tko, Sko, Ako] = IKslow(kslow1_param_ko, holding_p, holding_t, V(i), P1_t, Ek);
    [twt, Swt, Awt] = IKslow(kslow1_param_wt, holding_p, holding_t, V(i), P1_t, Ek);

    kslow1_peaks_ko(i) = max(Ako(:,5));
    kslow1_peaks_wt(i) = max(Awt(:,5));

    kslow1_ssa_ko(i) = Ako(end,1);
    kslow1_ssi_ko(i) = Ako(end,2);

    kslow1_ssa_wt(i) = Awt(end,1);
    kslow1_ssi_wt(i) = Awt(end,2);

    figure(1)
    hold on
        plot(tko, Ako(:,5))
    hold off
    title('I_{Kslow1} - KO')

    figure(2)
    hold on
        plot(twt, Awt(:,5))
    hold off
    title('I_{Kslow1} - WT')
end

figure(3)
hold on
    plot(V, kslow1_peaks_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow1_peaks_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
legend('MGAT1KO', 'WT')
title('I_{Kslow1} - Densities')

figure(4)
hold on
    plot(V, kslow1_ssa_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow1_ssa_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
legend('MGAT1KO', 'WT')
title('I_{Kslow1} - SSA')

% figure(5)
% hold on
%     plot(V, kslow1_ssi_ko, '-o', 'Color','red')
%     plot(V, kslow1_ssi_wt, '-o', 'Color','blue')
% hold off
% title('WT SS')

%% IKslow2
V = -70:10:50;
numv = length(V);

kslow2_peaks_ko = zeros(1, numv);
kslow2_ssa_ko = zeros(1, numv);
kslow2_ssi_ko = zeros(1, numv);
kslow2_peaks_wt = zeros(1, numv);
kslow2_ssa_wt = zeros(1, numv);
kslow2_ssi_wt = zeros(1, numv);
for i=1:length(V)
    [tko, Sko, Ako] = IKslow(kslow2_param_ko, holding_p, holding_t, V(i), P1_t, Ek);
    [twt, Swt, Awt] = IKslow(kslow2_param_wt, holding_p, holding_t, V(i), P1_t, Ek);

    kslow2_peaks_ko(i) = max(Ako(:,5));
    kslow2_peaks_wt(i) = max(Awt(:,5));

    kslow2_ssa_ko(i) = Ako(end,1);
    kslow2_ssi_ko(i) = Ako(end,2);

    kslow2_ssa_wt(i) = Awt(end,1);
    kslow2_ssi_wt(i) = Awt(end,2);

    figure(1)
    hold on
        plot(tko, Ako(:,5))
    hold off
    title('I_{Kslow2} - KO')

    figure(2)
    hold on
        plot(twt, Awt(:,5))
    hold off
    title('I_{Kslow2} - WT')
end

figure(3)
hold on
    plot(V, kslow2_peaks_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow2_peaks_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
legend('MGAT1KO', 'WT')
title('I_{Kslow2} - Densities')

figure(4)
hold on
    plot(V, kslow2_ssa_ko, '-o', 'Color','red', 'LineWidth',2)
    plot(V, kslow2_ssa_wt, '--s', 'Color','blue', 'LineWidth',2)
hold off
legend('MGAT1KO', 'WT')
title('I_{Kslow2} - SSA')

% figure(5)
% hold on
%     plot(V, kslow2_ssi_ko, '-o', 'Color','red')
%     plot(V, kslow2_ssi_wt, '-o', 'Color','blue')
% hold off
% title('WT SS')
