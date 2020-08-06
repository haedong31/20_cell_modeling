function [best_amps, best_taus, best_gens, best_chroms] = direct_decomp(nv_to, nv_Kslow, y, n0, n1, n2)
    global num_var;
    global N0;
    global N1;
    global N2;
    global sig_ctl;
    num_var = nv_to + nv_Kslow*2 + 1;
    N0 = n0;
    N1 = n1;
    N2 = n2;
    sig_ctl = [];
    best_fits = [];
    best_amps = [];
    best_taus = [];
    best_gens = [];
    best_chroms = [];
    
    amp_tol = 0.5;
    tau_tol = 0.5;

    % initialization
    low_to = [0.0, 0.0, 0.0, 20.0, 2.0, 0.2];
    high_to = [70.0, 70.0, 70.0, 70.0, 50.0, 0.6];
    low_Kslow = [-60.0, 2.0, -20.0, 2.0, 200.0, 0.05];
    high_Kslow = [80.0, 14.0, 80.0, 24.0, 5000.0, 0.3];

    % Ito initialization
    init_gen_to = zeros(N0, nv_to);
    for j=1:nv_to
        unif1 = makedist('Uniform', 'lower',low_to(j), 'upper',high_to(j));
        init_gen_to(:,j) = random(unif1, N0, 1);
    end
    
    % IKslow1 & IKslow2 initialization
    init_gen_Kslow1 = zeros(N0, nv_Kslow);
    init_gen_Kslow2 = zeros(N0, nv_Kslow);
    for j=1:nv_Kslow
        unif2 = makedist('Uniform', 'lower',low_Kslow(j), 'upper',high_Kslow(j));
        unif3 = makedist('Uniform', 'lower',low_Kslow(j), 'upper',high_Kslow(j));
        init_gen_Kslow1(:,j) = random(unif2, N0, 1);
        init_gen_Kslow2(:,j) = random(unif3, N0, 1);
    end
    
    % Iss initialization
    unif4 = makedist('Uniform', 'lower',2.08, 'upper',13.50);
    init_ss = random(unif4, N0, 1);
    
    % consolidated initial population
    init_gen = [init_gen_to, init_gen_Kslow1, init_gen_Kslow2, init_ss];
    
    cnt = 1;
    [fits, amp_dels, tau_dels] = eval_fn(init_gen, y);
    [bf, bf_idx] = min(fits);
    bamp = amp_dels(bf_idx);
    btau = tau_dels(bf_idx);
    bchrom = init_gen(bf_idx,:);

    fprintf('Initial best fit: %f|Amp: %f|Tau: %f \n', bf, amp_dels(bf_idx), tau_dels(bf_idx));
    disp(bchrom)

    best_cnt = 1;
    best_fits = [best_fits, bf];
    best_amps = [best_amps, bamp];
    best_taus = [best_taus, btau];
    best_gens = [best_gens, 1];
    best_chroms = [best_chroms; bchrom];
    
    new_gen = evolve(init_gen, amp_dels, tau_dels);
    
    while 1
        tic
        fprintf('\n Generation %i \n', cnt)
        cnt = cnt + 1;
        [fits, amp_dels, tau_dels] = eval_fn(new_gen, y);
        [bf, bf_idx] = min(fits);
        bamp = amp_dels(bf_idx);
        btau = tau_dels(bf_idx);
        bchrom = new_gen(bf_idx,:);
        
        % stopping tolerance
        if (bamp <= amp_tol) && (btau <= tau_tol)
            fprintf('Termination: %f|Amp: %f|Tau: %f \n', bf, bamp, btau);
            disp(bchrom)

            best_fits = [best_fits, bf];
            best_amps = [best_amps, bamp];
            best_taus = [best_taus, btau];
            best_gens = [best_gens, cnt];
            best_chroms = [best_chroms; bchrom];
            
            break
        end

        if (bf < best_fits(best_cnt))
            fprintf('Best fit is updated: %f|Amp: %f|Tau: %f \n', bf, bamp, btau)
            disp(bchrom)

            best_cnt = best_cnt + 1;
            best_fits = [best_fits, bf];
            best_amps = [best_amps, bamp];
            best_taus = [best_taus, btau];
            best_gens = [best_gens, cnt];
            best_chroms = [best_chroms; bchrom];
            % damp_factor = (bf-min_tol)/(max_tol-min_tol)*0.5 + 0.5;
            % fprintf('Damping factor %f \n', damp_factor)
        end 
        
        new_gen = evolve(new_gen, amp_dels, tau_dels);
        toc
    end 
end

function new_gen = evolve(chrom, amp_dels, tau_dels)
    global num_var;
    global N0;
    global N1;
    global N2;
    global sig_ctl;

    fits = amp_dels + tau_dels;
    [~, bf_idx] = min(fits);
    bamp = amp_dels(bf_idx);
    btau = tau_dels(bf_idx);
    
    [bfits, super_idx] = mink(fits, N1);
    elites = chrom(super_idx, :);
    sigs = std(elites);

    mean_elite = mean(elites, 1);
    elites(end,:) = mean_elite;

%     disp([bamp, btau])
%     disp(bfits)
%     disp(sigs)
%     disp(elites)

    sig_ctl = [sig_ctl; sigs];
    pooled_sigs = mean(sig_ctl,1);
%     disp(pooled_sigs)

    % breeding
    cnt = 1;
    new_gen = zeros(N0, num_var);
    new_gen(1:N1,:) = elites;
    for i=1:N1
        elite = elites(i,:);
        for j=1:N2
            offspring = elite + normrnd(0,pooled_sigs);
%             offspring(2) = abs(offspring(2));
%             offspring(4) = abs(offspring(4));
%             offspring(6) = abs(offspring(6));
            new_gen((N1+cnt),:) = offspring;
            cnt = cnt + 1;
        end
    end
end

function [fits, amp_dels, tau_dels] = eval_fn(chrom, y)
    global num_var;
    global N0;

    holding_p = -70; %mV
    holding_t = 450; %ms
    P1 = 50; %mV
    P1_t = 25*1000; % ms
    Ek = -91.1;
    
    fits = zeros(1, N0);
    amp_dels = zeros(1, N0);
    tau_dels = zeros(1, N0);
    wrn_idx = [];
    for i=1:N0
        try
            [t, ~, A] = IKsum(chrom(i,1:num_var), holding_p, holding_t, P1, P1_t, Ek);
            Ito = A(:,5);
            IKslow1 = A(:,10);
            IKslow2 = A(:,15);
            Iss = chrom(i,end);
            
            IKsum_trc = Ito + IKslow1 + IKslow2 + Iss;
            
            wrong_shape_iden = any(IKsum_trc < 0);
            [peak, peak_idx] = max(IKsum_trc);
            if (wrong_shape_iden == 1) || (t(peak_idx) < holding_t)
                % fprintf('Wrong shape at %i \n', i);
                % disp(chrom(i,:));
                    
                wrn_idx = [wrn_idx, i];
                amp_dels(i) = 15000;
                tau_dels(i) = 15000;
                fits(i) = amp_dels(i) + tau_dels(i);
            else
                amp_dels(i) = abs(peak - y(1));

                [~, tau_idx] = min(abs(peak*exp(-1) - IKsum_trc(peak_idx:end)));
                tau_dels(i) = abs(t(tau_idx)-y(2));
                
                fits(i) = amp_dels(i) + tau_dels(i);
            end
        catch
            % lastwarn
            % lasterr
            % fprintf('Error or warning at %i \n', i);
            % disp(chrom(i,:));

            wrn_idx = [wrn_idx, i];
            amp_dels(i) = 15000;
            tau_dels(i) = 15000;
            fits(i) = amp_dels(i) + tau_dels(i);
        end
    end
end
