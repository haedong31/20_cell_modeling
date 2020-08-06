function [best_amps, best_taus, best_gens, best_chroms] = IKslow_AGA(nv, y, n0, n1, n2)
    global num_var;
    global N0;
    global N1;
    global N2;
    global sig_ctl;
    num_var = nv;
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
    min_tol = amp_tol + tau_tol;

    % X0 = [22.5, 7.7, 45.2, 5.7, 2.058, 1200.0, 45.2, 5.7, 0.16]
    low = [-50.0, 2.0, -20.0, 2.0, 35000.0, 0.05];
    high = [70.0, 14.0, 80.0, 24.0, 500000.0, 0.3];
    init_gen = init_pop(low, high);
    
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
    
    max_tol = bf;
    new_gen = evolve(init_gen, amp_dels, tau_dels, 1);
    
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
        
        new_gen = evolve(new_gen, amp_dels, tau_dels, 1);
        toc
    end
end

function new_gen = evolve(chrom, amp_dels, tau_dels, damp_factor)
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
    sigs = damp_factor * std(elites);

    mean_elite = mean(elites, 1);
    elites(end,:) = mean_elite;

    disp([bamp, btau])
    disp(bfits)
    disp(sigs)
    disp(elites)

    sig_ctl = [sig_ctl; sigs];
    pooled_sigs = mean(sig_ctl,1);
    disp(pooled_sigs)

    % breeding
    cnt = 1;
    new_gen = zeros(N0, num_var);
    new_gen(1:N1,:) = elites;
    for i=1:N1
        elite = elites(i,:);
        for j=1:N2
            offspring = elite + normrnd(0,pooled_sigs);
            offspring(2) = abs(offspring(2));
            offspring(4) = abs(offspring(4));
            offspring(6) = abs(offspring(6));
            new_gen((N1+cnt),:) = offspring;
            cnt = cnt + 1;
        end
    end
end

function [fits, amp_dels, tau_dels] = eval_fn(chrom, y)
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
            [t, ~, A] = IKslow(chrom(i,:), holding_p, holding_t, P1, P1_t, Ek);
            trc = A(:,5);
            
            wrong_shape_iden = any(trc < 0);
            [peak, peak_idx] = max(trc);
            if (wrong_shape_iden == 1) || (t(peak_idx) < holding_t)
                % fprintf('Wrong shape at %i \n', i);
                % disp(chrom(i,:));
                    
                wrn_idx = [wrn_idx, i];
                amp_dels(i) = 15000;
                tau_dels(i) = 15000;
                fits(i) = amp_dels(i) + tau_dels(i);
            else
                amp_dels(i) = abs(peak - y(1));

                [~, tau_idx] = min(abs(peak*exp(-1) - trc(peak_idx:end)));
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

function init_gen = init_pop(low, high)
    global num_var;
    global N0;

    init_gen = zeros(N0, num_var);
    for j=1:num_var
        unif = makedist('Uniform', 'lower',low(j), 'upper',high(j));
        init_gen(:,j) = random(unif, N0, 1);
    end
end
