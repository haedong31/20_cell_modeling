function [best_fits, best_iters, best_chroms] = direct_decomp(nv_to, nv_Kslow, trc, n0, n1, n2)
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
    best_iters = [];
    best_chroms = [];

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
        init_gen_Kslow1(:,j) = random(unif2, N0, 1);
        init_gen_Kslow2(:,j) = random(unif2, N0, 1);
    end
    
    % Iss initialization
    unif3 = makedist('Uniform', 'lower',2.08, 'upper',13.50);
    init_ss = random(unif3, N0, 1);
    
    % consolidated initial population
    init_gen = [init_gen_to, init_gen_Kslow1, init_gen_Kslow2, init_ss];
    
    cnt = 1;
    fits = eval_fn(init_gen, trc);
    [bf, bf_idx] = min(fits);
    bchrom = init_gen(bf_idx,:);

    fprintf('Initial best fit: %f \n', bf);
    disp(bchrom)

    best_cnt = 1;
    best_fits = [best_fits, bf];
    best_iters = [best_iters, 1];
    best_chroms = [best_chroms; bchrom];
    
    new_gen = evolve(init_gen, fits);
    
    while 1
        tic
        cnt = cnt + 1;
        fprintf('\n Generation %i \n', cnt)
        
        fits = eval_fn(new_gen, trc);
        [bf, bf_idx] = min(fits);
        bchrom = new_gen(bf_idx,:);
        
        % stopping criterion
        if (cnt <= 1000) || (bf <= 1)
            fprintf('Termination: %f \n', bf);
            disp(bchrom)

            best_fits = [best_fits, bf];
            best_iters = [best_iters, cnt];
            best_chroms = [best_chroms; bchrom];
            
            break
        end

        if (bf < best_fits(best_cnt))
            fprintf('Best fit is updated: %f \n', bf)
            disp(bchrom)

            best_cnt = best_cnt + 1;
            best_fits = [best_fits, bf];
            best_iters = [best_iters, cnt];
            best_chroms = [best_chroms; bchrom];
        end 
        
        new_gen = evolve(new_gen, fits);
        toc
    end 
end


function new_gen = evolve(chrom, fits)
    global num_var;
    global N0;
    global N1;
    global N2;
    global sig_ctl;

    [~, super_idx] = mink(fits, N1);
    elites = chrom(super_idx, :);
    
    mean_elite = mean(elites, 1);
    elites(end,:) = mean_elite;
    
    sigs = std(elites);
    sig_ctl = [sig_ctl; sigs];
    pooled_sigs = mean(sig_ctl,1);

%     disp(elites)
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


function fits = eval_fn(chrom, trc)
    global num_var;
    global N0;

    holding_p = -70; %mV
    holding_t = 450; %ms
    P1 = 50; %mV
    P1_t = 25*1000; % ms
    Ek = -91.1;
    
    fits = zeros(1, N0);
    for i=1:N0
        try
            [t, ~, A] = IKsum(chrom(i,1:num_var), holding_p, holding_t, P1, P1_t, Ek);
            Ito = A(:,5);
            IKslow1 = A(:,10);
            IKslow2 = A(:,15);
            Iss = chrom(i,end);
            
            IKsum_trc = Ito + IKslow1 + IKslow2;
            
            wrong_shape_iden = any(IKsum_trc < 0);
            [~, peak_idx] = max(IKsum_trc);
            if (wrong_shape_iden == 1) || (t(peak_idx) < holding_t)
                % fprintf('Wrong shape at %i \n', i);
                % disp(chrom(i,:));
                fits(i) = 100000000;
            else
                % chop off the early phase
                IKsum_trc_rd = (IKsum_trc(peak_idx:end) + Iss);
                t_rd = t(peak_idx:end);
                t_rd = t_rd - t_rd(1);

                appx_time_idx = dsearchn(t_rd, trc.time);
                ds_IKsum = IKsum_trc_rd(appx_time_idx);

                fits(i) = sum((ds_IKsum - trc.current).^2);
            end
        catch
            % lastwarn
            % lasterr
            % fprintf('Error or warning at %i \n', i);
            % disp(chrom(i,:));
            fits(i) = 100000000;
        end
    end
end
