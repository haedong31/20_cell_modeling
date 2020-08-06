function [best_fits, best_gens, best_chroms] = tri_exp(nv, trc, n0, n1, n2)
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
    best_gens = [];
    best_chroms = [];
    
    % initialization
    % Iss, A1, A2, A3, Tau1, Tau2, Tau3
    % A1 - IKslow2; A2 - IKslow1; A3 - Ito
    low = [2.08, 1.81, 0.43, 2.13, 4476.04, 262.21,	84.47];
    high = [13.50, 6.35, 12.52,	60.99, 57941.35, 2000.97, 203.91];

    init_gen = zeros(N0, num_var);
    for j=1:num_var
        unif = makedist('Uniform', 'lower',low(j), 'upper',high(j));
        init_gen(:,j) = random(unif, N0, 1);
    end
    
    cnt = 1;
    fits = eval(init_gen, trc);
    [bf, bf_idx] = min(fits);
    bchrom = init_gen(bf_idx,:);
    
    fprintf('Initial best fit: %f \n', bf)
    disp(bchrom)
    
    best_cnt = 1;
    best_fits = [best_fits, bf];
    best_gens = [best_gens, 1];
    best_chroms = [best_chroms, bchrom];
    
    new_gen = breed(init_gen, fits);
    
    while 1
        cnt = cnt + 1;
        fits = eval(new_gen, trc);
        [bf, bf_idx] = min(fits);
        bchrom = init_gen(bf_idx,:);
        
        % stopping
        if (bf <= 1) || (cnt == 10000)
            fprintf('Termination: %f \n', bf);
            disp(bchrom)
            
            best_cnt = best_cnt + 1;
            best_fits = [best_fits, bf];
            best_gens = [best_gens, cnt];
            best_chroms = [best_chroms; bchrom];
            
            break
        end
        
        % update the best
        if (bf < best_fits(best_cnt))
%             fprintf('\n Generation %i \n', cnt)
%             fprintf('Best fit is updated: %f \n', bf)
%             disp(bchrom)
            
            best_cnt = best_cnt + 1;
            best_fits = [best_fits, bf];
            best_gens = [best_gens, cnt];
            best_chroms = [best_chroms; bchrom];
        end
        
        new_gen = breed(new_gen, fits);
    end
end


function new_gen = breed(chrom, fits)
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
    pooled_sigs = mean(sig_ctl, 1);
    
%     disp(super_fits)
%     disp(sigs)
%     disp(pooled_sigs)
    
    cnt = 1;
    new_gen = zeros(N0, num_var);
    new_gen(1:N1,:) = elites;
    for i=1:N1
        elite = elites(i,:);
        for j=1:N2
            offspring = elite + normrnd(0,pooled_sigs);
            offspring = abs(offspring);
            new_gen((N1+cnt),:) = offspring;
            cnt = cnt + 1;
        end
    end        
end


function fits = eval(chrom, trc)
    global N0;
    
    t = trc.time;
    fits = zeros(1, N0);
    for i=1:N0
        param = chrom(i,:);
        Iss = param(1);
        IKslow2 = exp_fn(t, param(2), param(5));
        IKslow1 = exp_fn(t, param(3), param(6));
        Ito = exp_fn(t, param(4), param(7));
        
        IKsum = Iss + IKslow2 + IKslow1 + Ito;
        fits(i) = sum((IKsum - trc.current).^2);
    end
end
