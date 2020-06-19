function [best_fits, best_chrom] = IKsum_AGA(trcy, nv, iters, N0, N1, N2)
    global num_var;
    num_var = nv;

    init_to = [-13.5655  128.4098  321.7877  127.2189   58.4796];
    init_Kslow1 = [-0.0613    0.0097    0.2070    0.0128    1.1628];
    init_Kslow1 = init_Kslow1*1000;
    init_Kslow2 = [-0.0717    0.0123    0.0245    0.0399    8.6985];
    init_Kslow2 = init_Kslow2*1000;
    init_param = [init_to init_Kslow1 init_Kslow2];

    % initial population
    init_gen = zeros(N0, num_var);
    init_gen(1,:) = init_param;
    for i=1:num_var
        init_gen(2:end,i) = init_param(i) + normrnd(0,1,[N0-1,1]);
    end

    best_fits = [];
    cnt = 1;
    fits = eval_fn(init_gen, trcy, N0);
    
    [bf, bf_idx] = min(fits);
    best_fits = [best_fits, bf];
    best_chrom = init_gen(bf_idx,:);
    fprintf('Initial best fit: %f \n', bf);
    disp(best_chrom)
    best_cnt = 1;

    new_gen = evolve(init_gen, fits, N0, N1, N2);
    while 1
        tic
        fprintf('\n Generation %i \n', cnt)
        
        fits = eval_fn(new_gen, trcy, N0);
        [bf, bf_idx] = min(fits);
        
        if (bf < best_fits(best_cnt))
            fprintf('Best fit is updated: %f \n', bf);
            best_fits = [best_fits, bf];
            best_chrom = new_gen(bf_idx,:);
            disp(best_chrom)
            best_cnt = best_cnt + 1;
        end

        if cnt == iters
            break
        end
        new_gen = evolve(new_gen, fits, N0, N1, N2);
        cnt = cnt + 1;
        toc
    end
end

function new_gen = evolve(chrom, fits, N0, N1, N2)
    global num_var;

    new_gen = zeros(N0, num_var);
    [~, super_idx] = mink(fits, N1);
    elites = chrom(super_idx, :);
    new_gen(1:N1,:) = elites;

    % breeding
    cnt = 1;
    sigs = std(elites);
    for i=1:N1
        elite = elites(i,:);
        for j=1:N2
            offspring = elite + normrnd(0,sigs);
            offspring(7) = abs(offspring(7));
            offspring(9) = abs(offspring(9));
            offspring(12) = abs(offspring(12));
            offspring(14) = abs(offspring(14));
            new_gen((N1+cnt),:) = offspring;
            cnt = cnt + 1;
        end
    end    
end

function fits = eval_fn(chrom, trcy, N0)
    holding_p = -70; %mV
    holding_t = 450; %ms
    P1 = 50; %mV
    P1_t = 25*1000; % ms
    Ek = -91.1;
    
    fits = zeros(1, N0);
    for i=1:N0
        try
            [t, ~, A, ~] = IKsum(chrom(i,:), holding_p, holding_t, P1, P1_t, Ek);
            trc = A(:,5) + A(:,10) + A(:,15);
            
            wrong_shape_iden = any(trc < 0);
            [peak, peak_idx] = max(trc);
            if (wrong_shape_iden == 1) || (t(peak_idx) < holding_t)
                fprintf('Wrong shape at %i \n', i);
                disp(chrom(i,:));
                fits(i) = 1.0e+5;
            else
                fits(i) = log(dtw(trcy, trc));
            end
        catch
            % lastwarn
            % lasterr
            fprintf('Error or warning at %i \n', i);
            disp(chrom(i,:));
            fits(i) = 1.0e+5;
        end
    end
end
