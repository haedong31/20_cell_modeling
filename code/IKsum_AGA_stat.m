function [best_amps, best_taus, best_gens, best_chroms] = IKsum_AGA(nv, y, n0, n1, n2)
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
    
    % initial population
    low = [0.0, 0.0, 0.0, 20.0, 2.0, 0.2,...
           -50.0, 2.0, -20.0, 2.0, 200.0, 0.05,...
           -50.0, 2.0, -20.0, 2.0, 1000.0, 0.05];
    high = [70.0, 70.0, 70.0, 70.0, 50.0, 0.6,...
            70.0, 14.0, 80.0, 24.0, 5000.0, 0.3,...
            70.0, 14.0, 80.0, 24.0, 10000.0, 0.3];
    init_gen = init_pop(low, high);

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
    for i=1:N0
        try
            [t, ~, A] = IKsum(chrom(i,:), holding_p, holding_t, P1, P1_t, Ek);
            trc = A(:,5) + A(:,10) + A(:,15);

            wrong_shape_iden = any(trc < 0);
            [peak, peak_idx] = max(trc);
            if (wrong_shape_iden == 1) || (t(peak_idx) < holding_t)
                % fprintf('Wrong shape at %i \n', i);
                % disp(chrom(i,:));

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
            fprintf('Error or warning at %i \n', i);
            % disp(chrom(i,:));

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
