function [best_fits, best_chroms] = custom_GA(nv, y, N0, N1, N2)
    global num_var;
    num_var = nv;

    low = [0.0, 0.0, 0.0, 20.0, 2.0];
    high = [70.0, 70.0, 70.0, 70.0, 50.0];
    init_gen = init_pop(low, high, N0);

    best_fits = [];
    best_chroms = [];
    
    cnt = 1;
    [fits, amp_dels, tau_dels] = eval_fn(init_gen, y, N0);
    [bf, bf_idx] = min(fits);
    
    fprintf('Initial best fit: %f|Amp: %f|Tai: %f \n', bf, amp_dels(bf_idx), tau_dels(bf_idx));
    best_fits = [best_fits, bf];
    best_chroms = [best_chroms; init_gen(bf_idx,:)];
    
    new_gen = evolve(init_gen, fits, N0, N1, N2);
    
    while 1
        fprintf('Generation %i \n', cnt)
        cnt = cnt + 1;
        [fits, amp_dels, tau_dels] = eval_fn(new_gen, y, N0);
        [bf, bf_idx] = min(fits);
        
        if (bf < best_fits(cnt-1))
            fprintf('Best fit is updated: %f|Amp: %f|Tai: %f \n', bf, amp_dels(bf_idx), tau_dels(bf_idx));
            best_fits = [best_fits, bf];
            best_chroms = [best_chroms; new_gen(bf_idx,:)];
            
            % stopping tolerance
            if bf < 0.1*2
                break
            end
        end 
        
        new_gen = evolve(new_gen, fits, N0, N1, N2);            
    end
end

function init_gen = init_pop(low, high, N0)
    global num_var;

    init_gen = zeros(N0, num_var);
    for j=1:num_var
        unif = makedist('Uniform', 'lower',low(j), 'upper',high(j));
        init_gen(:,j) = random(unif, N0, 1);
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
    for i=1:N1
        elite = elites(i,:);
        
        for j=1:N2
            offspring = elite + normrnd(0,1,[1,num_var]);
            new_gen((N1+cnt),:) = offspring;
            cnt = cnt + 1;
        end
    end    
end

function [fits, amp_dels, tau_dels] = eval_fn(chrom, y, N0)
    holding_p = -70; %mV
    holding_t = 450; %ms
    P1 = 50; %mV
    P1_t = 25*1000; % ms
    Ek = -80.3;
    
    fits = zeros(1, N0);
    amp_dels = zeros(1, N0);
    tau_dels = zeros(1, N0);
    wrn_idx = [];
    for i=1:N0
        try
            [t, ~, A, ~] = Ito(chrom(i,:), holding_p, holding_t, P1, P1_t, Ek);
            trc = A(:,5);
            
            wrong_shape_iden = any(trc < 0);
            if wrong_shape_iden == 1
                wrn_idx = [wrn_idx, i];
                amp_dels(i) = 15000;
                tau_dels(i) = 15000;
                fits(i) = amp_dels(i) + tau_dels(i);
            else
                [peak, ~] = max(trc);
                amp_dels(i) = abs(peak - y(1));

                [~, tau_idx] = min(abs(peak*exp(-1) - trc));
                tau_dels(i) = abs(t(tau_idx)-y(2));
                
                fits(i) = amp_dels(i) + tau_dels(i);
            end
        catch
            % lastwarn
            % lasterr
            fprintf('Error or warning at %i', i);
            wrn_idx = [wrn_idx, i];
            amp_dels(i) = 15000;
            tau_dels(i) = 15000;
            fits(i) = amp_dels(i) + tau_dels(i);
        end
    end
end
