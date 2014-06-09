using DynamicFactorModels
using Distributions


# Montecarlo-tests: replicates tables 2 and 3 of Breitung and Eickmeier 2011
### test with no structural break - meassures size of the test
function montecarlo_test_size(B=1000, T=100, N=60, r=1, break_period=int(ceil(T/2)))
    montecarlo_lrs = zeros(B, N);
    montecarlo_lms = zeros(B, N);
    montecarlo_walds = zeros(B, N);
    for b in 1:B
        println("b: ", b, "(T, N): ", (T, N))
        y, x, f, lambda, epsilon_x = factor_model_DGP(T, N, r; model="Breitung_Eickmeier_2011", b=0)
        w = reshape(ones(length(y)), (length(y), 1))
        dfm = DynamicFactorModel(y,w,x,1)  # here we could estimate the number of factors which gives strongly different results especially for small N and T
        montecarlo_lrs[b, :] = [LR_test(dfm, break_period, i) for i in 1:N]
        montecarlo_lms[b, :] = [LM_test(dfm, break_period, i) for i in 1:N]
        montecarlo_walds[b, :] = [Wald_test(dfm, break_period, i) for i in 1:N]
    end
    critical_value = quantile(Distributions.Chisq(r), 0.95)
    montecarlo_lr_size = mean(sum((montecarlo_lrs.>critical_value), 2)/N)
    montecarlo_lm_size = mean(sum((montecarlo_lms.>critical_value), 2)/N)
    montecarlo_wald_size = mean(sum((montecarlo_walds.>critical_value), 2)/N)
    return [montecarlo_lr_size, montecarlo_lm_size, montecarlo_wald_size]
end

### test with structural break - meassures power of the test
function montecarlo_test_power(B=1000, T=100, N=60, r=1, b=1, break_period=int(ceil(T/2)))
    montecarlo_lrs = zeros(B, N);
    montecarlo_lms = zeros(B, N);
    montecarlo_walds = zeros(B, N);
    for rep in 1:B
        y, x, f, lambda, epsilon_x = factor_model_DGP(T, N, r; model="Breitung_Eickmeier_2011", b=b)
        w = reshape(ones(length(y)), (length(y), 1))
        dfm = DynamicFactorModel(y,w,x,1)  # here we could estimate the number of factors which gives strongly different results especially for small N and T
        montecarlo_lrs[rep, :] = [LR_test(dfm, break_period, i) for i in 1:N]
        montecarlo_lms[rep, :] = [LM_test(dfm, break_period, i) for i in 1:N]
        montecarlo_walds[rep, :] = [Wald_test(dfm, break_period, i) for i in 1:N]
    end
    critical_value = quantile(Distributions.Chisq(r), 0.95)
    montecarlo_lr_size = mean(sum((montecarlo_lrs.>critical_value), 2)/N)
    montecarlo_lm_size = mean(sum((montecarlo_lms.>critical_value), 2)/N)
    montecarlo_wald_size = mean(sum((montecarlo_walds.>critical_value), 2)/N)
    return [montecarlo_lr_size, montecarlo_lm_size, montecarlo_wald_size]
end

### test with no structural break - meassures size of the test
function bootstrap_test_size(R=1000, B=1000, T=100, N=60, r=1, break_period=int(ceil(T/2)))
    bootstrap_lr_p_values = zeros(R, N)
    bootstrap_lm_p_values = zeros(R, N)
    bootstrap_wald_p_values = zeros(R, N)
    for rep in 1:R
        println("rep: ", rep, "(T, N): ", (T, N))
        y, x, f, lambda, epsilon_x = factor_model_DGP(T, N, r; model="Breitung_Eickmeier_2011", b=0)
        w = reshape(ones(length(y)), (length(y), 1))
        dfm = DynamicFactorModel(y,w,x,1)  # here we could estimate the number of factors which gives strongly different results especially for small N and T
        lr_stats = [LR_test(dfm, break_period, i) for i in 1:N]
        lm_stats = [LM_test(dfm, break_period, i) for i in 1:N]
        wald_stats = [Wald_test(dfm, break_period, i) for i in 1:N]
        bootstrap_lrs = apply(hcat, [residual_bootstrap(dfm, B, dfm->LR_test(dfm, break_period, i)) for i in 1:N])'
        bootstrap_lms = apply(hcat, [residual_bootstrap(dfm, B, dfm->LM_test(dfm, break_period, i)) for i in 1:N])'
        bootstrap_walds = apply(hcat, [residual_bootstrap(dfm, B, dfm->Wald_test(dfm, break_period, i)) for i in 1:N])'
        bootstrap_lr_p_values[rep, :] = [mean(bootstrap_lrs[i, :].>lr_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
        bootstrap_lm_p_values[rep, :] = [mean(bootstrap_lms[i, :].>lm_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
        bootstrap_wald_p_values[rep, :] = [mean(bootstrap_walds[i, :].>wald_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
    end
    avg_rejections_lr = mean(bootstrap_lr_p_values.<0.05)  # bootstrap_lr_p_values holds bootstrapped p_values for each variable per repetition, we average over repetitions (and variables to get something comparable to table 3)
    avg_rejections_lm = mean(bootstrap_lm_p_values.<0.05)
    avg_rejections_wald = mean(bootstrap_wald_p_values.<0.05)
    critical_value = quantile(Distributions.Chisq(r), 0.95)
    return [avg_rejections_lr, avg_rejections_lm, avg_rejections_wald]
end

### test with structural break - meassures power of the test
function bootstrap_test_power(R=1000, B=1000, T=100, N=60, r=1, b=1, break_period=int(ceil(T/2)))
    bootstrap_lr_p_values = zeros(R, N)
    bootstrap_lm_p_values = zeros(R, N)
    bootstrap_wald_p_values = zeros(R, N)
    for rep in 1:R
        println("rep: ", rep, "(T, N): ", (T, N))
        y, x, f, lambda, epsilon_x = factor_model_DGP(T, N, r; model="Breitung_Eickmeier_2011", b=b)
        w = reshape(ones(length(y)), (length(y), 1))
        dfm = DynamicFactorModel(y,w,x,1)  # here we could estimate the number of factors which gives strongly different results especially for small N and T
        lr_stats = [LR_test(dfm, break_period, i) for i in 1:N]
        lm_stats = [LM_test(dfm, break_period, i) for i in 1:N]
        wald_stats = [Wald_test(dfm, break_period, i) for i in 1:N]
        bootstrap_lrs = apply(hcat, [residual_bootstrap(dfm, B, dfm->LR_test(dfm, break_period, i)) for i in 1:N])'
        bootstrap_lms = apply(hcat, [residual_bootstrap(dfm, B, dfm->LM_test(dfm, break_period, i)) for i in 1:N])'
        bootstrap_walds = apply(hcat, [residual_bootstrap(dfm, B, dfm->Wald_test(dfm, break_period, i)) for i in 1:N])'
        bootstrap_lr_p_values[rep, :] = [mean(bootstrap_lrs[i, :].>lr_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
        bootstrap_lm_p_values[rep, :] = [mean(bootstrap_lms[i, :].>lm_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
        bootstrap_wald_p_values[rep, :] = [mean(bootstrap_walds[i, :].>wald_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
    end
    avg_rejections_lr = mean(bootstrap_lr_p_values.<0.05)  # bootstrap_lr_p_values holds bootstrapped p_values for each variable per repetition, we average over repetitions (and variables to get something comparable to table 3)
    avg_rejections_lm = mean(bootstrap_lm_p_values.<0.05)
    avg_rejections_wald = mean(bootstrap_wald_p_values.<0.05)
    critical_value = quantile(Distributions.Chisq(r), 0.95)
    return [avg_rejections_lr, avg_rejections_lm, avg_rejections_wald]
end



# replicate table 2 of Breitung, Eickmeier 2011
function table2_montecarlo()
    B = 1000
    r = 1
    Ts = [50, 100, 150, 200]
    Ns = [20, 50, 100, 150, 200]
    table = zeros(5*4, 3)
    for T_ind in 1:length(Ts), N_ind in 1:length(Ns)
        T = Ts[T_ind]; N = Ns[N_ind]
        table[(T_ind-1) * length(Ns) + N_ind, :] = montecarlo_test_size(B, T, N, r)
    end
    return table
end



# replicate table 3 of Breitung, Eickmeier 2011
function table3_montecarlo()
    B = 1000
    r = 1
    Ts = [50, 100, 150, 200]
    bs = [0.1, 0.2, 0.3, 0.5]
    N = 50
    table = zeros(4*4, 3)
    for T_ind in 1:length(Ts), b_ind in 1:length(bs)
        T = Ts[T_ind]
        table[(T_ind-1)*length(bs) + b_ind, :] = montecarlo_test_power(B, T, N, r, bs[b_ind])
    end
    return table
end



# Bootstrap: rerun monte-carlo of tables 2 and 3 but this time we bootstrap the tests hoping for performance increases with lower sample sizes
# replicate table 2 of Breitung, Eickmeier 2011 - Bootstrapped
function table2_bootstrap(B=100, R=1, r=1, Ts=[50, 100, 150, 200], Ns=[20, 50, 100, 150, 200])
    # replicate table 2 of Breitung, Eickmeier 2011 - Bootstrapped
    table = zeros(length(Ns)*length(Ts), 3)
    for T_ind in 1:length(Ts), N_ind in 1:length(Ns)
        T = Ts[T_ind]; N = Ns[N_ind]
        table[(T_ind-1) * length(Ns) + N_ind, :] = bootstrap_test_size(R, B, T, N, r)
    end
    return table
end

# replicate table 3 of Breitung, Eickmeier 2011 - Bootstrapped
function table3_bootstrap(B=100, R=1, r=1, Ts=[50, 100, 150, 200], bs=[0.1, 0.2, 0.3, 0.5])
    N = 50
    table = zeros(length(Ts)*length(bs), 3)
    for T_ind in 1:length(Ts), b_ind in 1:length(bs)
        T = Ts[T_ind]
        table[(T_ind-1) * length(bs) + b_ind, :] = bootstrap_test_power(R, B, T, N, r, bs[b_ind])
    end
    return table
end


model_estimation_speed(y,w,x,r=1) = begin  # returns number of estimated models per second
    tic();
    for i in 1:1000
        dfm = DynamicFactorModel(y,w,x,1)
    end
    toc()/1000
end

complexity_bootstrap_table2(Ns, Ts, B=1000, R=1000) = sum([R*(3*B*Ns[i] + 3*Ns[i]) for i in 1:length(Ns), j in 1:length(Ts)])  # number of models estimated
complexity_bootstrap_table3(Ts, bs, B=1000, R=1000) = sum([R*(3*B*50 + 3*50) for i in 1:length(bs), j in 1:length(Ts)])  # number of models estimated
