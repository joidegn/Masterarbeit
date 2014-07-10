using FactorModels
using Distributions

# Montecarlo-tests: used for replicating tables 2 and 3 of Breitung and Eickmeier 2011
### test with no structural break - meassures size of the test

function montecarlo_test_statistics(R=1000, T=100, N=60, r=1, break_period=int(ceil(T/2)), b=0)  # returns 3 arrays of RxN test statistics
    montecarlo_lrs = zeros(R, N)
    montecarlo_lms = zeros(R, N)
    montecarlo_walds = zeros(R, N)
    for rep in 1:R
        println("rep: ", rep, "(T, N): ", (T, N))
        y, x, f, lambda, epsilon_x = factor_model_DGP(T, N, r; model="Breitung_Eickmeier_2011", b=b)
        fm = FactorModel(x, 1)  # here we could estimate the number of factors which gives strongly different results especially for small N and T
        montecarlo_lrs[rep, :] = [LR_test(fm, break_period, i) for i in 1:N]
        montecarlo_lms[rep, :] = [LM_test(fm, break_period, i) for i in 1:N]
        montecarlo_walds[rep, :] = [Wald_test(fm, break_period, i) for i in 1:N]
    end
    return (montecarlo_lrs, montecarlo_lms, montecarlo_walds)
end


### test with structural break - meassures power of the test
function montecarlo_test_power(R=1000, T=100, N=60, r=1, b=1, break_period=int(ceil(T/2)))
    montecarlo_lrs, montecarlo_lms, montecarlo_walds = montecarlo_test_statistics(R, T, N, r, break_period, b)
    critical_value = quantile(Distributions.Chisq(r), 0.95)
    montecarlo_lr_power = mean(montecarlo_lrs.>critical_value)
    montecarlo_lm_power = mean(montecarlo_lms.>critical_value)
    montecarlo_wald_power = mean(montecarlo_walds.>critical_value)
    return [montecarlo_lr_power, montecarlo_lm_power, montecarlo_wald_power]
end

### test with no structural break - meassures size of the test
montecarlo_test_size(R=1000, T=100, N=60, r=1, break_period=int(ceil(T/2))) = montecarlo_test_power(R, T, N, r, 0, break_period)


function bootstrap_test_statistics(fm, B, break_period=int(ceil(T/2)))
    bootstrap_lrs = apply(hcat, [residual_bootstrap(fm, B, fm->LR_test(fm, break_period, i)) for i in 1:size(fm.x, 2)])'
    bootstrap_lms = apply(hcat, [residual_bootstrap(fm, B, fm->LM_test(fm, break_period, i)) for i in 1:size(fm.x, 2)])'
    bootstrap_walds = apply(hcat, [residual_bootstrap(fm, B, fm->Wald_test(fm, break_period, i)) for i in 1:size(fm.x, 2)])'
    return (bootstrap_lrs, bootstrap_lms, bootstrap_walds, )
end

function bootstrap_test_rejections(R=1000, B=1000, T=100, N=60, r=1, b=0, break_period=int(ceil(T/2)))
    println("bootstrapping test rejections using b=", b, ", B=", B, " and R=", R)
    bootstrap_lrs = Array(Float64, (R, N, B))
    bootstrap_lms = Array(Float64, (R, N, B))
    bootstrap_walds = Array(Float64, (R, N, B))
    bootstrap_lr_p_values = Array(Float64, (R, N))
    bootstrap_lm_p_values = Array(Float64, (R, N))
    bootstrap_wald_p_values = Array(Float64, (R, N))
    for rep in 1:R
        println("rep: ", rep, "(T, N): ", (T, N))
        y, x, f, lambda, epsilon_x = factor_model_DGP(T, N, r; model="Breitung_Eickmeier_2011", b=b)
        fm = FactorModel(x, r)  # here we could estimate the number of factors which gives strongly different results especially for small N and T
        lr_stats = [LR_test(fm, break_period, i) for i in 1:N]
        lm_stats = [LM_test(fm, break_period, i) for i in 1:N]
        wald_stats = [Wald_test(fm, break_period, i) for i in 1:N]
        bootstrap_lrs[rep, :, :], bootstrap_lms[rep, :, :], bootstrap_walds[rep, :, :] = bootstrap_test_statistics(fm, B, break_period)
        bootstrap_lr_p_values[rep, :] = [mean(bootstrap_lrs[rep, i, :].>lr_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
        bootstrap_lm_p_values[rep, :] = [mean(bootstrap_lms[rep, i, :].>lm_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
        bootstrap_wald_p_values[rep, :] = [mean(bootstrap_walds[rep, i, :].>wald_stats[i]) for i in 1:N]  # estimated p_value as in Davidson - MacKinnon p. 158 (1 - F_hat(tau_hat) = 1/B * sum(tau_start .> tau_hat))
    end
    bootstrap_lr_rejections = bootstrap_lr_p_values .< 0.05
    bootstrap_lm_rejections = bootstrap_lm_p_values .< 0.05
    bootstrap_wald_rejections = bootstrap_wald_p_values .< 0.05
    return (bootstrap_lr_rejections, bootstrap_lm_rejections, bootstrap_wald_rejections)
end

### test with structural break - meassures power of the test
function bootstrap_test_power(R=1000, B=1000, T=100, N=60, r=1, b=1, break_period=int(ceil(T/2)))
    bootstrap_lr_rejections, bootstrap_lm_rejections, bootstrap_wald_rejections = bootstrap_test_rejections(R, B, T, N, r, b, break_period)
    return [mean(bootstrap_lr_rejections), mean(bootstrap_lm_rejections), mean(bootstrap_wald_rejections)]
end

### test with no structural break - meassures size of the test
bootstrap_test_size(R=1000, B=1000, T=100, N=60, r=1, break_period=int(ceil(T/2))) = bootstrap_test_power(R, B, T, N, r, 0, break_period)

# replicate table 2 of Breitung, Eickmeier 2011
function table2_montecarlo(R=1000, r=1, Ts=[50, 100, 150, 200], Ns=[20, 50, 100, 150, 200])
    table = zeros(5*4, 3)
    for T_ind in 1:length(Ts), N_ind in 1:length(Ns)
        println("now doing: T=", Ts[T_ind], "N:", Ns[N_ind])
        T = Ts[T_ind]; N = Ns[N_ind]
        table[(T_ind-1) * length(Ns) + N_ind, :] = montecarlo_test_size(R, T, N, r)
    end
    return table
end
# get all test statistic values for table 2 (to check their distribution)
table2_stats_montecarlo(R=1000, r=1, Ts=[50, 100, 150, 200], Ns=[20, 50, 100, 150, 200]) = [montecarlo_test_statistics(R, T, N, r) for T in Ts, N in Ns]

# replicate table 3 of Breitung, Eickmeier 2011
function table3_montecarlo(R=1000, r=1, Ts=[50, 100, 150, 200], bs=[0.1, 0.2, 0.3, 0.5], N=50)
    table = zeros(4*4, 3)
    for T_ind in 1:length(Ts), b_ind in 1:length(bs)
        println("now doing: T=", Ts[T_ind], "b:", bs[b_ind])
        T = Ts[T_ind]
        table[(T_ind-1)*length(bs) + b_ind, :] = montecarlo_test_power(R, T, N, r, bs[b_ind])
    end
    return table
end
# get all test statistic values for table 3 (to check their distribution)
table3_stats_montecarlo(R=1000, r=1, Ts=[50, 100, 150, 200], bs=[0.1, 0.2, 0.3, 0.5], N=50) = [montecarlo_test_statistics(R, T, N, r) for T in Ts, b in bs]

# Bootstrap: rerun monte-carlo of tables 2 and 3 but this time we bootstrap the tests hoping for performance increases with lower sample sizes
# replicate table 2 of Breitung, Eickmeier 2011 - Bootstrapped
function table2_bootstrap(B=100, R=1000, r=1, Ts=[50, 100, 150, 200], Ns=[20, 50, 100, 150, 200])
    table = zeros(length(Ns)*length(Ts), 3)
    for T_ind in 1:length(Ts), N_ind in 1:length(Ns)
        println("now doing: T=", Ts[T_ind], "N:", Ns[N_ind])
        T = Ts[T_ind]; N = Ns[N_ind]
        table[(T_ind-1) * length(Ns) + N_ind, :] = bootstrap_test_size(R, B, T, N, r)
    end
    return table
end
table2_stats_bootstrap(R=1000, B=1000, r=1, Ts=[50, 100, 150, 200], Ns=[20, 50, 100, 150, 200]) = [bootstrap_test_statistics(R, B, T, N, r) for T in Ts, N in Ns]

# replicate table 3 of Breitung, Eickmeier 2011 - Bootstrapped
function table3_bootstrap(B=100, R=1000, r=1, Ts=[50, 100, 150, 200], bs=[0.1, 0.2, 0.3, 0.5], N=50)
    table = zeros(length(Ts)*length(bs), 3)
    for T_ind in 1:length(Ts), b_ind in 1:length(bs)
        println("now doing: T=", Ts[T_ind], "b:", bs[b_ind])
        T = Ts[T_ind]
        table[(T_ind-1) * length(bs) + b_ind, :] = bootstrap_test_power(R, B, T, N, r, bs[b_ind])
    end
    return table
end
table3_stats_bootstrap(R=1000, B=1000, r=1, Ts=[50, 100, 150, 200], bs=[0.1, 0.2, 0.3, 0.5], N=50) = [bootstrap_test_statistics(R, B, T, N, r, b) for T in Ts, b in bs]


#model_estimation_speed(y,w,x,r=1) = begin  # returns number of estimated models per second
#    tic();
#    for i in 1:1000
#        fm = FactorModel(x, 1)
#    end
#    toc()/1000
#end
#
## Note: meassures and actual time dont correspond perfectly
#complexity_bootstrap_table2(Ns, Ts, B=1000, R=1000) = sum([R*(3*B*Ns[i] + 3*Ns[i]) for i in 1:length(Ns), j in 1:length(Ts)])  # number of models estimated
#complexity_bootstrap_table3(Ts, bs, B=1000, R=1000) = sum([R*(3*B*50 + 3*50) for i in 1:length(bs), j in 1:length(Ts)])  # number of models estimated

Ts = [50, 100, 150, 200]
Ns = [20, 50, 100, 150, 200]
bs = [0.1, 0.2, 0.3, 0.5]

table2_monte = table2_montecarlo()
table2_stats_monte = table2_stats_montecarlo(5)
#table2_lrs = [table2_stats_monte[T_ind, N_ind][1] for T_ind in 1:length(Ts), N_ind in 1:length(Ns)]
table3_monte = table3_montecarlo()
#table3_stats_monte = table3_stats_montecarlo(1)
table2_boot = table2_bootstrap()
#table2_stats_boot = table2_stats_bootstrap()
table3_boot = table3_bootstrap(10, 1)
#table3_stats_boot = table3_stats_bootstrap()
