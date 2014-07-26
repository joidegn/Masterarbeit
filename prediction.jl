using FactorModels
using DataFrames
using Datetime
using DataStructures  # for ordered dict which makes reading easier


# Models:
# - AR(p)
# - Static FM with no breaks, number of factors according to ICp1
# - Static FM with 1 breaks, number of factors according to ICp1 on all subperiods
# ??- Static FM with more breaks?, number of factors according to ICp1 on all subperiods
# - Dynamic FM with no breaks, number of factors according to ICp1 AND number of factor lags according to Breitung, Eickmeier p. 75
# - Dynamic FM with breaks, number of factors according to ICp1 AND number of factor lags according to Breitung, Eickmeier p. 75 --> loop through and let both number or factors and factor lags converge?? (or get number of factors first and then factor lags)

#data = readtable("data/final_data.csv", separator=' ', header=true)  # this is how it used to be
#data_matrix = reduce(hcat, [convert(Array{Float64}, col) for col in data.columns[2:size(data.columns, 1)]])

file = readcsv("data/final_data.csv")
data, header, dates = convert(Array{Float64, 2}, file[2:end, 2:end]), file[1, 2:end], Date[date(dt) for dt in file[2:end, 1]]
#ids = map(string, names(data)[2:end])
#titles = series_titles(ids)  # TODO: does not work at the moment because ICU.jl and with it Requests.jl seems to be broken
y_ind = 1
x = convert(Array{Float64, 2}, data)  # x includes y
y = data[:, 1]

number_of_forecasts = 15

BIC(mse::Float64, n::Int64, T) = log(mse) + n * log(T)/T  # Bayesian information criterion
FPE(mse::Float64, n::Int64, T) = log(mse) + n * 2/T  # Fisher information criterion

include("choose_model.jl")


# Benchmark MSE
bench_forecasts, true_values = benchmark_forecasts(x, 1; num_predictions=number_of_forecasts, window=number_of_forecasts)
rmse_bench = RMSE(bench_forecasts, true_values)

different_ar_forecasts = [benchmark_ar(x[:, 1], i; num_predictions=number_of_forecasts)[1] for i in 1:10]
bench_ar_forecasts, true_values_ = benchmark_ar(x[:, 1], 4; num_predictions=number_of_forecasts)  # lazy way to get the true values
rmses = [RMSE(ar_forecasts, true_values) for ar_forecasts in different_ar_forecasts]
num_lags = indmin(rmses)
bench_ar_forecasts, true_values_ = benchmark_ar(x[:, 1], 6; num_predictions=number_of_forecasts)
rmse_bench_ar = RMSE(bench_ar_forecasts, true_values)



# Static model in-sample model selection
#static_model_in_sample_bic = choose_static_factor_model_in_sample(x, BIC)
#static_model_in_sample_bic_icp2 = choose_static_factor_model_in_sample(x, BIC, "ICp2")
#static_model_in_sample_fpe = choose_static_factor_model_in_sample(x, FPE)
#static_model_in_sample_fpe_icp2 = choose_static_factor_model_in_sample(x, FPE, "ICp2")
#static_model_in_sample_icp2 = choose_static_factor_model_in_sample(x, FPE, "ICp2"; bai_criterion=true, max_lags=3)  # if bai_criterion is true information criterion (FPE) plays no role
# Static model out-of-sample model selection
#static_model_out_of_sample = choose_static_factor_model_out_of_sample(x)
#static_model_out_of_sample_icp2 = choose_static_factor_model_out_of_sample(x, "ICp2")


# Dynamic model in-sample model selection
#dynamic_model_in_sample_bic = choose_dynamic_factor_model_in_sample(x, BIC; max_factors=10, max_factor_lags=4)  # if max factors is higher residuals will be tiny (non distibguishable from 0, clear overfit)
#dynamic_model_in_sample_bic_icp2 = choose_dynamic_factor_model_in_sample(x, BIC, "ICp2"; max_factors=10, max_factor_lags=4)
#dynamic_model_in_sample_fpe = choose_dynamic_factor_model_in_sample(x, FPE; max_factors=10, max_factor_lags=4)
#dynamic_model_in_sample_fpe_icp2 = choose_dynamic_factor_model_in_sample(x, FPE, "ICp2"; max_factors=10, max_factor_lags=4)
#dynamic_model_in_sample_icp2 = choose_dynamic_factor_model_in_sample(x, FPE, "ICp2"; bai_criterion=true, max_factors=10, max_factor_lags=4)  # if bai_criterion is true information_criterion plays no role
# Dynamic model out-of-sample model selection
#dynamic_model_out_of_sample = choose_dynamic_factor_model_out_of_sample(x; max_factors=10, max_factor_lags=4)
#dynamic_model_out_of_sample_ICp2 = choose_dynamic_factor_model_out_of_sample(x, "ICp2"; max_factors=10, max_factor_lags=4)





# bootstrapping the tests:
#residual_bootstrap(fm, 100, fm->LM_test(fm, t, i))
#tests = [quandt_andrews(fm, LM_test_gls, i) for i in 1:size(x, 2)]



# Structural break tests
function sup_test(fm::FactorModel, significance_level=0.05)
    tests = [quandt_andrews(fm, LM_test, i) for i in 1:size(fm.x, 2)]  # TODO: should use LM_test_gls here once it is bug free
    periods_tested = sort(tests[1][1])
    critical_value = quandt_andrews_critical_value(fm.number_of_factors, significance_level)
    highest_statistic = [test[2][1] for test in tests]  # extract most likely break_date statistic (aka the period with highest test statistic)
    highest_statistic_at = [test[1][1] for test in tests]  # extract most likely break_date (aka the period with highest test statistic)
    significant_break_indices = highest_statistic_at[highest_statistic.>critical_value]
    break_period_per_variable = highest_statistic_at .* (highest_statistic.>critical_value)
    #break_counts = OrderedDict(periods_tested, [count(stat->stat==t, highest_statistic_at) for t in periods_tested])
    significant_break_counts = OrderedDict(periods_tested, [count(stat->stat==i, significant_break_indices) for i in periods_tested])
    break_period_per_variable, significant_break_counts, periods_tested
end


using Gadfly
fm = FactorModel(x, "ICp2")
break_period_per_variable, significant_break_counts, periods_tested = sup_test(fm)

#plt = plot(
#    layer(x=periods_tested, y=[i for i in values(significant_break_counts)], Geom.point),
#    layer(x=periods_tested, y=[i for i in values(significant_break_counts)], Geom.line),
#    layer(xintercept=[14, 69], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=8)
#)
#draw(PNG("graphs/structural_breaks.png", 20cm, 15cm), plt)

# which variables have a break in period 25 and periods 62 to 65 which is when most variables have breaks
series_idx_with_break_at25 = find(Bool[break_period_per_variable[i].==25 for i in 1:size(fm.x, 2)])
series_id_with_break_at25 = header[series_idx_with_break_at25]
series_idx_with_break_62_to_65 = find(Bool[65 .>= break_period_per_variable[i] .>= 62 for i in 1:size(fm.x, 2)])
series_id_with_break_62_to_65 = header[series_idx_with_break_62_to_65]

# as an answer to the structural breaks we change the window of observations considered  (period 47 is 2004-07-01, period 61 is 2008-01-01)
# alternativel we can exclude the series with breaks around 2008
# we do this to that we have something to compare

x_window = x[:, Bool[!in(series_id, series_idx_with_break_62_to_65) for series_id in 1:size(x, 2)]]
fm_window = FactorModel(x_window, "ICp2")
break_period_per_variable_window, significant_break_counts_window, periods_tested_window = sup_test(fm_window, 0.05)

#plt = plot(
#    layer(x=periods_tested_window, y=[i for i in values(significant_break_counts_window)], Geom.point),
#    layer(x=periods_tested_window, y=[i for i in values(significant_break_counts_window)], Geom.line),
#    layer(xintercept=[minimum(periods_tested_window), maximum(periods_tested_window)], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=10)
#)
#draw(PNG("graphs/structural_breaks_less_variables_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with less variables
#dynamic_model_out_of_sample_ICp2_window = choose_dynamic_factor_model_out_of_sample(x_window, "ICp2")

#dynamic_model_out_of_sample_window = choose_dynamic_factor_model_out_of_sample(x_window, "ICp2")
#dynamic_model_out_of_sample_ICp2_window = choose_dynamic_factor_model_out_of_sample(x_window)


x_window2 = x[:, Bool[!in(series_id, series_idx_with_break_62_to_65) && !in(series_id, series_idx_with_break_at25) for series_id in 1:size(x, 2)]]
fm_window2 = FactorModel(x_window2, "ICp2")
break_period_per_variable_window2, significant_break_counts_window2, periods_tested_window2 = sup_test(fm_window2, 0.05)

#plt = plot(
#    layer(x=periods_tested_window2, y=[i for i in values(significant_break_counts_window2)], Geom.point),
#    layer(x=periods_tested_window2, y=[i for i in values(significant_break_counts_window2)], Geom.line),
#    layer(xintercept=[minimum(periods_tested_window2), maximum(periods_tested_window2)], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=10)
#)
#draw(PNG("graphs/structural_breaks_even_less_variables_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
#dynamic_model_out_of_sample_window2 = choose_dynamic_factor_model_out_of_sample(x_window2)
#dynamic_model_out_of_sample_ICp2_window2 = choose_dynamic_factor_model_out_of_sample(x_window2, "ICp2")


x_window3 = x[:, break_period_per_variable.==0]
fm_window3 = FactorModel(x_window3, "ICp2")
break_period_per_variable_window3, significant_break_counts_window3, periods_tested_window3 = sup_test(fm_window3, 0.05)

#plt = plot(
#    layer(x=periods_tested_window3, y=[i for i in values(significant_break_counts_window3)], Geom.point),
#    layer(x=periods_tested_window3, y=[i for i in values(significant_break_counts_window3)], Geom.line),
#    layer(xintercept=[minimum(periods_tested_window3), maximum(periods_tested_window3)], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=10)
#)
#draw(PNG("graphs/structural_breaks_all_breaks_removed_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
#dynamic_model_out_of_sample_window3 = choose_dynamic_factor_model_out_of_sample(x_window3)
#dynamic_model_out_of_sample_ICp2_window3 = choose_dynamic_factor_model_out_of_sample(x_window3, "ICp2")


x_window35 = x[:, break_period_per_variable.==0][:, break_period_per_variable_window3.>0]  # only testing if removing breaks again improves RMSE again
fm_window35 = FactorModel(x_window35, "ICp2")
break_period_per_variable_window35, significant_break_counts_window35, periods_tested_window35 = sup_test(fm_window35, 0.05)

#plt = plot(
#    layer(x=periods_tested_window35, y=[i for i in values(significant_break_counts_window35)], Geom.point),
#    layer(x=periods_tested_window35, y=[i for i in values(significant_break_counts_window35)], Geom.line),
#    layer(xintercept=[minimum(periods_tested_window35), maximum(periods_tested_window35)], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=10)
#)
#draw(PNG("graphs/structural_breaks_all_breaks_removed_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
#dynamic_model_out_of_sample_window35 = choose_dynamic_factor_model_out_of_sample(x_window35; max_factors=10, max_factor_lags=3, max_lags=10)
#dynamic_model_out_of_sample_ICp2_window35 = choose_dynamic_factor_model_out_of_sample(x_window35, "ICp2"; max_factors=10, max_factor_lags=3, max_lags=10)


x_window4 = x[26:end, Bool[!in(series_id, series_idx_with_break_62_to_65) for series_id in 1:size(x, 2)]]  # drop observations in and before first cluster and variables in second cluster
fm_window4 = FactorModel(x_window4, "ICp2")
break_period_per_variable_window4, significant_break_counts_window4, periods_tested_window4 = sup_test(fm_window4, 0.05)

#plt = plot(
#    layer(x=periods_tested_window4, y=[i for i in values(significant_break_counts_window4)], Geom.point),
#    layer(x=periods_tested_window4, y=[i for i in values(significant_break_counts_window4)], Geom.line),
#    layer(xintercept=[minimum(periods_tested_window4), maximum(periods_tested_window4)], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=10)
#)
#draw(PNG("graphs/structural_breaks_second_cluster_and_obs_removed.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
#dynamic_model_out_of_sample_window4 = choose_dynamic_factor_model_out_of_sample(x_window4)
#dynamic_model_out_of_sample_ICp2_window4 = choose_dynamic_factor_model_out_of_sample(x_window4, "ICp2")




# Robustnes check: same procedure but with predictors targeted for predicting y beforehand
# hard thresholding
x_targeted = x[:, targeted_predictors(1, x, 6, "hard", significance_level=0.01)]   # AR(6) is the best performing AR model so we use 6 lags here
fm_targeted = FactorModel(x_targeted, "ICp2")

#dynamic_model_out_of_sample_targeted = choose_dynamic_factor_model_out_of_sample(x_targeted)
#dynamic_model_out_of_sample_ICp2_targeted = choose_dynamic_factor_model_out_of_sample(x_targeted, "ICp2")

# soft thresholding. Be advised: this takes a while
soft_targeting_steps = 20:50
x_targeted_sets = [x[:, targeted_predictors(1, x, 6, "soft"; number_of_steps_in_lars=i)] for i in soft_targeting_steps]   # AR(6) is the best performing AR model so we use 6 lags here
fm_targeted_sets = [FactorModel(x_targeted, "ICp2") for x_targeted in x_targeted_sets]
x_targeted_sets = x_targeted_sets[find([fm.number_of_factors <= 7 for fm in fm_targeted_sets])]  # can only test breaks for less than 8 static factors
fm_targeted_sets = fm_targeted_sets[find([fm.number_of_factors <= 7 for fm in fm_targeted_sets])]  # can only test breaks for less than 8 static factors

breaks_period_per_variable_targeted_sets = [sup_test(fm_targeted, 0.05)[1] for fm_targeted in fm_targeted_sets]
breaks_period_per_variable_targeted_sets = [[0, breaks_period_per_variable_targeted[2:end]] for breaks_period_per_variable_targeted in breaks_period_per_variable_targeted_sets]  # first variables is not removed even if it has a break because its our variable of interest
#break_period_per_variable_targeted, significant_break_counts_targeted, periods_tested_targeted = sup_test(fm_targeted, 0.05)

#plt = plot(
#    layer(x=periods_tested, y=[i for i in values(significant_break_counts_targeted)], Geom.point),
#    layer(x=periods_tested, y=[i for i in values(significant_break_counts_targeted)], Geom.line),
#    layer(xintercept=[14, 69], Geom.vline(color="green")),
#    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
#    Scale.y_continuous(minvalue=0, maxvalue=8)
#)
#draw(PNG("graphs/structural_breaks_targeted_1percent.png", 20cm, 15cm), plt)


#dynamic_model_out_of_sample_targeted_sets = [choose_dynamic_factor_model_out_of_sample(x_targeted) for x_targeted in x_targeted_sets]
dynamic_model_out_of_sample_ICp2_targeted_sets = [choose_dynamic_factor_model_out_of_sample(x_targeted, "ICp2") for x_targeted in x_targeted_sets]

x_targeted_window_sets = [x_targeted_sets[i][:, breaks_period_per_variable_targeted_sets[i] .== 0] for i in 1:length(x_targeted_sets)]
# [size(x_targeted, 2) for x_targeted in x_targeted_sets]  # data sets are tiny!
# [size(x_targeted_window, 2) for x_targeted_window in x_targeted_window_sets]  # data sets are tiny!
dynamic_model_out_of_sample_ICp2_targeted_window_sets = [choose_dynamic_factor_model_out_of_sample(x_targeted_window, "ICp2"; max_factors=5, max_factor_lags=3, max_lags=5) for x_targeted_window in x_targeted_window_sets]
rmses_targeted_sets = [dynamic_model_out_of_sample_ICp2_targeted_sets[i][4] for i in 1:length(dynamic_model_out_of_sample_ICp2_targeted_sets)]
rmses_targeted_window_sets = [dynamic_model_out_of_sample_ICp2_targeted_window_sets[i][4] for i in 1:length(dynamic_model_out_of_sample_ICp2_targeted_window_sets)]
differences = rmses_targeted_sets - rmses_targeted_window_sets
forecasts_targeted_sets = [dynamic_model_out_of_sample_ICp2_targeted_sets[i][5] for i in 1:length(dynamic_model_out_of_sample_ICp2_targeted_sets)]
forecasts_targeted_window_sets = [dynamic_model_out_of_sample_ICp2_targeted_window_sets[i][5] for i in 1:length(dynamic_model_out_of_sample_ICp2_targeted_window_sets)]

#dynamic_model_out_of_sample_targeted = choose_dynamic_factor_model_out_of_sample(x_targeted_window)
#dynamic_model_out_of_sample_ICp2_targeted = choose_dynamic_factor_model_out_of_sample(x_targeted_window, "ICp2")


# now we can compare the two and calculate diebold mariano tests

diebold_results = [diebold_mariano(forecasts_targeted_sets[i][1], forecasts_targeted_window_sets[i][1], true_values, 0.05) for i in 1:length(forecasts_targeted_sets)]
p_values_one_sided = [db[2] for db in diebold_results]
diebold_results_two_sided = [diebold_mariano(forecasts_targeted_sets[i][1], forecasts_targeted_window_sets[i][1], true_values, 0.05; test_type="two sided") for i in 1:length(forecasts_targeted_sets)]
p_values_two_sided = [db[2] for db in diebold_results_two_sided]
dm_stats = [db[1] for db in diebold_results]
dm_rejections = [db[3] for db in diebold_results]
dm_rejections_two_sided = [db[3] for db in diebold_results_two_sided]
#diebold_results = [diebold_mariano(forecasts_targeted_window_sets[i][1], forecasts_targeted_sets[i][1], true_values, 0.10) for i in 1:length(forecasts_targeted_sets)]
#diebold_results = [diebold_mariano(forecasts_targeted_sets[i][1], forecasts_targeted_window_sets[i][1], true_values, 0.10; test_type="two sided") for i in 1:length(forecasts_targeted_sets)]
#diebold_mariano(bench_ar_forecasts, forecasts_targeted_sets[1][1], true_values)
#critical_value_improved_accuracy = quantile(TDist(number_of_forecasts-1), 0.90)
#critical_value_equal_accuracy = quantile(TDist(number_of_forecasts-1), 0.95)
import FactorModels.matrix_to_table
rnd(x) = round(x, 4)  # save me some letters typing
latex = matrix_to_table(convert(Array{Float64, 2}, hcat([int(size(x_targeted, 2)) for x_targeted in x_targeted_sets], [int(size(x_targeted_window, 2)) for x_targeted_window in x_targeted_window_sets], rnd(rmses_targeted_sets), rnd(rmses_targeted_window_sets), rnd(differences), rnd(dm_stats), rnd(p_values_one_sided), rnd(p_values_two_sided))))
#file = open("table.tex", "w")
#write(file, latex)
#close(file)


# compare x_window3 which results from removing all breaks from the original set and x_targeted_sets
x_window3_indices = break_period_per_variable.==0  # after all breaks have been removed
x_targeted_indices_sets = [targeted_predictors(1, x, 6, "soft"; number_of_steps_in_lars=i) for i in soft_targeting_steps]  # set of indices after soft thresholding by number of steps in lasso

same_indices = [[find(idx) in x_targeted_indices for idx in find(x_window3_indices)] for x_targeted_indices in x_targeted_indices_sets]

