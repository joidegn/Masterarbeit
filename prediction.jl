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

function pseudo_out_of_sample_forecasts(x::Array{Float64, 2}, y_index::Int64, number_of_lags::Int64, number_of_factors::Int64, number_of_factors_criterion::String=""; num_predictions::Int=10)
    # one step ahead pseudo out-of-sample forecasts
    # the number of factors is used for the forecasting equation not for the factor equation if a criterion is given
    T = size(x,1)
    predictions = zeros(num_predictions)
    to_predict = T-num_predictions+1:T
    for date_index in to_predict
        period = date_index - T+num_predictions
        println("out-of-sample period: $period")
        let y=x[1:date_index, y_index], x=x[1:date_index, [1:size(x,2)].!=y_index]  # y and x are updated so its easier for humans to read the next lines
            let newx=x[end, :], y=y[1:end-1], x=x[1:end-1, :]  # pseudo-one step ahead --> we predict only the last value in y using only information before that
                length(number_of_factors_criterion) > 0 ? fm = FactorModel(hcat(y, x), number_of_factors_criterion) : fm = FactorModel(hcat(y, x), number_of_factors)
                predictions[period] = predict(fm, y, 1, number_of_factors)[2]
            end
        end
    end
    return(predictions, x[to_predict, y_index])
end

function pseudo_out_of_sample_forecasts(x::Array{Float64, 2}, y_index::Int64, number_of_lags::Int64, number_of_factors::Int64, number_of_factor_lags::Int64, number_of_factors_criterion::String=""; num_predictions::Int64=10)
    # one step ahead pseudo out-of-sample forecasts for a dynamic model
    # the number of factors is used for the forecasting equation not for the factor equation if a criterion is given
    T = size(x,1)
    predictions = zeros(num_predictions)
    to_predict = T-num_predictions+1:T
    for date_index in to_predict
        period = date_index - T+num_predictions
        println("out-of-sample period: $period")
        let y=x[1:date_index, y_index], x=x[1:date_index, [1:size(x,2)].!=y_index]  # y and x are updated so its easier for humans to read the next lines
            let newx=x[end, :], y=y[1:end-1], x=x[1:end-1, :]  # pseudo-one step ahead --> we predict only the last value in y using only information before that
                length(number_of_factors_criterion) > 0 ? dfm = DynamicFactorModel((hcat(y, x), number_of_factors_criterion,), number_of_factor_lags) : dfm = DynamicFactorModel((hcat(x, y), number_of_factors, ), number_of_factor_lags)
                predictions[period] = predict(dfm, y, 1, number_of_lags, number_of_factors)[2]
            end
        end
    end
    return(predictions, x[to_predict, y_index])
end

# static factor model
function choose_static_factor_model_in_sample(x::Array{Float64, 2}, information_criterion::Function=BIC, number_of_factors_criterion::String=""; bai_criterion=false, max_factors=15, max_lags=10)  # if bai_criterion is true number_of_factors_criterion is used in forecasting equation as well
    # choice of parameters according to in-sample criteria
    if length(number_of_factors_criterion) == 0
        residuals = [predict(FactorModel(x, number_of_factors), x[:, y_ind], 1, number_of_lags)[1] for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
    else
        if bai_criterion
            residuals = [predict(FactorModel(x, number_of_factors_criterion), x[:, y_ind], 1, number_of_lags)[1] for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
        else
            residuals = [predict(FactorModel(x, number_of_factors_criterion), x[:, y_ind], 1, number_of_lags, number_of_factors)[1] for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
        end
    end
    if bai_criterion
        number_of_factors = FactorModel(x, number_of_factors_criterion).number_of_factors
        info_criteria = [information_criterion(sum(residuals[number_of_factors, number_of_lags].^2), number_of_factors+number_of_lags+1, size(x, 2)-number_of_lags-1) for number_of_lags in 1:max_lags]  # T is smaller because of lags and specificity of predict function
        number_of_lags = indmin(info_criteria)
    else
        info_criteria = [information_criterion(sum(residuals[number_of_factors, number_of_lags].^2), number_of_factors+number_of_lags+1, size(x, 2)-number_of_lags-1) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]  # T is smaller because of lags and specificity of predict function
        number_of_factors, number_of_lags = ind2sub(size(info_criteria), indmin(info_criteria))
    end
    rmse = apply(RMSE, pseudo_out_of_sample_forecasts(x, 1, number_of_lags, number_of_factors; num_predictions=number_of_forecasts))
    return(number_of_factors, number_of_lags, rmse)
end
function choose_static_factor_model_out_of_sample(x::Array{Float64, 2}, number_of_factors_criterion::String=""; max_factors=15, max_lags=10) 
    # choice of parameters according to out-of-sample criteria
    prediction_tuples = [pseudo_out_of_sample_forecasts(x, 1, number_of_lags, number_of_factors; num_predictions=number_of_forecasts) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
    rmses = [RMSE(prediction_tuples[number_of_factors, number_of_lags][1], prediction_tuples[number_of_factors, number_of_lags][2]) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
    number_of_factors, number_of_lags = ind2sub(size(rmses), indmin(rmses))
    return(number_of_factors, number_of_lags, minimum(rmses))
end

# dynamic factor model
# in-sample choice
function choose_dynamic_factor_model_in_sample(x::Array{Float64, 2}, information_criterion::Function=BIC, number_of_factors_criterion::String=""; bai_criterion=false, max_factors=15, max_factor_lags=3, max_lags=10)
    if length(number_of_factors_criterion) == 0
        residuals = [predict(DynamicFactorModel((x, number_of_factors), number_of_factor_lags), x[:, y_ind], 1, number_of_lags, number_of_factors)[1] for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
    else
        if bai_criterion
            residuals = [predict(DynamicFactorModel((x, number_of_factors_criterion), number_of_factor_lags), x[:, y_ind], 1, number_of_lags)[1] for number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
        else
            residuals = [predict(DynamicFactorModel((x, number_of_factors_criterion), number_of_factor_lags), x[:, y_ind], 1, number_of_lags, number_of_factors)[1] for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
        end
    end
    if bai_criterion
        number_of_factors = FactorModel(x, number_of_factors_criterion).number_of_factors  # number of factors in forecating equation is derrived using bai criterion and not via BIC criterion
        info_criteria = [information_criterion(sum(residuals[number_of_lags, number_of_factor_lags].^2), number_of_factors*(number_of_factor_lags+1)+number_of_lags, size(x, 2)-maximum([number_of_lags, number_of_factor_lags])-1) for number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
        number_of_lags, number_of_factor_lags = ind2sub(size(info_criteria), indmin(info_criteria))
    else
        info_criteria = [information_criterion(sum(residuals[number_of_factors, number_of_lags, number_of_factor_lags].^2), number_of_factors*(number_of_factor_lags+1)+number_of_lags, size(x, 2)-maximum([number_of_lags, number_of_factor_lags])-1) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]  # T is smaller because of lags and specificity of predict function
        number_of_factors, number_of_lags, number_of_factor_lags = ind2sub(size(info_criteria), indmin(info_criteria))
    end
    #number_of_factors, number_of_lags, number_of_factor_lags = ind2sub(size(fpes), indmin(fpes))
    forecasts, true_values1 = pseudo_out_of_sample_forecasts(x, 1, number_of_lags, number_of_factors, number_of_factor_lags; num_predictions=number_of_forecasts)
    rmse = RMSE(forecasts, true_values1)
    return (number_of_factors, number_of_lags, number_of_factor_lags, rmse)
end
function choose_dynamic_factor_model_out_of_sample(x::Array{Float64, 2}, number_of_factors_criterion::String=""; max_factors=15, max_factor_lags=3, max_lags=10)
    #out-of-sample choice
    prediction_tuples = [pseudo_out_of_sample_forecasts(x, 1, number_of_lags, number_of_factors, number_of_factor_lags, number_of_factors_criterion; num_predictions=number_of_forecasts) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
    rmses = [RMSE(prediction_tuples[number_of_factors, number_of_lags, number_of_factor_lags][1], prediction_tuples[number_of_factors, number_of_lags][2]) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
    number_of_factors, number_of_lags, number_of_factor_lags = ind2sub(size(rmses), indmin(rmses))
    return (number_of_factors, number_of_lags, number_of_factor_lags, minimum(rmses))
end





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
static_model_in_sample_bic = choose_static_factor_model_in_sample(x, BIC)
static_model_in_sample_bic_icp2 = choose_static_factor_model_in_sample(x, BIC, "ICp2")
static_model_in_sample_fpe = choose_static_factor_model_in_sample(x, FPE)
static_model_in_sample_fpe_icp2 = choose_static_factor_model_in_sample(x, FPE, "ICp2")
static_model_in_sample_icp2 = choose_static_factor_model_in_sample(x, FPE, "ICp2"; bai_criterion=true, max_lags=3)  # if bai_criterion is true information criterion (FPE) plays no role
# Static model out-of-sample model selection
static_model_out_of_sample = choose_static_factor_model_out_of_sample(x)
static_model_out_of_sample_icp2 = choose_static_factor_model_out_of_sample(x, "ICp2")


# Dynamic model in-sample model selection
dynamic_model_in_sample_bic = choose_dynamic_factor_model_in_sample(x, BIC; max_factors=10, max_factor_lags=4)  # if max factors is higher residuals will be tiny (non distibguishable from 0, clear overfit)
dynamic_model_in_sample_bic_icp2 = choose_dynamic_factor_model_in_sample(x, BIC, "ICp2"; max_factors=10, max_factor_lags=4)
dynamic_model_in_sample_fpe = choose_dynamic_factor_model_in_sample(x, FPE; max_factors=10, max_factor_lags=4)
dynamic_model_in_sample_fpe_icp2 = choose_dynamic_factor_model_in_sample(x, FPE, "ICp2"; max_factors=10, max_factor_lags=4)
dynamic_model_in_sample_icp2 = choose_dynamic_factor_model_in_sample(x, FPE, "ICp2"; bai_criterion=true, max_factors=10, max_factor_lags=4)  # if bai_criterion is true information_criterion plays no role
# Dynamic model out-of-sample model selection
dynamic_model_out_of_sample = choose_dynamic_factor_model_out_of_sample(x; max_factors=10, max_factor_lags=4)
dynamic_model_out_of_sample_ICp2 = choose_dynamic_factor_model_out_of_sample(x, "ICp2"; max_factors=10, max_factor_lags=4)





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

plt = plot(
    layer(x=periods_tested, y=[i for i in values(significant_break_counts)], Geom.point),
    layer(x=periods_tested, y=[i for i in values(significant_break_counts)], Geom.line),
    layer(xintercept=[14, 69], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=8)
)
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

plt = plot(
    layer(x=periods_tested_window, y=[i for i in values(significant_break_counts_window)], Geom.point),
    layer(x=periods_tested_window, y=[i for i in values(significant_break_counts_window)], Geom.line),
    layer(xintercept=[minimum(periods_tested_window), maximum(periods_tested_window)], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=10)
)
#draw(PNG("graphs/structural_breaks_less_variables_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with less variables
dynamic_model_out_of_sample_ICp2_window = choose_dynamic_factor_model_out_of_sample(x_window, "ICp2")

dynamic_model_out_of_sample_window = choose_dynamic_factor_model_out_of_sample(x_window, "ICp2")
dynamic_model_out_of_sample_ICp2_window = choose_dynamic_factor_model_out_of_sample(x_window)


x_window2 = x[:, Bool[!in(series_id, series_idx_with_break_62_to_65) && !in(series_id, series_idx_with_break_at25) for series_id in 1:size(x, 2)]]
fm_window2 = FactorModel(x_window2, "ICp2")
break_period_per_variable_window2, significant_break_counts_window2, periods_tested_window2 = sup_test(fm_window2, 0.05)

plt = plot(
    layer(x=periods_tested_window2, y=[i for i in values(significant_break_counts_window2)], Geom.point),
    layer(x=periods_tested_window2, y=[i for i in values(significant_break_counts_window2)], Geom.line),
    layer(xintercept=[minimum(periods_tested_window2), maximum(periods_tested_window2)], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=10)
)
#draw(PNG("graphs/structural_breaks_even_less_variables_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
dynamic_model_out_of_sample_window2 = choose_dynamic_factor_model_out_of_sample(x_window2)
dynamic_model_out_of_sample_ICp2_window2 = choose_dynamic_factor_model_out_of_sample(x_window2, "ICp2")


x_window3 = x[:, break_period_per_variable.==0]
fm_window3 = FactorModel(x_window3, "ICp2")
break_period_per_variable_window3, significant_break_counts_window3, periods_tested_window3 = sup_test(fm_window3, 0.05)

plt = plot(
    layer(x=periods_tested_window3, y=[i for i in values(significant_break_counts_window3)], Geom.point),
    layer(x=periods_tested_window3, y=[i for i in values(significant_break_counts_window3)], Geom.line),
    layer(xintercept=[minimum(periods_tested_window3), maximum(periods_tested_window3)], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=10)
)
#draw(PNG("graphs/structural_breaks_all_breaks_removed_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
dynamic_model_out_of_sample_window3 = choose_dynamic_factor_model_out_of_sample(x_window3)
dynamic_model_out_of_sample_ICp2_window3 = choose_dynamic_factor_model_out_of_sample(x_window3, "ICp2")


x_window35 = x[:, break_period_per_variable.==0][:, break_period_per_variable_window3.>0]  # only testing if removing breaks again improves RMSE again
fm_window35 = FactorModel(x_window35, "ICp2")
break_period_per_variable_window35, significant_break_counts_window35, periods_tested_window35 = sup_test(fm_window35, 0.05)

plt = plot(
    layer(x=periods_tested_window35, y=[i for i in values(significant_break_counts_window35)], Geom.point),
    layer(x=periods_tested_window35, y=[i for i in values(significant_break_counts_window35)], Geom.line),
    layer(xintercept=[minimum(periods_tested_window35), maximum(periods_tested_window35)], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=10)
)
#draw(PNG("graphs/structural_breaks_all_breaks_removed_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
dynamic_model_out_of_sample_window35 = choose_dynamic_factor_model_out_of_sample(x_window35; max_factors=10, max_factor_lags=3, max_lags=10)
dynamic_model_out_of_sample_ICp2_window35 = choose_dynamic_factor_model_out_of_sample(x_window35, "ICp2"; max_factors=10, max_factor_lags=3, max_lags=10)


x_window4 = x[26:end, Bool[!in(series_id, series_idx_with_break_62_to_65) for series_id in 1:size(x, 2)]]  # drop observations in and before first cluster and variables in second cluster
fm_window4 = FactorModel(x_window4, "ICp2")
break_period_per_variable_window4, significant_break_counts_window4, periods_tested_window4 = sup_test(fm_window4, 0.05)

plt = plot(
    layer(x=periods_tested_window4, y=[i for i in values(significant_break_counts_window4)], Geom.point),
    layer(x=periods_tested_window4, y=[i for i in values(significant_break_counts_window4)], Geom.line),
    layer(xintercept=[minimum(periods_tested_window4), maximum(periods_tested_window4)], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=10)
)
#draw(PNG("graphs/structural_breaks_second_cluster_and_obs_removed.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with even less variables
dynamic_model_out_of_sample_window4 = choose_dynamic_factor_model_out_of_sample(x_window4)
dynamic_model_out_of_sample_ICp2_window4 = choose_dynamic_factor_model_out_of_sample(x_window4, "ICp2")




# Robustnes check: same procedure but with predictors targeted for predicting y beforehand
targeted_predictors(1, x, 6; thresholding="hard")

fm = FactorModel(x, "ICp2")
break_period_per_variable, significant_break_counts, periods_tested = sup_test(fm)

plt = plot(
    layer(x=periods_tested, y=[i for i in values(significant_break_counts)], Geom.point),
    layer(x=periods_tested, y=[i for i in values(significant_break_counts)], Geom.line),
    layer(xintercept=[14, 69], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=8)
)
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

plt = plot(
    layer(x=periods_tested_window, y=[i for i in values(significant_break_counts_window)], Geom.point),
    layer(x=periods_tested_window, y=[i for i in values(significant_break_counts_window)], Geom.line),
    layer(xintercept=[minimum(periods_tested_window), maximum(periods_tested_window)], Geom.vline(color="green")),
    Guide.XLabel("period"), Guide.YLabel("number of breaks"),
    Scale.y_continuous(minvalue=0, maxvalue=10)
)
#draw(PNG("graphs/structural_breaks_less_variables_1percent.png", 20cm, 15cm), plt)
# redo the out-of-sample forecasts with less variables
dynamic_model_out_of_sample_ICp2_window = choose_dynamic_factor_model_out_of_sample(x_window, "ICp2")

dynamic_model_out_of_sample_window = choose_dynamic_factor_model_out_of_sample(x_window, "ICp2")
dynamic_model_out_of_sample_ICp2_window = choose_dynamic_factor_model_out_of_sample(x_window)

