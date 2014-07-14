using FactorModels
using DataFrames


# Models:
# - AR(p)
# - Static FM with no breaks, number of factors according to ICp1
# - Static FM with 1 breaks, number of factors according to ICp1 on all subperiods
# ??- Static FM with more breaks?, number of factors according to ICp1 on all subperiods
# - Dynamic FM with no breaks, number of factors according to ICp1 AND number of factor lags according to Breitung, Eickmeier p. 75
# - Dynamic FM with breaks, number of factors according to ICp1 AND number of factor lags according to Breitung, Eickmeier p. 75 --> loop through and let both number or factors and factor lags converge?? (or get number of factors first and then factor lags)

#data = readtable("data/final_data.csv", separator=' ', header=true)  # this is how it used to be
#data_matrix = reduce(hcat, [convert(Array{Float64}, col) for col in data.columns[2:size(data.columns, 1)]])

file = readdlm("data/final_data.csv")
data, header, dates = file[2:end, 2:end], file[1, 2:end], file[2:end, 1]
#ids = map(string, names(data)[2:end])
#titles = series_titles(ids)  # TODO: does not work at the moment because ICU.jl and with it Requests.jl seems to be broken
y_ind = 1
x = convert(Array{Float64, 2}, data)  # x includes y
y = data[:, 1]

BIC(mse::Float64, n::Int64, T) = log(mse) + n * log(T)/T  # Bayesian information criterion
FPE(mse::Float64, n::Int64, T) = log(mse) + n * 2/T  # Fisher information criterion

function pseudo_out_of_sample_forecasts(y_index::Int64, number_of_lags::Int64, number_of_factors::Int64, number_of_factors_criterion::String=""; num_predictions::Int=10)
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

function pseudo_out_of_sample_forecasts(y_index::Int64, number_of_lags::Int64, number_of_factors::Int64, number_of_factor_lags::Int64, number_of_factors_criterion::String=""; num_predictions::Int64=10)
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
                predictions[period] = predict(dfm, y, 1, number_of_factors)[2]
            end
        end
    end
    return(predictions, x[to_predict, y_index])
end

# static factor model
# we estimate

max_factors = 15
max_factor_lags = 3
max_lags = 10
number_of_forecasts = 15


# Benchmark MSE
bench_forecasts, true_values = benchmark_forecasts(x, 1; num_predictions=number_of_forecasts, window=number_of_forecasts)
rmse_bench = RMSE(bench_forecasts, true_values)

bench_ar_forecasts, true_values_ = benchmark_ar(x[:, 1]; num_predictions=number_of_forecasts)
rmse_bench_ar = RMSE(bench_ar_forecasts, true_values)



function choose_static_factor_model_in_sample(information_criterion::Function=BIC, number_of_factors_criterion::String=""; bai_criterion=false)  # if bai_criterion is true number_of_factors_criterion is used in forecasting equation as well
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
    rmse = apply(RMSE, pseudo_out_of_sample_forecasts(1, number_of_lags, number_of_factors; num_predictions=number_of_forecasts))
    return(number_of_factors, number_of_lags, rmse)
end
function choose_static_factor_model_out_of_sample(number_of_factors_criterion::String="")
    # choice of parameters according to out-of-sample criteria
    prediction_tuples = [pseudo_out_of_sample_forecasts(1, number_of_lags, number_of_factors; num_predictions=number_of_forecasts) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
    rmses = [RMSE(prediction_tuples[number_of_factors, number_of_lags][1], prediction_tuples[number_of_factors, number_of_lags][2]) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags]
    number_of_factors, number_of_lags = ind2sub(size(rmses), indmin(rmses))
    return(number_of_factors, number_of_lags, minimum(rmses))
end

# dynamic factor model
# in-sample choice
function choose_dynamic_factor_model_in_sample(information_criterion::Function=BIC, number_of_factors_criterion::String=""; bai_criterion=false)
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
    forecasts, true_values1 = pseudo_out_of_sample_forecasts(1, number_of_lags, number_of_factors, number_of_factor_lags; num_predictions=number_of_forecasts)
    rmse = RMSE(forecasts, true_values1)
    return (number_of_factors, number_of_lags, number_of_factor_lags, rmse)
end
function choose_dynamic_factor_model_out_of_sample(number_of_factors_criterion::String="")
    #out-of-sample choice
    prediction_tuples = [pseudo_out_of_sample_forecasts(1, number_of_lags, number_of_factors, number_of_factor_lags, number_of_factors_criterion; num_predictions=number_of_forecasts) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
    rmses = [RMSE(prediction_tuples[number_of_factors, number_of_lags, number_of_factor_lags][1], prediction_tuples[number_of_factors, number_of_lags][2]) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
    number_of_factors, number_of_lags, number_of_factor_lags = ind2sub(size(rmses), indmin(rmses))
    return (number_of_factors, number_of_lags, number_of_factor_lags, minimum(rmses))
end


# Static model in-sample model selection
static_model_in_sample_bic = choose_static_factor_model_in_sample(BIC)
static_model_in_sample_bic_icp2 = choose_static_factor_model_in_sample(BIC, "ICp2")
static_model_in_sample_fpe = choose_static_factor_model_in_sample(FPE)
static_model_in_sample_fpe_icp2 = choose_static_factor_model_in_sample(FPE, "ICp2")
static_model_in_sample_icp2 = choose_static_factor_model_in_sample(FPE, "ICp2"; bai_criterion=true)
# Static model out-of-sample model selection
static_model_out_of_sample = choose_static_factor_model_out_of_sample()
static_model_out_of_sample_ICp2 = choose_static_factor_model_out_of_sample("ICp2")


# Dynamic model in-sample model selection
dynamic_model_in_sample_bic = choose_dynamic_factor_model_in_sample(BIC)
dynamic_model_in_sample_bic_icp2 = choose_dynamic_factor_model_in_sample(BIC, "ICp2")
dynamic_model_in_sample_fpe = choose_dynamic_factor_model_in_sample(FPE)
dynamic_model_in_sample_fpe_icp2 = choose_dynamic_factor_model_in_sample(FPE, "ICp2")
dynamic_model_in_sample_icp2 = choose_dynamic_factor_model_in_sample(FPE, "ICp2"; bai_criterion=true)
# Dynamic model out-of-sample model selection
dynamic_model_out_of_sample = choose_dynamic_factor_model_out_of_sample()
dynamic_model_out_of_sample_ICp2 = choose_dynamic_factor_model_out_of_sample("ICp2")






#predictions_frame = DataFrame(
#    predictions=vcat(y[end-199:end], predictions, predictions_targeted, predictions_ols, predictions_average),
#    method=split(^("true value,", 200) * ^("Static Factor Model,", 200) * ^("targeted Static Factor Model,", 200) * ^("OLS,", 200) * ^("Average,", 200), ",", 0, false)
#)
#set_default_plot_format(:png)
#display(predictions_plot)




