
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
    rmses = [RMSE(prediction_tuples[number_of_factors, number_of_lags, number_of_factor_lags][1], prediction_tuples[number_of_factors, number_of_lags, number_of_factor_lags][2]) for number_of_factors in 1:max_factors, number_of_lags in 1:max_lags, number_of_factor_lags in 1:max_factor_lags]
    number_of_factors, number_of_lags, number_of_factor_lags = ind2sub(size(rmses), indmin(rmses))
    forecasts, true_values = prediction_tuples[number_of_factors, number_of_lags, number_of_factor_lags]
    return (number_of_factors, number_of_lags, number_of_factor_lags, minimum(rmses), forecasts, true_values)
end
