using FactorModels
using DataFrames


# Models:
# - AR(p)
# - Static FM with no breaks, number of factors according to ICp1
# - Static FM with 1 breaks, number of factors according to ICp1 on all subperiods
# ??- Static FM with more breaks?, number of factors according to ICp1 on all subperiods
# - Dynamic FM with no breaks, number of factors according to ICp1 AND number of factor lags according to Breitung, Eickmeier p. 75
# - Dynamic FM with breaks, number of factors according to ICp1 AND number of factor lags according to Breitung, Eickmeier p. 75 --> loop through and let both number or factors and factor lags converge?? (or get number of factors first and then factor lags)

data = readtable("/home/joi/Documents/Konstanz/Masterarbeit/data/1959-2014_normalized.csv")
data_matrix = reduce(hcat, [convert(Array{Float64}, col) for col in data.columns[2:size(data.columns, 1)]])
#ids = map(string, names(data)[2:end])
#titles = series_titles(ids)  # TODO: does not work at the moment because ICU.jl and with it Requests.jl seems to be broken
y_ind = 1
x = data_matrix  # x includes y

BIC(mse::Float64, n::Int64, T) = log(mse) + n * log(T)/T  # Bayesian information criterion
FPE(mse::Float64, n::Int64, T) = log(mse) + n * 2/T  # Fisher information criterion

function pseudo_out_of_sample_forecasts(y_index, number_of_lags::Int64, number_of_factors::Int64; num_predictions::Int=10)
    # one step ahead pseudo out-of-sample forecasts
    T = size(x,1)
    predictions = zeros(num_predictions)
    to_predict = T-num_predictions+1:T
    for date_index in to_predict
        month = date_index - T+num_predictions
        println("month: $month")
        let y=x[1:date_index, y_index], x=x[1:date_index, [1:size(x,2)].!=y_index]  # y and x are updated so its easier for humans to read the next lines
            let newx=x[end, :], y=y[1:end-1], x=x[1:end-1, :]  # pseudo-one step ahead --> we predict only the last value in y using only information before that
                res = FactorModel(hcat(x, y))
                predictions[month] = predict(res, y, 1, number_of_factors)[2]
            end
        end
    end
    return(predictions, x[to_predict, y_index])
end

function pseudo_out_of_sample_forecasts(y_index::Int64, number_of_lags::Int64, number_of_factors::Int64, number_of_factor_lags::Int64; num_predictions::Int64=10)
    # one step ahead pseudo out-of-sample forecasts
    T = size(x,1)
    predictions = zeros(num_predictions)
    to_predict = T-num_predictions+1:T
    for date_index in to_predict
        month = date_index - T+num_predictions
        println("month: $month")
        let y=x[1:date_index, y_index], x=x[1:date_index, [1:size(x,2)].!=y_index]  # y and x are updated so its easier for humans to read the next lines
            let newx=x[end, :], y=y[1:end-1], x=x[1:end-1, :]  # pseudo-one step ahead --> we predict only the last value in y using only information before that
                res = DynamicFactorModel((hcat(x, y),), number_of_factor_lags)
                predictions[month] = predict(res, y, 1, number_of_factors)[2]
            end
        end
    end
    return(predictions, x[to_predict, y_index])
end

# static factor model
residuals = [predict(FactorModel(x, number_of_factors), x[:, y_ind], 1, number_of_lags)[1] for number_of_factors in 1:int(ceil(size(x,2)/2)), number_of_lags in 1:10]
bics = [BIC(sum(residuals[number_of_factors, number_of_lags].^2), number_of_factors+number_of_lags, size(x, 2)-number_of_lags-1) for number_of_factors in 1:int(ceil(size(x,2)/2)), number_of_lags in 1:10]  # T is smaller because of lags and specificity of predict function
fpes = [FPE(sum(residuals[number_of_factors, number_of_lags].^2), number_of_factors+number_of_lags, size(x, 2)-number_of_lags-1) for number_of_factors in 1:int(ceil(size(x,2)/2)), number_of_lags in 1:10]
number_of_lags = int(ceil(indmin(bics)/size(bics, 1)))  # do indmin by hand because it doesnt seem to work with two-dimensional arrays
number_of_factors = int(indmin(bics)-(number_of_lags-1)*size(bics, 1))
forecasts, true_values1 = pseudo_out_of_sample_forecasts(1, number_of_lags, number_of_factors; num_predictions=100)
bench_forecasts, true_values2 = benchmark_forecasts(x, 1; num_predictions=100, window=100)
mse = MSE(forecasts, true_values1)
mse_bench = MSE(bench_forecasts, true_values2)

# dynamic factor model
residuals_dyn = [predict(DynamicFactorModel((x, number_of_factors), number_of_factor_lags), x[:, y_ind], 1, number_of_lags)[1] for number_of_factors in 1:int(ceil(size(x,2)/2)), number_of_lags in 1:10, number_of_factor_lags in 1:11]
bics_dyn = [BIC(sum(residuals[number_of_factors, number_of_lags].^2), number_of_factors+number_of_lags, size(x, 2)-number_of_lags-1) for number_of_factors in 1:int(ceil(size(x,2)/2)), number_of_lags in 1:10, number_of_factor_lags in 1:11]  # T is smaller because of lags and specificity of predict function
fpes_dyn = [FPE(sum(residuals[number_of_factors, number_of_lags].^2), number_of_factors+number_of_lags, size(x, 2)-number_of_lags-1) for number_of_factors in 1:int(ceil(size(x,2)/2)), number_of_lags in 1:10]
# do indmin by hand because it doesnt seem to work with 3-dimensional arrays
number_of_factor_lags_dyn = int(ceil(indmin(bics_dyn)/(size(bics_dyn, 1)*size(bics_dyn, 2))))  # third dimension we are in is index number divided by number of elements in first and second dimension combined (which is the number of elements "per third dimension")
number_of_lags_dyn = int(ceil((indmin(bics_dyn)-(number_of_factor_lags_dyn-1)*size(bics_dyn, 1)*size(bics_dyn, 2))/size(bics_dyn, 1)))  # second dimension is then (index number - elements left if we deduct the elements from the third dimension) / elements in first dimension
number_of_factors_dyn = int(indmin(bics_dyn)-(number_of_factor_lags_dyn-1)*size(bics_dyn, 2)*size(bics_dyn, 1)-(number_of_lags_dyn-1)*size(bics_dyn, 1))  # first dimension is then again what is left if we account for third and second
forecasts_dyn, true_values1 = pseudo_out_of_sample_forecasts(1, number_of_lags_dyn, number_of_factors_dyn, number_of_factor_lags_dyn; num_predictions=100)
mse_dyn = MSE(forecasts_dyn, true_values1)

#predictions_frame = DataFrame(
#    predictions=vcat(y[end-199:end], predictions, predictions_targeted, predictions_ols, predictions_average),
#    method=split(^("true value,", 200) * ^("Static Factor Model,", 200) * ^("targeted Static Factor Model,", 200) * ^("OLS,", 200) * ^("Average,", 200), ",", 0, false)
#)
#set_default_plot_format(:png)
#display(predictions_plot)




