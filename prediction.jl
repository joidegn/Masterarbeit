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
T = size(data_matrix, 1) - 4  # we include 4 lags
y = data_matrix[:,1]  # TODO: this is not something we actually want to predict...
data_matrix = data_matrix[:, 2:end]
lag1 = lag_vector(y)
lag2 = lag_vector(lag1)
lag3 = lag_vector(lag2)
lag4 = lag_vector(lag3)
y = y[5:end]
w = hcat(ones(T), array(lag1[5:end]), array(lag2[5:end]), array(lag3[5:end]), array(lag4[5:end]))
x = data_matrix[5:end, :]

function pseudo_out_of_sample_forecasts(y_index=1, num_predictions::Int=10)
    # one step ahead pseudo out-of-sample forecasts
    T = size(x,1)
    predictions = zeros(num_predictions)
    to_predict = T-num_predictions+1:T
    for date_index in to_predict
        month = date_index - T+num_predictions
        println("month: $month")
        let y=x[1:date_index, y_index], x=x[1:date_index, [1:size(x,2)].!=y_index]  # y and x are updated so its easier for humans to read the next lines
            let newx=x[end, :], y=y[1:end-1], x=x[1:end-1, :]  # pseudo-one step ahead --> we predict only the last value in y using only information before that
                res = FactorModel(x, "ICp1")
                predictions[month] = predict(res, y)
            end
        end
    end
    return(predictions, x[to_predict, y_index])
end

forecasts, true_values1 = pseudo_out_of_sample_forecasts(1, 100)
bench_forecasts, true_values2 = benchmark_forecasts(x, 1; num_predictions=100)
mse = MSE(forecasts, true_values1)
mse_bench = MSE(bench_forecasts, true_values2)

#predictions_frame = DataFrame(
#    predictions=vcat(y[end-199:end], predictions, predictions_targeted, predictions_ols, predictions_average),
#    method=split(^("true value,", 200) * ^("Static Factor Model,", 200) * ^("targeted Static Factor Model,", 200) * ^("OLS,", 200) * ^("Average,", 200), ",", 0, false)
#)
#set_default_plot_format(:png)
#display(predictions_plot)




