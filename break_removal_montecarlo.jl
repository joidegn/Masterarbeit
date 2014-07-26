using FactorModels

R=10
bs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
deltas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
#bs = [0.1, 0.3, 0.5, 0.7, 1]
#deltas = [0.1, 02, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]

T=50
N=20
r=1
number_of_forecasts = 10
include("choose_model.jl")

results = Array(Any, R)  # each results holds an (length(bs)) times (length(deltas)) array of result tuples
results_reduced = Array(Any, R)
diebold_results = Array(Any, R)
diebold_results_two_sided = Array(Any, R)

for rep in 1:R
    data_sets = [factor_model_DGP(T, N, r; model="single_break", b=b, delta=delta) for b in bs, delta in deltas]
    x_sets = [data_set[1] for data_set in data_sets]
    optimal_args = [choose_dynamic_factor_model_out_of_sample(x, "ICp2"; max_factors=6, max_factor_lags=3, max_lags=5) for x in x_sets]
    rmses = Float64[args[4] for args in optimal_args]
    forecasts_sets, true_values_sets = [args[5] for args in optimal_args], [args[6] for args in optimal_args]
    reduced_x_sets = [x_set[:, 1:end-1] for x_set in x_sets]
    optimal_args_reduced = [choose_dynamic_factor_model_out_of_sample(x, "ICp2"; max_factors=6, max_factor_lags=3, max_lags=5) for x in reduced_x_sets]
    rmses_reduced = Float64[optimal_args[4] for optimal_args in optimal_args_reduced]
    forecasts_reduced_sets, true_values_reduced_sets = [args[5] for args in optimal_args_reduced], [args[6] for args in optimal_args_reduced]
    results[rep] = reshape([(rmses[i], forecasts_sets[i], true_values_sets[i]) for i in 1:length(data_sets)], (length(bs), length(deltas)))
    results_reduced[rep] = reshape([(rmses_reduced[i], forecasts_reduced_sets[i], true_values_reduced_sets[i]) for i in 1:length(data_sets)], (length(bs), length(deltas)))  # we have lost the shape above for convenience
    diebold_results[rep] = reshape([diebold_mariano(forecasts_sets[i], forecasts_reduced_sets[i], true_values_sets[i]) for i in 1:length(data_sets)], (length(bs), length(deltas)))
    diebold_results_two_sided[rep] = reshape([diebold_mariano(forecasts_sets[i], forecasts_reduced_sets[i], true_values_sets[i]; test_type="two sided") for i in 1:length(data_sets)], (length(bs), length(deltas)))
end

rmses = Float64[results[rep][b_ind, delta_ind][1] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
rmses_reduced = Float64[results_reduced[rep][b_ind, delta_ind][1] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
diff_rmses = rmses - rmses_reduced
mean_diff = mean(diff_rmses, 3)[:, :, 1]

diebold_rejections_improved = Bool[diebold_results[rep][b_ind, delta_ind][3] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
diebold_rejections_different = Bool[diebold_results[rep][b_ind, delta_ind][3] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
