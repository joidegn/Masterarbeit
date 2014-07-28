using FactorModels

R=50
bs = [0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1]
deltas = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
#bs = [0.1, 0.5]
#deltas = [0.5, 0.9]

T=50
N=20
r=1
number_of_forecasts = 10
include("choose_model.jl")

results = Array(Any, R)  # each results holds an (length(bs)) times (length(deltas)) array of result tuples
results_reduced = Array(Any, R)
diebold_results = Array(Any, R)
diebold_results_two_sided = Array(Any, R)

tic()
for rep in 1:R
    data_sets = [factor_model_DGP(T, N, r; model="single_break", b=b, delta=delta) for b in bs, delta in deltas]
    x_sets = [data_set[1] for data_set in data_sets]
    optimal_args = [choose_static_factor_model_out_of_sample(x, "ICp2"; max_factors=6, max_lags=5) for x in x_sets]
    rmses = Float64[args[3] for args in optimal_args]
    forecasts_sets, true_values_sets = [args[4] for args in optimal_args], [args[5] for args in optimal_args]
    reduced_x_sets = [x_set[:, 1:end-1] for x_set in x_sets]
    optimal_args_reduced = [choose_static_factor_model_out_of_sample(x, "ICp2"; max_factors=6, max_lags=5) for x in reduced_x_sets]
    rmses_reduced = Float64[optimal_args[3] for optimal_args in optimal_args_reduced]
    results[rep] = reshape(optimal_args, (length(bs), length(deltas)))
    results_reduced[rep] = reshape(optimal_args_reduced, (length(bs), length(deltas)))
end
toc()  # 50 repetitions take about 5100 seconds which is 85 minutes  --> 1000 repetitions would take about 45 hours, 500 would take ~22 hours

diebold_results = [[diebold_mariano(results[rep][b_ind, delta_ind][2], results_reduced[rep][b_ind, delta_ind][2], results[rep][b_ind, delta_ind][3], 0.1) for b_ind in 1:length(bs), delta_ind in 1:length(deltas)] for rep in 1:R]
diebold_results_two_sided = [[diebold_mariano(results[rep][b_ind, delta_ind][2], results_reduced[rep][b_ind, delta_ind][2], results[rep][b_ind, delta_ind][3], 0.1; test_type="two sided") for b_ind in 1:length(bs), delta_ind in 1:length(deltas)] for rep in 1:R]

rmses = Float64[results[rep][b_ind, delta_ind][1] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
rmses_reduced = Float64[results_reduced[rep][b_ind, delta_ind][1] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
diff_rmses = rmses - rmses_reduced
mean_diff = mean(diff_rmses, 3)[:, :, 1]



diebold_stats = Float64[diebold_results[rep][b_ind, delta_ind][1] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
diebold_rejections_improved = Bool[diebold_results[rep][b_ind, delta_ind][3] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
diebold_rejections_different = Bool[diebold_results_two_sided[rep][b_ind, delta_ind][3] for b_ind in 1:length(bs), delta_ind in 1:length(deltas), rep in 1:R]
avg_rejections_improved = sum(diebold_rejections_improved, 3)[:, :, 1]./R
avg_rejections_different = sum(diebold_rejections_different, 3)[:, :, 1]./R

file = open("table_single_break.tex", "w")
write(file, matrix_to_table(hcat(bs, avg_rejections_improved)))
close(file)
