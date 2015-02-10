function [] = opt_elastic_net_for_R()

% perfromes the optimization step
% output is vectors of parameters and minimized optimial value

% XX is a matrix of patients data
% yy is the response vector (+1 -1 for non-diseased and diseased patients)
% alpha_var controls the sparseness of the results
% lambda is the penalty calculated during the algorithm
% noise is a vector of noise values
% noise_scale is the multiplicative scale of the noise

load('cvx_data.mat');

MM = length(XX(1, :)) ; % number of features
NN = length(XX(:, 1)) ; % number of patients


% Solve optimization problem
cvx_begin
    variable bbeta(MM);
    minimize ((-yy(:)' * (XX * bbeta(:)) + sum_log(1 + exp( (XX * bbeta(:))'))) / NN  + lambda / 2 * (1 - alpha_var) * sum_square_abs(bbeta(2:MM, 1)) + lambda * alpha_var *  sum(abs(bbeta(2:end, 1))) + noise_scale * dot(noise(2:MM), bbeta(2:MM, 1)));
cvx_end

csvwrite('cvx_results.csv', bbeta)

end
