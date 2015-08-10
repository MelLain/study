% This function implements EM-algorithm. All details, as well as
% descriptions of input and output parameters, can be found in
% docs/task.pdf file.
%
% Author: Murat Apishev (great-mel@yandex.ru)

function [pi, mu, L, rho_nk_last] = EM(X, options)
[no_objects, no_features] = size(X);
% set option variables default values
a = 1;
b = 1;
alpha = 0.001;
K = 50;
max_iter = 500;
tol = 1e-3;
no_start = 1;

% try to set option variables user defined values
if (isfield(options, 'a'))
    a = options.a;
end
if (isfield(options, 'b'))
    b = options.b;
end
if (isfield(options, 'alpha'))
    alpha = options.alpha;
end
if (isfield(options, 'K'))
    K = options.K;
end
if (isfield(options, 'max_iter'))
    max_iter = options.max_iter;
end
if (isfield(options, 'tol'))
    tol = options.tol;
end
if (isfield(options, 'n_start'))
    no_start = options.n_start;
end

best_value = -Inf;
best_mu = zeros(K, no_features);
best_pi = zeros(K, 1);
best_L = [];

best_rho_nk = [];

for attempt = 1 : no_start
    % generate random rho_nk
    rho_nk =  rand(no_objects, K);
    rho_nk = rho_nk ./ repmat(sum(rho_nk, 2), 1, K);

    cur_iter = 1;
    cur_values_diff = +Inf;
    cur_value = -Inf;
    no_iters_to_count = 10;

    mu = zeros(K, no_features);
    pi = repmat(1 / K, K, 1);
    L = [];
    
    while (cur_iter <= max_iter && cur_values_diff > tol)
        % proceed E-step
        hat_a = repmat(a, size(mu)) + rho_nk' * X;
        hat_b = repmat(b, size(mu)) + rho_nk' * (1 - X);
        mu = hat_a ./ (hat_a + hat_b);
        log_mu = psi(hat_a) - psi(hat_a + hat_b);
        log_1_minus_mu = psi(hat_b) - psi(hat_a + hat_b);

        log_pi = log(pi);
        log_pi(log_pi == -Inf) = 0;
        under_exp_with_pi = bsxfun(@plus, log_pi', X * log_mu' + (1 - X) * log_1_minus_mu');
        pre_rho_nk = exp(under_exp_with_pi);
        rho_nk = bsxfun(@rdivide, pre_rho_nk, sum(pre_rho_nk, 2));

        % proceed M-step
        pi = (sum(rho_nk)' + alpha - 1) ./ (no_objects + K * (alpha - 1));
        
        if (sum(pi(pi < 0)) ~= 0)
            pi(pi < 0) = 0;
            pi = pi ./ sum(pi);
        end
        non_zero_pi = pi > 0;
        pi = pi(non_zero_pi);
        pi = pi ./ sum(pi);
        rho_nk = rho_nk(:, non_zero_pi);
        rho_nk = rho_nk ./ repmat(sum(rho_nk, 2), 1, length(pi));
        mu = mu(non_zero_pi, :);
        
        % count value of functional L(q)
        if (mod(cur_iter, no_iters_to_count) == 0)
            log_mu = log_mu(non_zero_pi, :);
            log_1_minus_mu = log_1_minus_mu(non_zero_pi, :);
            hat_a = hat_a(non_zero_pi, :);
            hat_b = hat_b(non_zero_pi, :);
            log_pi = log(pi);
            log_pi(log_pi == -Inf) = 0; 
            under_exp_with_pi = bsxfun(@plus, log_pi', X * log_mu' + (1 - X) * log_1_minus_mu');
            part_1 = sum(sum(under_exp_with_pi .* rho_nk));
            part_2 = log(gamma(K * alpha));
            part_3 = - K * log(gamma(alpha));
            part_4 = sum((alpha - 1) .* log_pi);
            part_5 = - K * no_features * betaln(a, b);
            if (part_5 == -Inf)
                part_5 = 0;
            end
            part_6 = sum(sum((a - 1) .* log_mu + (b - 1) .* log_1_minus_mu));
            part_7 = - sum(sum(rho_nk(rho_nk > 0) .* log(rho_nk(rho_nk > 0))));
            part_8 = - sum(sum((hat_a - 1) .* log_mu + (hat_b - 1) .* log_1_minus_mu)) + ...
                sum(sum(betaln(hat_a, hat_b)));
            new_value = part_1 + part_2 + part_3 + part_4 + part_5 + part_6 + part_7 + part_8;            

            L(end + 1) = new_value;
            
            % if you'd like to replace absolute error with comparative one,
            % just uncomment next string and comment one after it
            %cur_values_diff = abs(new_value - cur_value) / abs(new_value);
            cur_values_diff = abs(new_value - cur_value);
            cur_value = new_value;
        end
        cur_iter = cur_iter + 1;
    end
    
    if (cur_value > best_value)
        best_value = cur_value;
        disp('last attempt iter count')
        disp(cur_iter - 1);
        disp('new best value')
        disp(cur_value)
        best_pi = pi;
        best_mu = mu;
        best_L = L;
        
        best_rho_nk = rho_nk;
    end
end

pi = [best_pi; zeros(K - length(best_pi), 1)];
mu = [best_mu; repmat(a / (a + b), K - length(best_pi), no_features)];
L = best_L;

rho_nk_last = best_rho_nk;

end

