% Author: Murat Apishev (great-mel@yandex.ru)
function [labels, energy, time] = alphaExpansionGridPotts(unary, vertC, horC, metric, options)

% get number of labels and verticies
[N, M, K] = size(unary);

% set default values to parameters
maxIter = 500;
display = false;
numStart = 1;
randOrder = false;

% check existence of 'options' variable and update parameters, if variable exists
if (nargin > 4)
  if isfield(options, 'maxIter')
    maxIter = options.maxIter;
  end
  if isfield(options, 'display')
    display = options.display;
  end
  if isfield(options, 'numStart')
    numStart = options.numStart;
  end
  if isfield(options, 'randOrder')
    randOrder = options.randOrder;
  end
end

% min difference between energy values on neighbour iterarions to convergence
eps = 1e-8;
% generate the initialization of labels
labels = randi(K, N, M);
% set default order of labels processing
order = 1 : K;

% variables for storing intermediate energy and time values
intermediate_energy = [];
intermediate_time = [];

% start processing
for attempt = 1 : numStart
  iter = 0;
  while iter <= maxIter
    % if it necessary, generate a random order of labels processing
    if randOrder
      order = randperm(K);
    end
    for index = 1 : K
      tic;
      % get current label
      alpha = order(index);
      iter = iter + 1;

      % Initialize terminal weights. The first column is psi_i(1) =
      % theta_i(alpha), and the second column is psi_i(0) = theta_i(x_i), so
      % for it we need to create a mask.
      x_i_mask = (repmat(labels, [1, 1, K]) == repmat(reshape(1 : K, [1, 1, K]), N, M));
      term_weights = [reshape(unary(:, :, alpha), N * M, 1), reshape(sum(unary .* x_i_mask, 3), N * M, 1)];

      % The next part of code is the reparameterization. At first let's get
      % the values of metric for all pairs of labels in 'labels' variable,
      % from the left to the right, from the top to the bottom (for 
      % horizontal edges).
      temp = (reshape([labels(:, 1), reshape(repmat(labels(:, 2 : (end - 1)), 2, 1), N, 2 * (M - 2)), ...
                      labels(:, end)]', 2, N * (M - 1)))';
      pairwise_dist = (reshape(metric(temp(:, 1) + (temp(:, 2) - 1) * K), M - 1, N))';

      % define a = psi_ij(0,0), b = psi_ij(1,1), c = psi_ij(0,1)
      a = reshape(horC .* pairwise_dist, [], 1);
      b = reshape(horC, [], 1) .* metric(alpha, alpha);
      c = reshape(horC, [], 1) .* metric(labels(:, 1 : (end - 1)), alpha);
      
      % proceed all three steps of reparameterization
      i_indices = 1 : (N * (M - 1));
      j_indices = i_indices + N;
      % psi_i(0) += a
      term_weights(i_indices, 2) = term_weights(i_indices, 2) + a;
      % psi_j(1) += c - a
      term_weights(j_indices, 1) = term_weights(j_indices, 1) + c - a;
      % psi_i(1) += b - c + a
      term_weights(i_indices, 1) = term_weights(i_indices, 1) + b - c + a;
      
      % prepare the only parameter for adding to psi_ij(1,0) for horizontal
      % edges (c + d - a - b)
      delta_horisontal = ...
                       + horC .* reshape(metric(labels(:, 1 : (end - 1)), alpha), N, M - 1) ...
                       + horC .* reshape(metric(alpha, labels(:, 2 : end)), N, M - 1) ... 
                       - horC .* pairwise_dist ...
                       - horC .* metric(alpha, alpha);

      % Let's now get the values of metric for all pairs of labels
      % in 'labels' variable, from the top to the bottom, from the left to
      % the right (for vertical edges).
      temp = (reshape([labels(1, :); (reshape((repmat(labels(2 : (end - 1), :), 1, 2))', M, 2 * (N - 2)))'; ...
                       labels(end, :)], 2, []))';
      pairwise_dist = reshape(metric(temp(:, 1) + (temp(:, 2) - 1) * K), N - 1, M);

      % define a = psi_ij(0,0), b = psi_ij(1,1), c = psi_ij(0,1)
      a = reshape(vertC .* pairwise_dist, [], 1);
      b = reshape(vertC, [], 1) .* metric(alpha, alpha);
      c = reshape(vertC, [], 1) .* metric(labels(1 : (end - 1), :), alpha);     

      % proceed all three steps of reparameterization
      i_indices = 1 : (N * M);
      i_indices = i_indices(mod(i_indices, N) > 0)';
      j_indices = i_indices + 1;
      % psi_i(0) += a
      term_weights(i_indices, 2) = term_weights(i_indices, 2) + a;
      % psi_j(1) += c - a
      term_weights(j_indices, 1) = term_weights(j_indices, 1) + c - a;
      % psi_i(1) += b - c + a
      term_weights(i_indices, 1) = term_weights(i_indices, 1) + b - c + a;

      % prepare the only parameter for adding to psi_ij(1,0) for vertical
      % edges (c + d - a - b)      
      delta_vertical = ...
                     + vertC .* reshape(metric(labels(1 : (end - 1), :), alpha), N - 1, M) ...
                     + vertC .* reshape(metric(alpha, labels(2 : end, :)), N - 1, M) ...
                     - vertC .* pairwise_dist ...
                     - vertC .* metric(alpha, alpha);

      % Now let's prepare 'edge_weights' variable (4 columns) for graph
      % cut. Here's the passage from 'GraphCut()' description:
      % ... edgeWeights(i, 3) connects node #edgeWeights(i, 1) to node #edgeWeights(i, 2)
  		%	edgeWeights(i, 4) connects node #edgeWeights(i, 2) to node #edgeWeights(i, 1) ...
      from_up = 1 : (N * M);
      from_up = from_up(mod(from_up, N) ~= 0)';
      from_down = (1 : (N * (M - 1)))';
      edge_weights = [from_up, from_up + 1, zeros(size(from_up)), reshape(delta_vertical, [], 1); ...
                      from_down, from_down + N, zeros(size(from_down)), reshape(delta_horisontal, [], 1)];
      
      % unary potential conversion with delta = min(psi_i(0), psi_i(1))
      delta = min(term_weights, [], 2);
      term_weights = term_weights - repmat(delta, 1, 2);
      % run graph cut to get alpha expansion
      [cut, labels_mask] = graphCutMex(term_weights, edge_weights);
      % apply the expanation to 'labels'
      labels(logical(reshape(labels_mask, N, M))) = alpha;
      
      % save current energy value
      intermediate_energy(end + 1) = cut + sum(delta);
      % display intermediate results, if necessary
      if display
        fprintf('Attempt #%d, Iteration #%d, Alpha label = %d, Energy = %f\n', attempt, ...
                iter, alpha, intermediate_energy(end));
      end
      % save current time value
      intermediate_time(end + 1) = toc;
    end
    
    if (all(abs(intermediate_energy((end - K + 1) : end) - intermediate_energy(end)) < eps))
      % the algorithm has converged
      break;
    end
  end
  
  % choose the final result due to best energy value
  if (attempt == 1) || (intermediate_energy(end) < energy(end))
    energy = intermediate_energy;
    time = intermediate_time;
  end
end
time = cumsum(time);
end