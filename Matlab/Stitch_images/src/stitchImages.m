% Author: Murat Apishev (great-mel@yandex.ru)
function [resultImage, resultMask] = stitchImages(images, seeds)

% get number of labels
K = size(images, 1);
% get number of verticies
[N, M, ~] = size(images{1});

% prepare unary potentials as euqlid distance to nearest granule
unary = zeros(N, M, K);
for label = 1 : K
  unary(:, :, label) = bwdist(seeds{label});
end

% set coefficint c_ij to ones
vertC = ones(N - 1, M);
horC = ones(N, M - 1);

% set metric to ones for different labels and to zeros for same ones
metric = 1 - eye(K);

% display the information about processing
% here new option cases can be added
options.display = true;
% run alpha-expansion algorithm
resultMask = alphaExpansionGridPotts(unary, vertC, horC, metric, options);

% compose final image using mask and base images
resultImage = uint8(zeros([N, M, 3]));
for alpha = 1 : K
  alpha_index = repmat((resultMask == alpha), [1, 1, 3]);
  resultImage(alpha_index) = images{alpha}(alpha_index);
end
end