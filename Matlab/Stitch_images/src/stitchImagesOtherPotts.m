% Author: Murat Apishev (great-mel@yandex.ru)
function [resultImage, resultMask] = stitchImagesOtherPotts(images, seeds)

% get number of labels
K = size(images, 1);
% get number of verticies
[N, M, ~] = size(images{1});

% prepare unary potentials
unary = zeros(N, M, K);
mask = zeros(N, M);
for i = 1 : K
    mask = mask | seeds{i};
end
for i = 1 : K
    unary(:, :, i) = xor(mask, seeds{i});
end
unary(unary ~= 0) = 1e+6 * unary(unary ~= 0);

% set coefficint c_ij as differencies between values of pixels of images
% in img_gray_scales scale
img_gray_scales = zeros(N, M, K);
for i = 1 : K
    img_gray_scales(:, :, i) = rgb2gray(images{i});
end
delta_gray = max(img_gray_scales, [], 3) - min(img_gray_scales, [], 3);
vertC = delta_gray(1 : (N - 1), :);
horC  = delta_gray(:, 1 : (M - 1));

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