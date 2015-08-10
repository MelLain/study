% This function implements counting of log(L). All details, as 
% well as descriptions of input and output parameters, can be 
% found in docs/task.pdf file.
%
% Author: Murat Apishev (great-mel@yandex.ru)

function LL = log_likelihood(X, mu, pi)

% use log(exp(...)) to avoid counting errors
p = X * log(mu)' + (1 - X) * log(1 - mu)';
LL = sum(log(exp(p) * pi));

end

