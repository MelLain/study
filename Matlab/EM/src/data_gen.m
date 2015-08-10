a = 2;
b = 2;
alpha = 0.01;
K = 20;
no_features = 50;
no_objects = 1000;

X = zeros(no_objects, no_features);
pi = gamrnd(repmat(alpha, K, 1), 1, K, 1);
pi = pi ./ sum(pi);
pi(pi < 1e-5) = 0;
mu = betarnd(a, b, K, no_features);

indices = rand(no_objects);
cum_pi = [0 cumsum(pi)'];

for n = 1 : no_objects
    for i = 1 : length(pi)
        if (indices(n) > cum_pi(i)) && (indices(n) <= cum_pi(i + 1))
            X(n,:) = binornd(1, mu(i,:));
            break;
        end
    end
end