% Practicum, Task #3, 'Compositions of algorithms'.
%
% FUNCTION:
% gradient_boosting_experiment (X, y, num_iterations, base_algorithm, loss, ...
%                                   model_complexity, learning_rate)
%
% DESCRIPTION:
% This function describes the experiment with gradient boosting and builds plots.
%
% INPUT:
% X --- matrix of objects, N x K double matrix, N --- number of objects, 
%       K --- number of features.
% y --- vector of answers, N x 1 double vector, N --- number of objects. y
%       can have only two values --- +1 and -1.
% num_iterations --- the number ob algorithms in composition, scalar.
% base_algorithm --- the base algorithm, string. Can have one of two
%                    values: 'regression_tree' or 'epsilon_svr'.
% loss --- the loss function, string. Can have one of two values: 
%          'logistic' (for classification) or 'absolute' (for regression).
% model_complexity --- int value in [1 3], discrete parameter, that allows 
% to control the complexity of model. 0 --- speshial case for stump usage.
% learning_rate --- parameter of gradient boosting, scalar.
% early_stopping --- parameter, that controls the early stop of train, bool.
%
% AUTHOR: 
% Murat Apishev (great-mel@yandex.ru)
%

function [] = gradient_boosting_experiment (X, y, num_iterations, base_algorithm, loss, ...
                                  model_complexity, learning_rate, early_stopping)

    switch model_complexity
        case 0 
            if strcmp(base_algorithm, 'regression_tree')
                min_parent = floor(size(X, 1) / 5) * 4 - 1;
            else
                error('Incorrect type of algorithm!');
            end            
        case 1
            if strcmp(base_algorithm, 'epsilon_svr')
                gamma = 0.1;
                C = 1;
                epsilon = 0.001;
            elseif strcmp(base_algorithm, 'regression_tree')
                min_parent = 50;
            else
                error('Incorrect type of algorithm!');
            end
        case 2
            if strcmp(base_algorithm, 'epsilon_svr')
                gamma = 0.001;
                C = 2;
                epsilon = 0.001;
            elseif strcmp(base_algorithm, 'regression_tree')
                min_parent = 15;
            else
                error('Incorrect type of algorithm!');
            end           
        case 3
            if strcmp(base_algorithm, 'epsilon_svr')
                gamma = 0.0001;
                C = 5;
                epsilon = 0.001;
            elseif strcmp(base_algorithm, 'regression_tree')
                min_parent = 5;
            else
                error('Incorrect type of algorithm!');
            end                  
        otherwise
            error('Incorrect complexity of model!');
    end   
    
    if early_stopping == true
        num_iterations = floor(num_iterations / 3);
    end

    no_objects = length(y);
    fifth = floor(no_objects / 5);
    error_train = zeros([num_iterations 1]);
    error_test = zeros([num_iterations 1]);
    
    for fold = 1 : 5
        mask_test = zeros([1 no_objects]);
        mask_test((fold - 1) * fifth + 1 : fifth * fold) = 1;
        mask_test = logical(mask_test);
        train_set = X(~mask_test,:);
        train_ans = y(~mask_test);
        test_set = X(mask_test,:);
        test_ans = y(mask_test);

        if strcmp(base_algorithm, 'epsilon_svr')
            model = gradient_boosting_train(train_set, train_ans, num_iterations, ...
                                        base_algorithm, loss, 'learning_rate', learning_rate, ...
                                        'epsilon', epsilon, 'gamma', gamma, 'C', C);
        elseif strcmp(base_algorithm, 'regression_tree')
            model = gradient_boosting_train(train_set, train_ans, num_iterations, ...
                                        base_algorithm, loss, 'learning_rate', learning_rate, ...
                                        'min_parent', min_parent);
        else
            error('Incorrect type of algorithm');
        end
        [~, error_train_loop] = gradient_boosting_predict(model, train_set, train_ans);
        [~, error_test_loop] = gradient_boosting_predict(model, test_set, test_ans);        
        error_train = error_train + error_train_loop;
        error_test = error_test + error_test_loop;
    end
    error_train = error_train / 5;
    error_test = error_test / 5;

    hold on;
    xlabel('Number of iterations');
    ylabel('Error');
%     title('Gradient boosting experiment'); 
    plot(1 : num_iterations, error_train, 'color', 'r', 'linewidth', 2);
    plot(1 : num_iterations, error_test, 'color', 'b', 'linewidth', 2);
    legend('Error on train', 'Error on test', 'Location', 'NorthEast');
    hold off;
end    