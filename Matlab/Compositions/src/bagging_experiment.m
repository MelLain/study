% Practicum, Task #3, 'Compositions of algorithms'.
%
% FUNCTION:
% bagging_experiment (X, y, num_iterations, base_algorithm, model_complexity)
%
% DESCRIPTION:
% This function describes the experiment with bagging and builds plots.
%
% INPUT:
% X --- matrix of objects, N x K double matrix, N --- number of objects, 
%       K --- number of features.
% y --- vector of answers, N x 1 double vector, N --- number of objects. y
%       can have only two values --- +1 and -1.
% num_iterations --- the number ob algorithms in composition, scalar.
% base_algorithm --- the base algorithm, string. Can have one of two
%                    values: 'classification_tree' or 'svm'.
% model_complexity --- int value in [1 3], discrete parameter, that allows 
% to control the complexity of model.
%
% AUTHOR: 
% Murat Apishev (great-mel@yandex.ru)
%

function [] = bagging_experiment (X, y, num_iterations, base_algorithm, model_complexity)    

    switch model_complexity
        case 1
            if strcmp(base_algorithm, 'svm')
                gamma = 0.1;
                C = 1;
            elseif strcmp(base_algorithm, 'classification_tree')
                min_parent = 20;
            else
                error('Incorrect type of algorithm!');
            end
        case 2
            if strcmp(base_algorithm, 'svm')
                gamma = 0.01;
                C = 2;
            elseif strcmp(base_algorithm, 'classification_tree')
                min_parent = 10;
            else
                error('Incorrect type of algorithm!');
            end           
        case 3
            if strcmp(base_algorithm, 'svm')
                gamma = 0.001;
                C = 10;
            elseif strcmp(base_algorithm, 'classification_tree')
                min_parent = 2;
            else
                error('Incorrect type of algorithm!');
            end                   
        otherwise
            error('Incorrect complexity of model!');
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

        if strcmp(base_algorithm, 'svm')
            model = bagging_train(train_set, train_ans, num_iterations, base_algorithm, ...
                                'gamma', gamma, 'C', C);
        elseif strcmp(base_algorithm, 'classification_tree')
            model = bagging_train(train_set, train_ans, num_iterations, base_algorithm, ...
                                'min_parent', min_parent);            
        else
            error('Incorrect type of algorithm');
        end
        [~, error_train_loop] = bagging_predict(model, train_set, train_ans);
        [~, error_test_loop] = bagging_predict(model, test_set, test_ans);        
        error_train = error_train + error_train_loop;
        error_test = error_test + error_test_loop;
    end
    error_train = error_train / 5;
    error_test = error_test / 5;

    hold on;
    xlabel('Number of iterations');
    ylabel('Error');
%     title('Bagging experiment'); 
    plot(1 : num_iterations, error_train, 'color', 'r', 'linewidth', 2);
    plot(1 : num_iterations, error_test, 'color', 'b', 'linewidth', 2);
    legend('Error on train', 'Error on test', 'Location', 'NorthEast');
    hold off;
end         