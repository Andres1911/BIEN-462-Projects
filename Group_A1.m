%% Part i)

N = 250;

x = linspace(-1,1,N);

%The signal

y = sin(pi*x);

% Adding noise to the signal

t1 = zeros(1, size(x,2));

t2 = zeros(1, size(x,2));

t3 = zeros(1, size(x,2));

t4 = zeros(1, size(x,2));

for point = 1:size(x,2)

    t1(point) = y(point) + 0.10*randn(1); % sigma = 0.10

    t2(point) = y(point) + 0.20*randn(1); % sigma = 0.15

    t3(point) = y(point) + 0.30*randn(1); % sigma = 0.20

    t4(point) = y(point) + 0.40*randn(1); % sigma = 0.25

end

figure

hold on

plot(x,y,'k')

plot(x,t1,'b')

plot(x,t2,'g')

plot(x,t3,'m')

plot(x,t4,'c')

legend('Noise-free', 'Noise (\sigma = 0.10)', 'Noise (\sigma = 0.20)', 'Noise (\sigma = 0.30)', 'Noise (\sigma = 0.40)', 'Location','SouthEast')

title('Data from f(x) = sin(\pi x) with N = 250 observations')

hold off









%% Part ii)

% Fitting a polynomial to the t1 (sine signal with noise sigma = 0.25)

t_all = [t1; t2; t3; t4];

error_all = zeros(4,10);

sigmas = [10 15 20 25];

for iter = 1:4

    t = t_all(iter, :);

    M = 1:1:10;

    p = zeros(10,10);

    y = zeros(10, N);

    error = zeros(1,10);

    % Generating the polynomials for orders 1 to 10

    for n = 1:10

        p_temp = polyfit(x,t,n);

        i = 1;

        for coeff = p_temp

            p(n, i) = coeff;

            i = i + 1;

        end

        % Evaluating the polynomial

        y(n,:) = polyval(p(n,1:n+1), x);

        % Calculating the least squares error for the fit

        error(n) = sum((y(n,:) - t(:,:)).^2);

    end

    % Plotting all polynomial fits with the original signal data

    subplot(2,2,iter)

    hold on

    plot(x,t,'o')

    for n = 1:10

        plot(x,y(n,:))

    end

    str = sprintf('Polynomial Fit to Sine Function data with noise (sigma = 0.%d)', fix(sigmas(iter)));

    title(str)

    hold off

    error_all(iter, :) = error;

end



% Plotting the error to complexity graph

figure

for iter = 1:4

    subplot(2,2,iter)

    hold on

    plot(M, error_all(iter, :))

    str = sprintf('Error to Model Complexity (Polynomial Order M) (sigma = 0.%d)', fix(sigmas(iter)));

    title(str)

    xlabel('Complexity (M)')

    ylabel('Least-Squares Error')

    hold off

end

%% Part 3

testData = zeros(1, length(x));
f = sin(pi*x);

for i = 1:length(x)

    testData(i) = f(i) + 0.4*randn(1);

end

figure;
plot(x, testData);
sumOfSquares = zeros(1,length(M));

for i = 1:length(M)

    sumOfSquares(1,i) = sum((y(i,:) - testData(1,:)).^2);

end


figure;
plot(M, sumOfSquares)
hold on
plot(M, error_all(4,:))
legend('Test Error', 'Training Error')
title('Error to Model Complexity (Polynomial Order M) (\sigma = 0.25)')
xlabel('Complexity (M)')
ylabel('Least-Squares Error')

%% Part 4
clear;

N = 250; % Data Points

p = 3; % pol order

variance = 0.15; % variance in the data

lambdaAll = linspace(0,3,20); % Create array of 20 lamdas for testing

errorLambdas = zeros(1, length(lambdaAll)); % error for the training data
errorTest = zeros(1, length(lambdaAll)); % error for the test data

for i = 1:length(lambdaAll)

    [errorLambdas(i), errorTest(i)] = regularization(N, p, lambdaAll(i), variance); % Insert error in vectors

end

figure; % Plot vectors (error vs log(lambda)

logLambdas = log(lambdaAll);

hold on
plot(logLambdas, errorLambdas);
plot(logLambdas, errorTest);

title('Obtained Error vs different regularization coefficients')
ylabel('Error')
xlabel('ln \lambda')
legend('Training', 'Test')
hold off


[s, v, x, t, y] = regularization(N, p, exp(-1), variance);

% Given lambda = 1/e, plot the graph

figure;
hold on
plot(x, t, 'o')
plot(x, y)
title('Linear Fit of t = sin(\pi x), 3rd degree polynomial, \lambda = 1/e')

hold off

function [sumOfSquares, valSquares, x, testData, y] = regularization(N, p, lambda, variance) % Returns the sum of squares of the training data and the sum of sqaures of testing data

x = linspace(-1,1,N);

testData = zeros(1, N);

f = sin(pi*x);  % Generate Sin function

for i = 1:N

    testData (i) = f(i) + sqrt(variance)*randn(1);  % Add noise to sin function

end

ridgeMatrix = ones(N, p + 1); % Build regression matrix

for j = 1:p

    ridgeMatrix(:, j+1) = x.^j; % Generate values of x

end

w = (lambda*eye(p + 1) + ridgeMatrix'* ridgeMatrix)\ridgeMatrix'*testData'; % solve given formula

y = ridgeMatrix * w; % find data points for model function
y = y';

sumOfSquares = sum((y(:) - testData(:)).^2); % find sum of sqaures (expected - actual)

% Build validation data

valData = zeros(1, N);

f = sin(pi*x);  % Generate Sin function

for i = 1:N

   valData (i) = f(i) + sqrt(variance)*randn(1);  % Add noise to sin function

end

valSquares = sum((y(:) - valData(:)).^2);

end









