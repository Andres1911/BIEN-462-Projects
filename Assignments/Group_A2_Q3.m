%% Question 3a
clear;
x1 = -sqrt(12) + 2 * sqrt(12)*rand(1,5000); % Uniformly distributed white noise
x2 = randn(1,5000); %Gaussian uniformly distributed white noise with unit variance

[r1,lag] = xcorr(x1,64, 'unbiased');
[r2,lag] = xcorr(x2,64, 'unbiased');

figure
subplot(1,2,1)
stem(lag,r1)
title("Autocorrelation of a uniformly distributed random signal")
xlabel('lag')


subplot(1,2,2)
stem(lag,r2)
title("Autocorrelation of a normally distributed random signal")
xlabel('lag')

% The theoretical autocorrelation of a white noise signal is equal to the variance squared times the dirac delta function.
% For the first signal, the variance was found to be ~ 4, while the second has a variance of 1. 
% Accordingly, the peak of the first function is four times greater than the second. 

%% Question 3b
%H(z) = bz/z-a or b/1-az^-1
% Numerator coefficients: b
% Denominator coefficients: 1, -a

b = 0.3;
a = 0.8;

y1 = filter(b,[1, -a],x1);
y2 = filter(b,[1, -a],x2);


r3 = xcorr(y1,64, 'unbiased');
r4 = xcorr(y2,64, 'unbiased');

% How do they compare?

n = -64:64;
f = zeros(1, length(n));
variance = 1;

for i = 1:length(n)
    if n(i) < 0
        f(1, i) = 0.25*variance*1.25^(n(i));
    else
        f(1, i) = 0.25*variance*0.8^(n(i));
    end
end
f = f';


n = -64:64;
f2 = zeros(1, length(n));
variance = 4;

for i = 1:length(n)
    if n(i) < 0
        f2(1, i) = 0.25*variance*1.25^(n(i));
    else
        f2(1, i) = 0.25*variance*0.8^(n(i));
    end
end
f2 = f2';

figure;
subplot(2,2,1)
stem(lag,r3)
title("Autocorrelation of filtered uniformly distributed random signal, \sigma = 4")
xlabel('lag')

subplot(2,2,2)
stem(lag,r4)
title("Autocorrelation of filtered normally dist. random signal, \sigma = 1")
xlabel('lag')

subplot(2,2,3)
stem(n,f2)
title("True autocorrelation, \sigma = 4")
xlabel('lag')

subplot(2,2,4)
stem(n,f)
title("True autocorrelation, \sigma = 1")
xlabel('lag')


%% Question 3c

figure 
subplot(2,1,1)
histogram(y1)

subplot(2,1,2)
histogram(y2)
