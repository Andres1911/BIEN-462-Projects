% BIEN 462 Assignment 3 part ii)
clear;
time_step = 0.001;
tspan = 0:time_step:300;
n = length(tspan);

a = [0.02, 0.02, 0.1, 0.1];
b = [0.2, 0.2, 0.2, 0.25];
c = [-65, -55, -65, -65];
d = [8, 4, 2, 2];
current = 20;
V = -70*ones(1,n);
U = [zeros(1,n)];

currentTypes = {'impulse', 'random', 'sinusoidal', 'ramp'};

for i = 1:length(a)
    subplot(2,2,i)
    [tspan, V] = DiscreteModel(a(i), b(i), c(i), d(i), current, currentTypes(4)); % change number index for currentType 
    % and generate plots with altered titles
    plot(tspan, V)
    xlabel('Time (ms)')
    ylabel('Membrane potential (mV)')
    if i == 1
        title('Regular Spiking, Ramp Input Type')
    elseif i == 2
        title('Intrinsically Bursting, Ramp Input Type')
    elseif i == 3
        title('Fast Spiking, Ramp Input Type')
    else 
        title('Low-Threshold Spiking, Ramp Input Type')
    end
    hold off
end




function [tspan, V] = DiscreteModel(a, b, c, d, current, currentType)

time_step = 0.001;
tspan = 0:time_step:300;
n = length(tspan);
V = -70*ones(1,n); % Initial membrane potential
U = zeros(1,n); % Initial recovery variable

% Adjust current based on type

if strcmp(currentType,'impulse')
    I = zeros(1,n);
    I(10000:10010) = 20; % Short impulsive input
elseif strcmp(currentType,'random')
    I = -0.1 + (0.2).*rand(1,n); % Random values
elseif strcmp(currentType, 'sinusoidal')
    I = 5 * sin(2 * pi * (1/60) * (1:n)); % Sinusoidal input
elseif strcmp(currentType,'ramp')
    I = linspace(0, 20, n); % Ramp input
else
    I = zeros(1,n); % Default to no input
end

% Simulation loop
for i = 2:n
    V(i) = V(i - 1) + time_step * (0.04 * V(i - 1)^2 + 5 * V(i - 1) + 140 - U(i - 1) + I(i - 1));
    U(i) = U(i - 1) + time_step * (a * (b * V(i - 1) - U(i - 1)));

    if V(i) >= 30
        V(i) = c;
        U(i) = U(i) + d;
    end
end

end
