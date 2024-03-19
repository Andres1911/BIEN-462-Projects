% BIEN 462 Assignment 3 part i)
clear;
time_step = 0.001;
tspan = 0:time_step:200;

a = [0.02, 0.02, 0.1, 0.1];
b = [0.2, 0.2, 0.2, 0.25];
c = [-65, -55, -65, -65];
d = [8, 4, 2, 2];
current = 20;


for i = 1:length(a)

    subplot(2,2, i)
    [tspan, V] = DiscreteModel(a(i), b(i), c(i), d(i), current);
    
    hold on 
    plot(tspan, V)
    xlabel('Time (ms)')
    ylabel('Membrane potential (mV)')
    
    if i == 1
        title('Regular Spiking')
    elseif i == 2
        title('Intrinsically Bursting')
    elseif i == 3
        title('Fast Spiking')
    else 
        title('Low-Threshold Spiking')
    end
    hold off
end


function [tspan, V] = DiscreteModel(a, b, c, d, current)

time_step = 0.001;
tspan = 0:time_step:150;

pulseWidth = 0.9;
currentAmplitude = current;
n = length(tspan);
first = (1-pulseWidth)*n;

I = zeros(1,n);
I(first:end) = currentAmplitude;

V = -70*ones(1,n);
U = [zeros(1,n)];


for i = 2:length(tspan)
    V(i) = V(i - 1) + time_step*(0.04*V(i - 1)^2 + 5*V(i - 1) + 140 - U(i - 1) + I(i - 1));
    U(i) = U(i - 1) + time_step*(a*(b*V(i - 1) - U(i - 1)));

    if V(i) >= 30
        V(i) = c;
        U(i) = U(i) + d;
    end 
end

plot(tspan, V)

end 