% BIEN 462 Question 4 Group Component

% 1 to 1 Input/Output system 
% Input: Arterial Blood Pressure (ABP) in mmHg
% Output: Cerebral Blood Flow Velocity (CBFV) in cm/s
% Sampling frequency of 1Hz -> 1 data point per sec
close all;

%% Part i)
bloodPressure = load('pres_data.txt');
flowVelocity = load('flow_data.txt');

bloodPressure_d = detrend(bloodPressure);
flowVelocity_d = detrend(flowVelocity);

% Segment for first 6 minutes
bp1 = bloodPressure_d(1:360); 
fv1 = flowVelocity_d(1:360);

% Segment for 30-36 minute interval
bp2 = bloodPressure_d(1800:2160);
fv2 = flowVelocity_d(1800:2160);

% Power Spectral Densities
[pbp, fbp] = pwelch(bloodPressure_d, [],[],[],1);
[pbp1, fbp1] = pwelch(bp1, [],[],[],1);
[pbp2, fbp2] = pwelch(bp2, [],[],[],1);
[pfv, ffv] = pwelch(flowVelocity_d, [],[],[],1);
[pfv1, ffv1] = pwelch(fv1, [],[],[],1);
[pfv2, ffv2] = pwelch(fv2, [],[],[],1);

figure
subplot(3,2,1)
hold on
title('Power Spectral Density of  Blood Pressure')
plot(fbp, pbp);
xlabel('Frequency (Hz)');
ylabel('PSD of Blood Pressure')
hold off

subplot(3,2,2)
hold on
title('Power Spectral Density of Flow Velocity')
plot(ffv, pfv);
xlabel('Frequency (Hz)');
ylabel('PSD of Flow Velocity')
hold off

subplot(3,2,3)
hold on
title('Power Spectral Density of  Blood Pressure for Time 0 to 6 mins')
plot(fbp1, pbp1);
xlabel('Frequency (Hz)');
ylabel('PSD of Blood Pressure')
hold off

subplot(3,2,4)
hold on
title('Power Spectral Density of Flow Velocity for Time 0 to 6 mins')
plot(ffv1, pfv1);
xlabel('Frequency (Hz)');
ylabel('PSD of Flow Velocity')
hold off

subplot(3,2,5)
hold on
title('Power Spectral Density of  Blood Pressure for Time 30 to 36 mins')
plot(fbp2, pbp2);
xlabel('Frequency (Hz)');
ylabel('PSD of Blood Pressure')
hold off

subplot(3,2,6)
hold on
title('Power Spectral Density of Flow Velocity for Time 30 to 36 mins')
plot(ffv2, pfv2);
xlabel('Frequency (Hz)');
ylabel('PSD of Flow Velocity')
hold off

% Comparing the full dataset case pwelch vs periodogram methods
[perio_pbp, perio_fbp] = periodogram(bloodPressure_d, [],[],1);
[perio_pfv, perio_ffv] = periodogram(flowVelocity_d, [],[],1);

figure
subplot(2,2,1)
hold on
title('PSD of  Blood Pressure with Welch')
plot(fbp, pbp)
xlabel('Frequency (Hz)');
ylabel('PSD of Blood Pressure')
hold off

subplot(2,2,2)
hold on
title('PSD of  Blood Pressure with Periodogram ')
plot(perio_fbp, perio_pbp)
xlabel('Frequency (Hz)');
ylabel('PSD of Blood Pressure')
hold off

subplot(2,2,3)
hold on
title('PSD of Flow Velocity with Welch')
plot(ffv, pfv)
xlabel('Frequency (Hz)');
ylabel('PSD of Flow Velocity')
hold off


subplot(2,2,4)
hold on
title('PSD of Flow Velocity with Periodogram')
plot(perio_ffv, perio_pfv);
xlabel('Frequency (Hz)');
ylabel('PSD of Flow Velocity')
hold off

%% Part ii)

% Coherence between both signals for the 3 time intervals
[c, f] = mscohere(bloodPressure_d, flowVelocity_d, [], [], [], 1);
[c1, f1] = mscohere(bp1, fv1, [], [], [], 1);
[c2, f2] = mscohere(bp2, fv2, [], [], [], 1);

figure
subplot(3,1,1)
hold on 
plot(f, c)
title('Coherence between ABP and CBFV')
xlabel('Frequency (Hz)')
ylabel('Coherence')
hold off

subplot(3,1,2)
hold on 
plot(f1, c1)
title('Coherence between ABP and CBFV for Time 0 to 6 mins')
xlabel('Frequency (Hz)')
ylabel('Coherence')
hold off

subplot(3,1,3)
hold on 
plot(f2, c2)
title('Coherence between ABP and CBFV for Time 30 to 36 mins')
xlabel('Frequency (Hz)')
ylabel('Coherence')
hold off


%% Q4 c

bloodPressure_d = detrend(bloodPressure(:));
flowVelocity_d = detrend(flowVelocity(:));

% Segment for first 6 minutes
bp1 = bloodPressure_d(1:360); 
fv1 = flowVelocity_d(1:360);

% Segment for 30-36 minute interval
bp2 = bloodPressure_d(1800:2160);
fv2 = flowVelocity_d(1800:2160);

% Entire segment of data (making sure of same length)
bp = bloodPressure_d(1:2398);
fv = flowVelocity_d(1:2398);

data = iddata(fv, bp, 1); % 2400s is the sampling time

% Estimate the impulse response model, using impulseest
sys = impulseest(data);

% Plotting the sytem's impulse response 

figure;
[impulseResponse, time] = impulse(sys);
plot(time*60, impulseResponse)

title('Estimated Impulse Response');
xlabel('Time(s)');
ylabel('Amplitude');

% Plotting frequency response:
figure; 
bode(sys);
title('Frequency Response of Estimated Impulse Response');

%% Q4 d

%% Q4 d

[H, f] = tfestimate(bp, fv, [], [], [], 1); % The last argument is the sampling frequency

% Plotting magnitude response:
figure;
f_rad_s = f * 2 * pi;
semilogx(f_rad_s, 20*log10(abs(H)));
title('System Magnitude Response');
xlabel('Frequency (rad/s)');
ylabel('Magnitude (dB)');


% Plotting phase response:
figure;
H_phase_deg = angle(H)*180 / pi;
semilogx(f_rad_s, H_phase_deg);

title("System Phase Response");
xlabel("Frequency (rad/s)");
ylabel('Phase (deg)');









