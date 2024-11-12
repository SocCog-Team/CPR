% Parameters
numFrames = 200;            % Total number of frames
numDots = 200;              % Number of dots
coherence = 0.5;              % Coherence level (e.g., 50%)
coherentDirection = pi/4;   % Nominal direction of coherent motion (e.g., 45 degrees)
frameRate = 100;            % Frame rate in Hz (100 Hz)

% Pre-allocate arrays for motion energy and directional deviations
motionEnergy = zeros(1, numFrames);
directionDeviation = zeros(1, numFrames);


% Initialize motion energy array to store values for each frame
motionEnergy = zeros(1, numFrames);
overallMotionEnergy = zeros(1, numFrames); % Overall motion energy

for frame = 1:numFrames
    % Generate random directions for all dots
    randomDirections = 2 * pi * rand(numDots, 1);
    
    % Assign coherence fraction
    coherentDots = round(coherence * numDots);
    randomDots = numDots - coherentDots;
    
    % Set coherent direction for a fraction of dots
    randomDirections(1:coherentDots) = coherentDirection;
    
    % Calculate x and y velocities for all dots
    xVelocities = cos(randomDirections);
    yVelocities = sin(randomDirections);
    
    % Calculate net x and y movement for the frame
    netX = sum(xVelocities);
    netY = sum(yVelocities);
    
    % Compute actual motion direction for this frame
    actualDirection = atan2(netY, netX);
    
    % Calculate signed deviation from the coherent direction in radians
    directionDeviation(frame) = actualDirection - coherentDirection;
    
    % Wrap deviation to the range [-pi, pi] for meaningful interpretation
    if directionDeviation(frame) > pi
        directionDeviation(frame) = directionDeviation(frame) - 2 * pi;
    elseif directionDeviation(frame) < -pi
        directionDeviation(frame) = directionDeviation(frame) + 2 * pi;
    end
    
    % Convert deviation to degrees
    directionDeviation(frame) = rad2deg(directionDeviation(frame));
    
    % Motion energy for this frame (projecting onto coherent direction)
    % This is the dot product of (netX, netY) with the coherent direction vector
    motionEnergy(frame) = netX * cos(coherentDirection) + netY * sin(coherentDirection);
    
    
    % Calculate overall motion energy as the magnitude of the net movement vector
    % (without projecting onto any specific direction)
    overallMotionEnergy(frame) = sqrt(netX^2 + netY^2);
    
    
end


% Plot all measures of fluctuation

figure;

% 1. Plot the signed directional deviation in degrees over time
subplot(2,1,1);
plot((1:numFrames) / frameRate, directionDeviation, 'b');
xlabel('Time (s)');
ylabel('Deviation (degrees)');
title('Signed Directional Deviations Over Time');
grid on

% 2. Plot the motion energy over time
subplot(2,1,2);
plot((1:numFrames) / frameRate, motionEnergy, 'b'); hold on
plot((1:numFrames) / frameRate, overallMotionEnergy, 'r:'); hold on
legend('Coherent','Overall');

xlabel('Time (s)');
ylabel('Motion Energy');
title('Motion Energy Over Time');
grid on

% % Power spectrum for projected motion energy
% figure;
% subplot(2,1,1);
% pspectrum(motionEnergy-mean(motionEnergy), frameRate);  % frameRate is the sampling frequency in Hz
% title('Power Spectrum of Projected Motion Energy');
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');
%
% % Power spectrum for overall motion energy
% subplot(2,1,2);
% pspectrum(overallMotionEnergy-mean(overallMotionEnergy), frameRate);  % frameRate is the sampling frequency in Hz
% title('Power Spectrum of Overall Motion Energy');
% xlabel('Frequency (Hz)');
% ylabel('Power/Frequency (dB/Hz)');


N = length(motionEnergy);

% Calculate the Fourier Transform of the motion energy signal
motionEnergyFFT = fft(motionEnergy-mean(motionEnergy));
[motionEnergySpectrum,f] = pspectrum(motionEnergy-mean(motionEnergy),frameRate);

frequencies = (0:N-1) * (frameRate / N); % Frequency axis in Hz

% Limit to Nyquist frequency (first half of the FFT results)
halfN = floor(N / 2);
validFrequencies = frequencies(1:halfN);             % Frequencies up to Nyquist

% Plot the power spectrum
figure;
plot(validFrequencies, abs(motionEnergyFFT(1:halfN)/N).^2); hold on;
plot(f, abs(motionEnergySpectrum),'r');
xlabel('Frequency (Hz)');
ylabel('Power');
title('Power Spectrum of Motion Energy');
xlim([0, frameRate / 2]); % Limit x-axis to Nyquist frequency

% Perform Continuous Wavelet Transform
[cwtCoeffs, freq] = cwt(motionEnergy, 'amor',frameRate); % Using Morlet wavelet

% Plot the Wavelet Power Spectrum
figure;
surface(1:numFrames, freq, abs(cwtCoeffs));
axis tight;
shading interp;
xlabel('Time (frames)');
ylabel('Frequency (Hz)');
title('Wavelet Power Spectrum of Motion Energy');
colorbar;

% Define window sizes for different time scales (in frames)
windowSizes = [5, 10, 20, 50];

% Calculate standard deviation of motion energy over different window sizes
stdDevOverTime = cell(size(windowSizes));
for i = 1:length(windowSizes)
    windowSize = windowSizes(i);
    stdDevOverTime{i} = movstd(motionEnergy, windowSize);
end

% Plot the results for each window size
figure;
hold on;
for i = 1:length(windowSizes)
    plot(stdDevOverTime{i});
end
legend(arrayfun(@(w) sprintf('%d frames', w), windowSizes, 'UniformOutput', false));
xlabel('Frame');
ylabel('Standard Deviation of Motion Energy');
title('Motion Energy Fluctuations at Different Time Scales');
hold off;

