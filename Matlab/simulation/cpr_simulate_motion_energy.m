if 0
    
% Parameters
numFrames = 200;            % Total number of frames
numDots = 200;              % Number of dots
coherence = 0.5;              % Coherence level (e.g., 50%)
coherentDirection = pi/4;   % Nominal direction of coherent motion (e.g., 45 degrees)
frameRate = 100;            % Frame rate in Hz (100 Hz)

% Pre-allocate arrays for motion energy and directional deviations
directionDeviation = zeros(1, numFrames);
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
end

% Parameters
numFrames = 300;           % Total number of frames
numDots = 100;              % Number of dots in the RDK
coherence = 0.5;            % Coherence level (0 to 1 for 0-100%)
coherentDirection = pi/4;   % Direction of coherent motion (e.g., 45 degrees in radians)
frameRate = 100;            % Frame rate in Hz (e.g., 100 Hz)
dotLifetime = 30;           % Dot lifetime in frames (e.g., 30 frames)
makeMovie = 0;

% Initialize dot positions and lifetimes
dotX = rand(numDots, 1) * 2 - 1;          % X positions (normalized between -1 and 1)
dotY = rand(numDots, 1) * 2 - 1;          % Y positions (normalized between -1 and 1)
dotLifetimes = randi([1 dotLifetime], numDots, 1);  % Random initial lifetime for each dot

% Initialize arrays for motion energy (projected) and overall motion energy
motionEnergy = zeros(1, numFrames);       % Projected motion energy onto coherent direction
overallMotionEnergy = zeros(1, numFrames); % Overall motion energy
directionDeviation = zeros(1, numFrames);


if makeMovie
    % Set up video writer
    videoFile = 'rdk_motion_video.avi';
    v = VideoWriter(videoFile, 'Motion JPEG AVI'); % Choose codec; can also use 'MPEG-4' for .mp4
    v.FrameRate = frameRate;   % Set video frame rate
    open(v);
end

for frame = 1:numFrames
    % 1. Update dot lifetimes and reset dots that reach their lifetime
    dotLifetimes = dotLifetimes - 1;          % Decrement lifetime for each dot
    respawnDots = dotLifetimes <= 0;          % Find dots whose lifetime has ended
    dotX(respawnDots) = rand(sum(respawnDots), 1) * 2 - 1; % Reset x position
    dotY(respawnDots) = rand(sum(respawnDots), 1) * 2 - 1; % Reset y position
    dotLifetimes(respawnDots) = dotLifetime;  % Reset lifetime

    % 2. Generate random directions for all dots
    randomDirections = 2 * pi * rand(numDots, 1);  % Random directions in radians

    % 3. Assign coherence: set a fraction of dots to move in the coherent direction
    coherentDots = round(coherence * numDots);     % Number of coherent dots
    randomDots = numDots - coherentDots;           % Number of random dots

    % Set coherent direction for coherent dots
    randomDirections(1:coherentDots) = coherentDirection;
    
    % 4. Calculate x and y components of velocity for each dot
    xVelocities = cos(randomDirections);   % x component of each dot's velocity
    yVelocities = sin(randomDirections);   % y component of each dot's velocity

    % 5. Update positions of all dots
    dotX = dotX + xVelocities * 0.01;  % Update x position (0.01 is an arbitrary speed factor)
    dotY = dotY + yVelocities * 0.01;  % Update y position (same speed factor)
    
    % Wrap-around: Dots that move out of bounds reappear on the opposite side
    dotX(dotX > 1) = dotX(dotX > 1) - 2;
    dotX(dotX < -1) = dotX(dotX < -1) + 2;
    dotY(dotY > 1) = dotY(dotY > 1) - 2;
    dotY(dotY < -1) = dotY(dotY < -1) + 2;

    % 6. Sum x and y velocities across all dots to get net x and y movement
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
    
    
    % 7. Project net movement onto the coherent direction to get projected motion energy
    % This is the dot product of (netX, netY) with the coherent direction vector
    motionEnergy(frame) = netX * cos(coherentDirection) + netY * sin(coherentDirection);
    
    % 8. Calculate overall motion energy as the magnitude of the net movement vector
    % (without projecting onto any specific direction)
    overallMotionEnergy(frame) = sqrt(netX^2 + netY^2);
    
    if makeMovie
        % 9. Plot the dots for the current frame
        figure(1); clf;
        plot(dotX, dotY, 'ko', 'MarkerSize', 4, 'MarkerFaceColor', 'k');
        xlim([-1, 1]); ylim([-1, 1]);
        axis square;
        title(['Frame ' num2str(frame)]);
        drawnow;

        % Capture the plot as a frame in the video
        frameData = getframe(gcf);
        writeVideo(v, frameData);
    end

end

if makeMovie
    % Close the video writer
    close(v);
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
windowLengths = [5, 10, 20, 50];

% Calculate standard deviation of motion energy over different window sizes
stdDevOverTime = cell(size(windowLengths));
for i = 1:length(windowLengths)
    windowSize = windowLengths(i);
    stdDevOverTime{i} = movstd(motionEnergy, windowSize);
end

% Plot the results for each window size
figure;
hold on;
for i = 1:length(windowLengths)
    plot(stdDevOverTime{i});
end
legend(arrayfun(@(w) sprintf('%d frames', w), windowLengths, 'UniformOutput', false));
xlabel('Frame');
ylabel('Standard Deviation of Motion Energy');
title('Motion Energy Fluctuations at Different Time Scales');
hold off;


% Define parameters for smoothing and shifting
windowLengths = [5, 10, 20, 50]; % Running average window lengths in frames
shiftFrames = [20 40];      % Number of frames to shift the smoothed signals

% Remove mean
motionEnergy = motionEnergy - mean(motionEnergy);

% Initialize matrix to store smoothed signals
smoothedSignals = zeros(length(windowLengths), numFrames);

% Apply running average filter with different window lengths
for i = 1:length(windowLengths)
    window = windowLengths(i);
    % Use 'movmean' function to smooth the motion energy
    smoothedSignals(i, :) = smoothdata(motionEnergy, 'movmean', window);
end

% Plot original and smoothed signals for comparison
time = (1:numFrames) / frameRate; % Time in seconds

figure;
subplot(length(windowLengths) + 1, 1, 1);
plot(time, motionEnergy, 'k');
title('Original Motion Energy');
xlabel('Time (s)');
ylabel('Energy');

for i = 1:length(windowLengths)
    subplot(length(windowLengths) + 1, 1, i + 1);
    plot(time, smoothedSignals(i, :));
    title(['Smoothed with Window Length = ' num2str(windowLengths(i)) ' frames']);
    xlabel('Time (s)');
    ylabel('Energy');
end



% Calculate cross-correlation between original and shifted smoothed signals
maxLag = 75; % Max lag for cross-correlation in frames
lags = -maxLag:maxLag; % Range of lags to compute
scaleopt = 'none';

figure;
for i = 1:length(windowLengths)
    % Original smoothed signal cross-correlation
    [xcorrValuesOriginal, lags] = xcorr(motionEnergy, smoothedSignals(i, :), maxLag, scaleopt);
    
    subplot(length(windowLengths), 1, i);
    plot(lags / frameRate, xcorrValuesOriginal, 'DisplayName', 'No Shift'); % Convert lag to seconds for x-axis
    hold on;
    
    % Calculate cross-correlation for shifted versions
    for j = 1:length(shiftFrames)
        % Shift the smoothed signal by specified frames, filling shifted regions with NaN
        shiftedSignal = circshift(smoothedSignals(i, :), shiftFrames(j));
        
        % Calculate cross-correlation with the original motion energy
        [xcorrValuesShifted, ~] = xcorr(motionEnergy, shiftedSignal, maxLag, scaleopt);
        
        % Plot cross-correlation for each shift
        plot(lags / frameRate, xcorrValuesShifted, 'DisplayName', ['Shift ' num2str(shiftFrames(j)) ' frames']);
    end
    
    title(['Cross-Correlation (Window Length = ' num2str(windowLengths(i)) ' frames)']);
    xlabel('Lag (s)');
    ylabel('Correlation Coefficient');
    legend;
    hold off;
end


