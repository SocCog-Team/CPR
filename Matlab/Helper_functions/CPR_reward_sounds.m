%% CPR reward sound generation

% Params
fs              = 44100; % sampling rate
snd_dur         = .2; % sound duration [s]
ramp_dur        = .05; % on/off ramp duration [s]
t               = 0:1/fs:snd_dur-1/fs; % time vector
frq             = [261.63 293.66 329.63 349.23 392.00 440.00 493.88 523.25]*2; % frequencies
wav_name        = {'1','2','3','4','5','6','7','8',}; % file names
dest_dir        = '/Users/fschneider/Desktop/sound/'; % destination directory

% Reward sounds: pure tones
for iFrq = 1:length(frq)
    pt(iFrq,:)	= sin(2*pi*frq(iFrq)*t);
end

% Error sound: noise burst 
nsamples        = floor(snd_dur*fs);
nse             = bandpass(randn(nsamples,1),[500 8000],fs); % bandpass-filtered noise
nse             = nse';
nse             = nse/max(nse); % normalisation

% Create on/off ramps
rampSamps       = floor(fs*ramp_dur);
window          = hanning(2*rampSamps)';
w1              = window(1:ceil((length(window))/2));
w2              = window(ceil((length(window))/2)+1:end);
w_on            = [w1 ones(1,length(nse)-length(w1))];
w_off           = [ones(1,length(nse)-length(w2)) w2];

% Ramp all stimuli
for iFrq = 1:length(frq)
    pt_r(iFrq,:)= pt(iFrq,:).*w_on.*w_off;
%     plot(pt_r(iFrq,:))
%     sound(pt_r(iFrq,:),fs)
%     pause(1)
end

nse_r           = nse.*w_on.*w_off;

% Write files - generate sound for left & right speaker
for iFrq = 1:length(frq)
    audiowrite([dest_dir wav_name{iFrq} 'R.wav'],[zeros(nsamples,1) pt_r(iFrq,:)'],fs)
    audiowrite([dest_dir wav_name{iFrq} 'L.wav'],[pt_r(iFrq,:)' zeros(nsamples,1)],fs)
end

audiowrite([dest_dir 'noise_R.wav'],[zeros(nsamples,1) nse_r'],fs)
audiowrite([dest_dir 'noise_L.wav'],[nse_r' zeros(nsamples,1)],fs)
