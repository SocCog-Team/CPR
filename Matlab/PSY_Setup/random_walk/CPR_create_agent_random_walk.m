function AGNT = CPR_create_agent_random_walk(STIM,AGNT)

% Create coherence-dependent confidence
temp = [];
for iCoh = 1:length(STIM.RDP_coherence)
    clear coh_block
    coh_block           = normrnd(STIM.RDP_coherence(iCoh),AGNT.str_sigma, [1,STIM.RDP_coherence_smple(2)]);
    temp              	= [temp coh_block];
end

% Zero-pad, add noise, smooth
pad_vec                 = zeros(1,AGNT.lag);
AGNT.dir_nse            = [pad_vec+randi(359) normrnd(STIM.RDP_direction_deg,AGNT.dir_sigma, [1,length(STIM.RDP_direction_deg)])];
AGNT.dir_smooth         = smoothdata(AGNT.dir_nse, AGNT.smooth_kernel, AGNT.win); %%% problematic! smoothing results in jumps over 0/359deg point

x_smooth                = smoothdata(cos(deg2rad(AGNT.dir_nse)), AGNT.smooth_kernel, AGNT.win);
y_smooth                = smoothdata(sin(deg2rad(AGNT.dir_nse)), AGNT.smooth_kernel, AGNT.win);
AGNT.dir_smooth       	= rad2deg(atan2(y_smooth,x_smooth));
AGNT.dir_smooth       	= mod(AGNT.dir_smooth,360);

AGNT.str_nse            = [pad_vec+rand(1) temp];
AGNT.str_smooth         = smoothdata(AGNT.str_nse, AGNT.smooth_kernel, AGNT.win);
AGNT.str_smooth(AGNT.str_smooth > 1) = 1;
AGNT.str_smooth(AGNT.str_smooth < 0) = 0;

figure;plot(AGNT.str_nse);hold on;plot(AGNT.str_smooth);shg
figure;plot(STIM.RDP_direction_deg);hold on;plot(AGNT.dir_smooth);shg

end
