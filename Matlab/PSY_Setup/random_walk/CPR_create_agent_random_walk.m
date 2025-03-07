function AGNT = CPR_create_agent_random_walk(STIM,AGNT)

% Scale by coherence
tmp = [];
nsamples_coh_block = unique(diff(STIM.RDP_coherence_smple));
for iCoh = 1:length(STIM.RDP_coherence_smple)
    tmp = [tmp repmat(STIM.RDP_coherence(iCoh),1,nsamples_coh_block)];
end

% Zero-pad, add noise, smooth
pad_vec                 = zeros(1,AGNT.lag);
AGNT.dir_nse            = [pad_vec+randi(359) normrnd(STIM.RDP_direction_deg,AGNT.dir_sigma, [1,length(STIM.RDP_direction_deg)])];
AGNT.str_nse            = [pad_vec+rand(1) normrnd(tmp, AGNT.str_sigma, [1,length(STIM.RDP_direction_deg)])];
AGNT.dir_smooth         = smoothdata(AGNT.dir_nse, AGNT.smooth_kernel, AGNT.win);
AGNT.str_smooth         = smoothdata(AGNT.str_nse, AGNT.smooth_kernel, AGNT.win);

end
