%% CPR main analysis
% This script processes and analyses CPR electrophysiology data.
% For each recording, it requires
%   .MWK2 Mworks data file
%   .MAT data files that were generated through the ULTRAsort routine
%
% Felix Schneider, CNL
%
% Version history
%   1.0     (fxs 2025-09-05) Initial version.

% TO DO:
% (1) Implement automised loop over sorted recordings
% (2) Add a quality label if channel is worth processing further:
% ---> Add Drift analysis to check is FR remains stable over time

clear all
close all

addpath /Users/cnl/Desktop/CPR/code

cfg_pth         = '/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
source_dir      = '/Users/cnl/Documents/DATA/Nilan/';
import_flag     = false;

rec_lst           = {%'20250807_nil_CPR_block1_phy4_rec045_ann';  % Onset transients before 0?
                    '20250924_nil_CPR_block1_phy4_rec059_fxs',...
                    '20250925_nil_CPR_block1_phy4_rec060_fxs',...
                    '20250926_nil_CPR_block1_phy4_rec061_fxs'};

for iRec = 1%:length(rec_lst)

    rec_info        = split(rec_lst{iRec},'_');
    dest_dir        = ['/Users/cnl/Documents/DATA/Nilan/spike_sorting/' rec_info{1} '_'  rec_info{6} '_'  rec_info{4} '/'];

    %% Import stimulus parameters and neural responses for stimulus cycles
    phy             = PHY_preprocessing(rec_lst{iRec},source_dir,dest_dir,cfg_pth,import_flag);

    %% Do quality checks here...
    phy             = PHY_quality_assessment(phy);

    % Save summary file
    save([dest_dir '/summary_' rec_lst{iRec} '.mat'], 'phy', '-v7.3')

    %% Sort by individual stimulus state
    state           = PHY_sort_spikes_by_state(phy);

    % Save summary file
    save([dest_dir '/state_responses_' rec_lst{iRec} '.mat'], 'state', '-v7.3')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function state = PHY_sort_spikes_by_state(in)
clear state
state_cnt       = 0;
stim            = in.stim;
joy             = in.joy;

for iCyc = 1:length(stim.rdp_dir)
    for iState = 2:length(stim.rdp_dir{iCyc}) % discard last direction from previous cycle
        state_cnt = state_cnt +1;

        % Task
        state.cpr_solo(state_cnt)   = stim.cpr_solo{iCyc};

        % Stimulus
        state.rdp_dir(state_cnt)    = stim.rdp_dir{iCyc}(iState);
        state.rdp_dir_ts(state_cnt) = stim.rdp_dir_ts{iCyc}(iState);
        state.rdp_coh(state_cnt)    = stim.rdp_coh{iCyc}(find(state.rdp_dir_ts(state_cnt) > stim.rdp_coh_ts{iCyc},1,'last'));

        if iState == length(stim.rdp_dir{iCyc})
            state.dur_s(state_cnt)  = double(stim.cpr_cyle(iCyc,2) - stim.rdp_dir_ts{iCyc}(iState)) / 1e6;
        else
            state.dur_s(state_cnt)  = double(stim.rdp_dir_ts{iCyc}(iState+1) - stim.rdp_dir_ts{iCyc}(iState)) /1e6;
        end
        
        %%% Extract target sample ID and outcome here %%%

        % Joystick
        if iState == length(stim.rdp_dir{iCyc})
            js_monk_idx     = joy.js_monk_ts{iCyc} > stim.rdp_dir_ts{iCyc}(iState) & joy.js_monk_ts{iCyc} < double(stim.cpr_cyle(iCyc,2));
            js_hum_idx      = joy.js_hum_ts{iCyc} > stim.rdp_dir_ts{iCyc}(iState) & joy.js_hum_ts{iCyc} < double(stim.cpr_cyle(iCyc,2));
        else
            js_monk_idx     = joy.js_monk_ts{iCyc} > stim.rdp_dir_ts{iCyc}(iState) & joy.js_monk_ts{iCyc} < stim.rdp_dir_ts{iCyc}(iState+1);
            js_hum_idx      = joy.js_hum_ts{iCyc} > stim.rdp_dir_ts{iCyc}(iState) & joy.js_hum_ts{iCyc} < stim.rdp_dir_ts{iCyc}(iState+1);
        end

        state.js_monk_tlt{state_cnt}    = joy.js_monk_tlt{iCyc}(js_monk_idx);
        state.js_monk_dir{state_cnt}    = joy.js_monk_dir{iCyc}(js_monk_idx);
        state.js_hum_dir{state_cnt}     = joy.js_hum_tlt{iCyc}(js_hum_idx);
        state.js_hum_tlt{state_cnt}     = joy.js_hum_dir{iCyc}(js_hum_idx);

        % Brain
        chan = fieldnames(in.brain.CPR.spks);
        chan(strcmp(chan, 'include')) = [];

        for iChan = 1:length(chan)
            unit_label = fieldnames(in.brain.CPR.spks.(chan{iChan}));

            for iUnit = 1:length(unit_label)
                % Skip unsorted spikes
                if strcmp(unit_label{iUnit}, 'unit0')
                    continue
                end

                % Unit included based on onset response?
                tmp_id = [chan{iChan} '_' unit_label{iUnit}];
                unit_idx = in.brain.CPR.spks.include.unit_ID == tmp_id;

                if ~in.brain.CPR.spks.include.inclusion_flag(unit_idx)
                    continue
                end

                clear dat
                dat = in.brain.CPR.spks.(chan{iChan}).(unit_label{iUnit});

                if iState == length(stim.rdp_dir{iCyc})
                    spk_idx = dat{iCyc} > stim.rdp_dir_ts{iCyc}(iState)-double(stim.cpr_cyle(iCyc,1)) & dat{iCyc} < double(stim.cpr_cyle(iCyc,2))-double(stim.cpr_cyle(iCyc,1));
                else
                    spk_idx = dat{iCyc} > stim.rdp_dir_ts{iCyc}(iState)-double(stim.cpr_cyle(iCyc,1)) & dat{iCyc} < stim.rdp_dir_ts{iCyc}(iState+1)-double(stim.cpr_cyle(iCyc,1));
                end

                unit_id                             = [chan{iChan} '_' unit_label{iUnit}];
                state.spk_ts.(unit_id){state_cnt}   = dat{iCyc}(spk_idx) - (stim.rdp_dir_ts{iCyc}(iState) - double(stim.cpr_cyle(iCyc,1)));
                state.spk_n.(unit_id)(state_cnt)    = length(dat{iCyc}(spk_idx));
                state.cIdx(state_cnt) = iCyc;
                % add baseline

                if in.brain.CPR.spks.include.cyc_id{unit_idx}(iCyc) == 0
                    state.include.(unit_id)(state_cnt)  = false;
                else
                    state.include.(unit_id)(state_cnt)  = true;
                end
            end
        end
    end
end
end

function gauss = fit_gaussian(spks, time)

sigma               = .005;                     % Width of gaussian/window [s]

% For every spike
for iSpk = 1:length(spks)

    % Center gaussian at spike time
    mu              = spks(iSpk);

    % Calculate gaussian
    p1              = -.5 * ((time - mu)/sigma) .^ 2;
    p2              = (sigma * sqrt(2*pi));
    gauss(iSpk,:)   = exp(p1) ./ p2;
end
end

function [all, sdf_g, sdf_a] = FR_estimation(spike_times, time, plot_flag)

sdf                 = [];
all                 = [];

for iTrial = 1:length(spike_times)

    spks            = double(spike_times{iTrial}') ./1e6;  	% Get all spikes of respective trial

    all             = [all spks];                   % Concatenate spikes of all trials
    xspikes         = repmat(spks,3,1);             % Replicate array
    yspikes      	= nan(size(xspikes));           % NaN array

    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                   	% Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end

    % Plot trial raster
    if plot_flag
        pl           = plot(xspikes, yspikes, 'Color', 'k', 'LineWidth',1.25);
    end
    % Spike density function
    if isempty(spks)
    else

        % Fit function to spikes
        gauss_kern      = fit_gaussian(spks, time); % problem --> backward smearing
        alpha_kern      = fit_alpha(spks, time); % problem: values too small - fix bug!

        % Sum over all distributions to get spike density function
        sdf_g(iTrial,:)	= sum(gauss_kern,1);
        sdf_a(iTrial,:)	= sum(alpha_kern,1);
    end

end
end

function alpha_kern = fit_alpha(spks, time)
    tau = 0.005; % [s]

    nSpks = length(spks);
    nTime = length(time);
    alpha_kern = zeros(nSpks, nTime);

    for iSpk = 1:nSpks
        t_rel = time - spks(iSpk); % time relative to spike

        % Only consider t >= 0 (causal)
        k = zeros(1, nTime);
        idx = t_rel >= 0;

        k(idx) = (t_rel(idx) ./ tau) .* exp(1 - t_rel(idx) ./ tau);

        % Normalize to unit area (like Gaussian)
        k = k / (exp(1) * tau);

        alpha_kern(iSpk, :) = k;
    end
end

function out = test_onset_response(spks_cyc)

% Crop spike trains for efficiency
for iCyc = 1:length(spks_cyc)
    spks_cyc{iCyc}(spks_cyc{iCyc}>1e6) = [];
end

% Calculate spike density function
tstep           = .001;
time            = [-.3:tstep:1];
% figure;hold on
[~, sdf_g, sdf_a]   = FR_estimation(spks_cyc, time, false);
msdf                 = mean(sdf_a); % Average time course

if isempty(sdf_a)
    out.mean_fr = nan;
    out.lat     = nan;
    out.p       = nan;
    out.h       = nan;
    out.stats.zval = nan;
    return
end

% Define time windos of interest
BLidx           = 1:200;
STIMidx         = 301:500;
out.mean_fr     = mean(msdf(STIMidx));

% figure
% plot(msdf)

% Test average activity in time windows
mBL             = mean(sdf_a(:,BLidx),2);  % Avg response to baseline
mStimOn         = mean(sdf_a(:,STIMidx),2); % Avg response to stimulus
[out.p,out.h,out.stats] = signrank(mStimOn, mBL);

if ~isfield(out.stats, 'zval')
    out.stats.zval = nan;
end

% figure;
% boxplot([mBL;mStimOn],[zeros(length(mBL),1); ones(length(mStimOn),1)], 'Notch', 'on')
% ylabel('Firing rate [norm]'); set(gca, 'xticklabels', {'BL', 'STIM'});

%%% Check SNR threshold and latency %%%
nStd            = 3;        % No. standard deviations
nBins           = 9;        % No. bins

% Define significance threshold based on baseline response
thresh    	= [mean(msdf(BLidx)) + (std(msdf(BLidx)) * nStd), ...
    mean(msdf(BLidx)) - (std(msdf(BLidx)) * nStd)];

% For response in defined time window...
exc = [];
for s = STIMidx(1:end-nBins)
    % Check at each timepoint if all bins in 10ms window exceed threshhold
    if sum(msdf(s:s+nBins) >= thresh(1)) == (nBins+1) || sum(msdf(s:s+nBins) <= thresh(2)) == (nBins+1)
        exc(s) = s; % Save sample location
    end
end

if ~isempty(exc)
    out.lat = find(exc, 1,'first') - STIMidx(1); % Extract latency of first significant bin
else
    out.lat = nan;
end
end

function  phy = PHY_quality_assessment(phy)
% Define table fields
var     = {'unit_ID','string';          % unit ID
    'n_Cycles', 'double'; ...           % Number of usable CPR cycles
    'mean_FR', 'double'; ...            % Average firing rate
    'signrank_p', 'double'; ...         % signrank p-value
    'signrank_z', 'double'; ...         % signrank z-value
    'thresh_latency_ms', 'double';...;  % latency of threshold crossing
    'inclusion_flag', 'logical'};       % inclusion flag

% Initialse table
t       = table('Size',[1000, size(var,1)],...
    'VariableTypes',var(:,2),...
    'VariableNames',var(:,1));

% Reset counter
cc = 0;

% Define inclusion criteria
min_num_cycles          = 30;   % Min Number of stimulus cycles
min_cycle_dur_sec       = 5;    % Min cycle duraction [seconds]
min_fr                  = 1;    % Min firing rate
p_thresh                = .05;  % Significance threshold
min_lat_ms              = 35;   % Min onset response latency
max_lat_ms              = 200;  % Max onset response latency

chan_lst                                = fieldnames(phy.brain.CPR.spks);
chan_lst(strcmp(chan_lst, 'include'))   = [];

for iChan = 1:length(chan_lst)
    unit_lst            = fieldnames(phy.brain.CPR.spks.(chan_lst{iChan}));

    for iUnit = 1:length(unit_lst)
        clear spk_times_cycle spk_times_cycle_filt test_struct cyc_duration_sec avg_fr_cycle
        spk_times_cycle         = phy.brain.CPR.spks.(chan_lst{iChan}).(unit_lst{iUnit});

        % Remove cycles with a low firing rate and shorter than 1 second.
        cyc_duration_sec        = (phy.stim.cpr_cyle(:,2) - phy.stim.cpr_cyle(:,1)) ./1e6;
        avg_fr_cycle            = cellfun(@length, spk_times_cycle)' ./ cyc_duration_sec;
        cyc_lst                 = avg_fr_cycle > min_fr & cyc_duration_sec > min_cycle_dur_sec;
        spk_times_cycle_filt    = spk_times_cycle(avg_fr_cycle > min_fr & cyc_duration_sec > min_cycle_dur_sec);

        % Test onset response of units with sufficient repetitions
        cc                      = cc+1;
        t.unit_ID(cc)           = [(chan_lst{iChan}) '_' (unit_lst{iUnit})];

        if length(spk_times_cycle_filt) >= min_num_cycles
            % Test onset response against baseline activity
            test_struct             = test_onset_response(spk_times_cycle_filt);

            % Save results to table
            t.n_Cycles(cc)          = length(spk_times_cycle_filt); % Number of remaining stimulus cycles for analysis
            t.cyc_id{cc}            = cyc_lst; % Cycle numbers
            t.mean_FR(cc)           = test_struct.mean_fr; % Avg. firing rate
            t.signrank_p(cc)        = test_struct.p; % P-value baseline vs onset response
            t.signrank_z(cc)        = test_struct.stats.zval; % Z-value baseline vs onset response
            t.thresh_latency_ms(cc) = test_struct.lat; % Latency of threshold crossing

            % Include if stable, above threshold response with expected latency
            if (test_struct.lat > min_lat_ms && test_struct.lat < max_lat_ms) && test_struct.p < p_thresh && test_struct.mean_fr > min_fr && length(spk_times_cycle_filt) > min_num_cycles
                t.inclusion_flag(cc) = true;
            else
                t.inclusion_flag(cc) = false;
            end
        else
            t.n_Cycles(cc)          = nan;
            t.cyc_id{cc}            = nan;
            t.mean_FR(cc)           = nan;
            t.signrank_p(cc)        = nan;
            t.signrank_z(cc)        = nan;
            t.thresh_latency_ms(cc) = nan;
            t.inclusion_flag(cc)    = false;
        end
    end
end

% Crop table
t(ismissing(t.unit_ID),:)           = [];
phy.brain.CPR.spks.include          = t;
end