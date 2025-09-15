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
% (2) Save data per recording and channel in appropriate format with ID
% (3) Find efficient way to import and process saved data for further
% analysis between sessions and conditions
% (4) Add a quality label if channel is worth processing further:
% ---> Add Drift analysis to check is FR remains stable over time
% ---> Check if trial numbers are sufficient after cleanup
% (5) Outsource plotting routines

clear all
close all
addpath /Users/cnl/Desktop/CPR/code

cfg_pth         = '/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
source_dir      = '/Users/cnl/Documents/DATA/Nilan/';
% dest_dir        = '/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250807_rec045/'; % Onset transients before 0?
% fname           = '20250807_nil_CPR_block1_phy4_rec045_ann';
dest_dir       = '/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250903_rec050/';
fname           = '20250903_nil_CPR_block1_phy4_rec050_ann';
import_flag     = false;

% Import stimulus parameters and neural responses for stimulus cycles
phy             = PHY_preprocessing(fname,source_dir,dest_dir,cfg_pth,import_flag);

%% Do quality checks here...
phy             = PHY_quality_assessment(phy);

% Save summary file
save([dest_dir '/summary_' fname '.mat'], 'phy', '-v7.3')

%% Sort by individual stimulus state
state           = PHY_sort_spikes_by_state(phy);

% Save summary file
save([dest_dir '/state_responses_' fname '.mat'], 'state', '-v7.3')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = PHY_sort_spikes_by_state(in)
clear state
state_cnt = 0;
cyc = in.cyc;

for iCyc = 1:length(in.cyc.rdp_dir)
    for iState = 2:length(cyc.rdp_dir{iCyc}) % discard last direction from previous cycle
        state_cnt                   = state_cnt+1;

        % Stimulus
        state.rdp_dir(state_cnt)    = cyc.rdp_dir{iCyc}(iState);
        state.rdp_dir_ts(state_cnt) = cyc.rdp_dir_ts{iCyc}(iState);
        state.rdp_coh(state_cnt)    = cyc.rdp_coh{iCyc}(find(state.rdp_dir_ts(state_cnt) > cyc.rdp_coh_ts{iCyc},1,'last'));

        if iState == length(in.cyc.rdp_dir{iCyc})
            state.dur_s(state_cnt)  = double(cyc.cEnd(iCyc) - cyc.rdp_dir_ts{iCyc}(iState)) / 1e6;
        else
            state.dur_s(state_cnt)  = double(cyc.rdp_dir_ts{iCyc}(iState+1) - cyc.rdp_dir_ts{iCyc}(iState)) /1e6;
        end

        % Joystick
        if iState == length(in.cyc.rdp_dir{iCyc})
            js_idx = cyc.js_ts{iCyc} > cyc.rdp_dir_ts{iCyc}(iState) & cyc.js_ts{iCyc} < double(cyc.cpr_cyle(iCyc,2));
        else
            js_idx = cyc.js_ts{iCyc} > cyc.rdp_dir_ts{iCyc}(iState) & cyc.js_ts{iCyc} < cyc.rdp_dir_ts{iCyc}(iState+1);
        end
        state.js_tlt{state_cnt}     = cyc.js_tlt{iCyc}(js_idx);
        state.js_dir{state_cnt}     = cyc.js_dir{iCyc}(js_idx);

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

                clear dat
                dat = in.brain.CPR.spks.(chan{iChan}).(unit_label{iUnit});

                if iState == length(in.cyc.rdp_dir{iCyc})
                    spk_idx = dat{iCyc} > cyc.rdp_dir_ts{iCyc}(iState) & dat{iCyc} < double(cyc.cpr_cyle(iCyc,2));
                else
                    spk_idx = dat{iCyc} > cyc.rdp_dir_ts{iCyc}(iState) & dat{iCyc} < cyc.rdp_dir_ts{iCyc}(iState+1);
                end

                unit_id = [chan{iChan} '_' unit_label{iUnit}];
                state.spk_ts.(unit_id){state_cnt} = dat{iCyc}(spk_idx) - cyc.rdp_dir_ts{iCyc}(iState);
                state.spk_n.(unit_id)(state_cnt) = length(dat{iCyc}(spk_idx));
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

function [all, sdf] = FR_estimation(spike_times, time, plot_flag)

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

        % Fit gaussian to spikes
        gauss           = fit_gaussian(spks, time);

        % Sum over all distributions to get spike density function
        sdf(iTrial,:)	= sum(gauss,1);
    end

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
[~, sdf]        = FR_estimation(spks_cyc, time, false);

% Defien time windos of interest
BLidx           = 1:200;
STIMidx         = 301:500;

% Test average activity in time windows
mBL             = mean(sdf(:,BLidx),2);  % Avg response to baseline
mStimOn         = mean(sdf(:,STIMidx),2); % Avg response to stimulus
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
msdf            = mean(sdf);% Average time course

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

chan_lst = fieldnames(phy.brain.CPR.spks);
chan_lst(strcmp(chan_lst, 'include')) = [];

for iChan = 1:length(chan_lst)
    unit_lst = fieldnames(phy.brain.CPR.spks.(chan_lst{iChan}));
    for iUnit = 1:length(unit_lst)
        spk_times_cycle = phy.brain.CPR.spks.(chan_lst{iChan}).(unit_lst{iUnit});
        
        % Remove cells that are empty or cycles shorter than 1 second.
        % We expect at least 1 spike per second
        % That way, we also exclude empty trials due to spike sorting
        cyc_duration_sec = (phy.cyc.cpr_cyle(:,2) - phy.cyc.cpr_cyle(:,1)) ./1e6;
        spk_times_cycle_filt = spk_times_cycle(~cellfun('isempty', spk_times_cycle)' & cyc_duration_sec > 1);

        % Test onset response against baseline activity
        cc                      = cc+1;
        test_struct             = test_onset_response(spk_times_cycle_filt);
        
        % Save results to table
        t.unit_ID(cc)           = [(chan_lst{iChan}) '_' (unit_lst{iUnit})];
        t.n_Cycles(cc)          =  length(spk_times_cycle_filt);
        t.signrank_p(cc)        = test_struct.p;
        t.signrank_z(cc)        = test_struct.stats.zval;
        t.thresh_latency_ms(cc) = test_struct.lat;

        % Include if reasonable, above threshold response  
        if (test_struct.lat > 35 && test_struct.lat < 150) && test_struct.p < .05 && length(spk_times_cycle_filt) > 50
            t.inclusion_flag(cc) = true;
        else
            t.inclusion_flag(cc) = false;
        end

    end
end

% Crop table
t(ismissing(t.unit_ID),:)       = [];
phy.brain.CPR.spks.include      = t;
end