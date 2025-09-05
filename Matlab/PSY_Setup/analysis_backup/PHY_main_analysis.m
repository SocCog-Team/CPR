
cfg_pth         = '/Users/cnl/Desktop/CPR/code/felix_nhp_solo.cfg';
source_dir      = '/Users/cnl/Documents/DATA/Nilan/';
dest_dir        = '/Users/cnl/Documents/DATA/Nilan/spike_sorting/20250807_rec045/';
fname           = '20250807_nil_CPR_block1_phy4_rec045_ann';
import_flag     = false;

phy             = PHY_preprocessing(fname,source_dir,dest_dir,cfg_pth,import_flag);
state           = sort_spikes_by_state(phy);

%% Plot

plot_RF(phy.brain.RF.stim_id, phy.brain.RF.stim_pos, phy.brain.RF.ch001_neg.nSpikes)
plot_CPRtuning(state.spk_n.ch001_neg_3, deg2rad(state.rdp_dir))

%% Plot CPR raster

close all
tstep = .001;
time =  [0:tstep:1];

% Correct by cycle onset
for iCyc = 1:length(phy.cyc.cpr_cyle)
    spk_times{iCyc} = phy.brain.CPR.spks.cyc.ch001_neg.unit1{iCyc} - phy.cyc.cpr_cyle(iCyc,1);
end

% [dir_val, index] = sort(state.rdp_dir)

figure;hold on
[all, sdf] = FR_estimation(spk_times, time, true);

figure;hold on
plot(mean(sdf))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function state = sort_spikes_by_state(in)
clear state
state_cnt = 0;
cyc = in.cyc;

for iCyc = 1:length(in.cyc.rdp_dir)
    for iState = 2:length(cyc.rdp_dir{iCyc}) % discard last direction
        state_cnt                   = state_cnt+1;
        state.rdp_dir(state_cnt)    = cyc.rdp_dir{iCyc}(iState);

        if iState == length(in.cyc.rdp_dir{iCyc})
            state.dur_ms(state_cnt) = cyc.cEnd(iCyc) - cyc.rdp_dir_ts{iCyc}(iState) / 1e3;
        else
            state.dur_ms(state_cnt) = cyc.rdp_dir_ts{iCyc}(iState+1) - cyc.rdp_dir_ts{iCyc}(iState) / 1e3;
        end

        chan = fieldnames(in.brain.CPR.spks.cyc);

        for iChan = 1:length(chan)
            for iUnit = 1:length(fieldnames(in.brain.CPR.spks.cyc.(chan{iChan})))

                clear dat
                dat = in.brain.CPR.spks.cyc.(chan{iChan}).(['unit' num2str(iUnit)]);

                if iState == length(in.cyc.rdp_dir{iCyc})
                    spk_idx = dat{iCyc} > cyc.rdp_dir_ts{iCyc}(iState) & dat{iCyc} < double(cyc.cpr_cyle(iCyc,2));
                else
                    spk_idx = dat{iCyc} > cyc.rdp_dir_ts{iCyc}(iState) & dat{iCyc} < cyc.rdp_dir_ts{iCyc}(iState+1);
                end

                unit_id = [chan{iChan} '_' num2str(iUnit)];
                state.spk_ts.(unit_id){state_cnt} = dat{iCyc}(spk_idx);
                state.spk_n.(unit_id)(state_cnt) = length(dat{iCyc}(spk_idx));
            end
        end
    end
end
end

function plot_RF(stim_id, stim_pos, nSpikes)

for iUnit = 1:size(nSpikes,1)
    f                   	= figure;
    mat_sum                 = nan(size(stim_id));

    for iPos = 1:max(max(stim_id))
        mat_sum(stim_id == iPos) = sum(nSpikes(iUnit, stim_pos == iPos));
    end

    % hm                      = imagesc(mat_sum);
    hm                      = imagesc(imgaussfilt(mat_sum,.75));
    
    x_position              = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
    y_position              = fliplr([-12 -9 -6 -3 0 3 6 9 12]);

    ax.XTick                = [1 5 9];
    ax.XTickLabel           = {num2str(x_position(1)) num2str(x_position(5)) num2str(x_position(9))};
    ax.YTick                = [1 5 9];
    ax.YTickLabel           = {num2str(y_position(1)) num2str(y_position(5)) num2str(y_position(9))};
    colormap(jet(256)); colorbar
end
end

function plot_CPRtuning(r, theta) 

figure
h = polaraxes;                    
polarscatter(h,theta,r)

h.ThetaZeroLocation = 'top';         
h.ThetaDir = 'clockwise';            
hold on

nbins = 360/30;                       
edges = linspace(0, 2*pi, nbins+1); 
r_mean_bin = zeros(1, nbins);
theta_bin = zeros(1, nbins);

for i = 1:nbins
    idx = theta >= edges(i) & theta < edges(i+1); 
    r_mean_bin(i) = mean(r(idx));             
    theta_bin(i) = (edges(i)+edges(i+1))/2;       
end
polarplot(h, theta_bin, r_mean_bin, 'r-.', 'LineWidth', 2)

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
