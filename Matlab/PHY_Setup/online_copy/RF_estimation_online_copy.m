function [out, nRep, stim_pos] = RF_estimation_online(d,f,ff,ax_lst,in,nRep,stim_pos)

if nargin < 7
    stim_pos                = [];
end

if nargin < 6
    nRep                    = 0;
end

if nargin < 5
    in                      = [];
end

%% Process

idx.tStart                  = d.event == 'TRIAL_start';
idx.tEnd                    = d.event == 'TRIAL_end';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.type                    = d.event == 'TRIAL_type';
idx.spikes                  = d.event == 'IO_spikes';

trl.start                   = d.time(idx.tStart);
trl.end                     = d.time(idx.tEnd);

channel                     = 1:64;

x_position                  = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_position                  = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
win_offset                  = [50e3 50e3];
stim_id                     = reshape(1:(length(x_position)*length(y_position)),[length(y_position),length(x_position)]);

% Trial index
trlIdx                      = d.time >= trl.start(end) & d.time <= trl.end(end);

% Trial outcome
outcome                     = getTrialData(d.value, trlIdx, idx.outcome);
ttype                       = getTrialData(d.value, trlIdx, idx.type);
rdp_time                    = getTrialData(d.time, trlIdx, idx.type);

% Abort if fixation break
if strcmp(outcome, 'fixation break')
    out = in;
    return
end

% Baseline spiking
blIdx                   = [];
blIdx                   = d.time >= trl.start & d.time <= rdp_time(1);
bl_spike_time           = getTrialData(d.time, blIdx, idx.spikes);
bl_spike_id          	= getTrialData(d.value, blIdx, idx.spikes);
bl_spike_id             = floor(bl_spike_id); % channel-wise for now

% Stimulus position loop
nSpikes                 = [];
cnt                     = 0;

for iStim = 2:length(ttype)
    % Stimulus index - offset by certain lag
    sIdx                    = [];
    sIdx                    = d.time >= rdp_time(iStim-1)+win_offset(1) & d.time <= rdp_time(iStim)+win_offset(2);
    cnt                     = cnt+1;

    % Stimulus position
    pos                     = split(ttype{iStim-1},'_');
    rdp_x                   = str2num(pos{1}(2:end));
    rdp_y                   = str2num(pos{2}(2:end));
    nRep                    = nRep+1;
    stim_pos(nRep)          = stim_id(rdp_y == y_position, rdp_x == x_position);

    % Spike count
    spike_time              = getTrialData(d.time, sIdx, idx.spikes);
    spike_id                = getTrialData(d.value, sIdx, idx.spikes);
    spike_id                = floor(spike_id); % channel-wise for now

    % Unit loop
    for iUnit = 1:length(channel)
        bl_spikes          	= length(bl_spike_time(bl_spike_id == channel(iUnit)));
        nSpikes(iUnit,cnt) = length(spike_time(spike_id == channel(iUnit))) - bl_spikes;
    end
end

out = [in nSpikes];

%% Plot %%%
spk_sum = nan(size(stim_id));

% Every N-th trial...
if mod(cell2mat(d.value(idx.tStart)),10) == 0
    %     close all
    %     f                               = figure('units','normalized','position',[0 0 1 1]);
    figure(f);

    for iUnit = 1:32

        %             % Continue if no spikes detected
        %             if sum(sum(nansum(out{iUnit},3))) == 0
        %                 continue
        %             end

        % Assign axis
        axes(ax_lst(iUnit));

        % Plot data
        for iPos = 1:stim_id(end,end)
            spk_sum(logical(stim_id == iPos)) = sum(out(iUnit, stim_pos == iPos));
        end

        im                          = imagesc(spk_sum);         % Maybe use: contour(X,Y,Z,[1,1])

        ax                          = gca;
        ax.FontSize                 = 16;
        ax.XLabel.String            = ['[Channel ' num2str(iUnit) '] dva'];
        ax.YLabel.String            = 'dva';
        ax.XTick                    = [2:2:8];
        ax.XTick                    = [1 5 9];
        ax.XTickLabel               = {num2str(x_position(1)) num2str(x_position(5)) num2str(x_position(9))};
        ax.YTick                    = [1 5 9];
        ax.YTickLabel               = {num2str(y_position(1)) num2str(y_position(5)) num2str(y_position(9))};

        cb                          = colorbar('northoutside');
        cb.Label.String           	= 'BL-corr. spike count';
        cb.FontSize                 = 10;
    end
    
    cmap                        = [0 0 0;jet(256)]; colormap(cmap); shg;

    figure(ff);

    for iUnit = 33:64

        %             % Continue if no spikes detected
        %             if sum(sum(nansum(out{iUnit},3))) == 0
        %                 continue
        %             end

        % Assign axis
        ax_lst
        axes(ax_lst(iUnit));

        % Plot data
        for iPos = 1:stim_id(end,end)
            spk_sum(logical(stim_id == iPos)) = sum(out(iUnit, stim_pos == iPos));
        end

        im                          = imagesc(spk_sum);         % Maybe use: contour(X,Y,Z,[1,1])

        ax                          = gca;
        ax.FontSize                 = 16;
        ax.XLabel.String            = ['[Channel ' num2str(iUnit) '] dva'];
        ax.YLabel.String            = 'dva';
        ax.XTick                    = [2:2:8];
        ax.XTick                    = [1 5 9];
        ax.XTickLabel               = {num2str(x_position(1)) num2str(x_position(5)) num2str(x_position(9))};
        ax.YTick                    = [1 5 9];
        ax.YTickLabel               = {num2str(y_position(1)) num2str(y_position(5)) num2str(y_position(9))};

        cb                          = colorbar('northoutside');
        cb.Label.String           	= 'BL-corr. spike count';
        cb.FontSize                 = 10;
    end
        
    cmap                        = [0 0 0;jet(256)]; colormap(cmap); shg;
end
