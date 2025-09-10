function [out, nRep] = RF_estimation_online(d,f,ax_lst,in,nRep)

if nargin < 5
    nRep                    = 0;
end

if nargin < 4
    for iChannel = 1:32
        in{iChannel}        = nan(9, 10, 1800);
    end
end

%% Process

idx.tStart                  = d.event == 'TRIAL_start';
idx.tEnd                    = d.event == 'TRIAL_end';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.type                    = d.event == 'TRIAL_type';
idx.spikes                  = d.event == 'IO_spikes';

trl.start                   = d.time(idx.tStart);
trl.end                     = d.time(idx.tEnd);

channel                     = 1:32;

x_position                  = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_position                  = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
win_offset                  = [50e3 250e3];

% Trial index
trl.start
trlIdx                      = d.time >= trl.start(end) & d.time <= trl.end(end);
% trlIdx                      = d.time >= trl.start & d.time <= trl.end;

% Trial outcome
outcome                     = getTrialData(d.value, trlIdx, idx.outcome);
ttype                       = getTrialData(d.value, trlIdx, idx.type);
rdp_time                    = getTrialData(d.time, trlIdx, idx.type);

% Abort if fixation break
if strcmp(outcome, 'fixation break')
    out = in;
    return
end

nRep                    = nRep+1;

% Baseline spiking
blIdx                   = [];
blIdx                   = d.time >= trl.start(end) & d.time <= rdp_time(1);
bl_spike_time           = getTrialData(d.time, blIdx, idx.spikes);
bl_spike_id          	= getTrialData(d.value, blIdx, idx.spikes);
bl_spike_id             = floor(bl_spike_id); % channel-wise for now

% Stimulus position loop
for iStim = 1:length(ttype)
    
    % Stimulus index - offset by certain lag
    sIdx                = [];
    sIdx            	= d.time >= rdp_time(iStim)+win_offset(1) & d.time <= rdp_time(iStim)+win_offset(2);
    %     dur                 = (rdp_time(iStim)+win_offset(2) - rdp_time(iStim)+win_offset(1)) / 1e6;
    
    % Stimulus position
    pos                 = split(ttype{iStim},'_');
    rdp_x               = str2num(pos{1}(2:end));
    rdp_y               = str2num(pos{2}(2:end));
    
    % Spike count
    spike_time        	= getTrialData(d.time, sIdx, idx.spikes);
    spike_id          	= getTrialData(d.value, sIdx, idx.spikes);
    spike_id            = floor(spike_id); % units pooled - channel-wise for now
    
    % Unit loop
    for iChan = 1:length(channel)
        bl_spikes       = length(bl_spike_time(bl_spike_id == channel(iChan)));
        nSpikes         = length(spike_time(spike_id == channel(iChan)));
        spk{iChan}      = spike_time(spike_id == channel(iChan));
        
        % Stimulus surface represents baseline-corrected spike count
        mat                 = zeros(length(y_position), length(x_position));
        mat(rdp_y == y_position, rdp_x == x_position) = nSpikes - bl_spikes;
        %         mat(rdp_y == y_position, rdp_x == x_position) = ( (nSpikes - bl_spikes) / dur);
        %         in{iUnit}          = cat(3,in{iUnit}, mat);
        in{iChan}(:,:,nRep)	= mat;
    end
end

out = in;

%% Plot %%%

% Every N-th trial...
if mod(cell2mat(d.value(idx.tStart)),5) == 0
    %     close all
    %     f                               = figure('units','normalized','position',[0 0 1 1]);
    figure(f);
    
    for iChan = 1:size(out,2)
        
        if sum(sum(nansum(out{iChan},3))) == 0
            continue
        end
        
        ax = ax_lst(iChan);
        
        for iUnit = 1:size(out,2)
            
            % Continue if no spikes detected
            if sum(sum(nansum(out{iUnit},3))) == 0
                continue
            end
            
            % Assign axis
            axes(ax_lst(iUnit));
            
            % Plot data
            spk_sum                     = nansum(out{iUnit},3);
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
end
