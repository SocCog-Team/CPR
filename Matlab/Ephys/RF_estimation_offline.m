%% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing/

clear all
close all

%% Import
pth                         = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_electrophysiology/Nilan/';
fname                       = '20250312_nil_CPRsolo_block1_phy4.mwk2';
fname_splt                  = split(fname,'_');
date                        = fname_splt{2}(4:end);
dest_dir                    = '/Users/fschneider/Desktop/EPHYS/RF/';

var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'RDP_',...
    'IO_fixation_flag',...
    'IO_spikes',...
    'EYE_x_dva',...
    'EYE_y_dva',...
    'IO_syncWord',...
    '#stimDisplay'};

% d                           = CPR_import_mwk2([pth '/mwk2' fname], var_import, true,'/Users/fschneider/Documents/GitHub/CPR/Matlab/PHY_Setup/server_code/felix_nhp_solo.cfg');
    d                       = MW_readH5([pth 'h5/' fname(1:end-5) '.h5']);
    
%% Process RF Mapping

idx.task                    = d.event == 'INFO_task';
idx.tStart                  = d.event == 'TRIAL_start';
idx.tEnd                    = d.event == 'TRIAL_end';
idx.outcome                 = d.event == 'TRIAL_outcome';
idx.type                    = d.event == 'TRIAL_type';
idx.spikes                  = d.event == 'IO_spikes';

trl.start                   = d.time(idx.tStart);
trl.end                     = d.time(idx.tEnd);

tmp_task                    = d.value(idx.task);
tmp_task_ts                 = d.time(idx.task);
for iTrl = 1:length(trl.start)
    trl.task{iTrl}          = tmp_task{find(tmp_task_ts > trl.start(iTrl),1,'first')};
end

spike_ID                   	= cell2mat(d.value(idx.spikes));
% units                       = unique(spike_ID);
% units                       = unique(floor(spike_ID)); % channel-wise for now
units                       = 1:32;

x_position                  = [-24 -21 -18 -15 -12 -9 -6 -3 0 3];
y_position                  = fliplr([-12 -9 -6 -3 0 3 6 9 12]);
win_offset                  = [50e3 100e3];
out                         = cell(1,length(units));
stim_id                     = reshape(1:(length(x_position)*length(y_position)),[length(y_position),length(x_position)]);
stim_cnt                    = 0;

% Trial loop
for iTrl = 1:length(trl.start)
    
    if ~strcmp(trl.task{iTrl}, 'RF_mapping')
        continue
    end
    
    % Trial index
    trlIdx                      = [];
    trlIdx                      = d.time >= trl.start(iTrl) & d.time <= trl.end(iTrl);
    
    % Trial outcome
    outcome                     = getTrialData(d.value, trlIdx, idx.outcome);
    ttype                       = getTrialData(d.value, trlIdx, idx.type);
    rdp_time                    = getTrialData(d.time, trlIdx, idx.type);
    
    if ~iscell(ttype)
        continue
    end
    
    if strcmp(outcome, 'fixation break')
        nStim                   = length(ttype) - 2; % Exclude last stimulus
    else
        nStim                   = length(ttype);
    end
    
    % Baseline spiking
    blIdx                     	= [];
    blIdx                       = d.time >= trl.start(iTrl) & d.time <= rdp_time(1);
    bl_spike_time               = getTrialData(d.time, blIdx, idx.spikes);
    bl_spike_id                 = getTrialData(d.value, blIdx, idx.spikes);
    bl_spike_id                 = floor(bl_spike_id); % channel-wise for now
    
    % Stimulus position loop
    for iStim = 2:nStim
        
        % Stimulus index - offset by certain lag
        sIdx                    = [];
        sIdx                    = d.time >= rdp_time(iStim-1)+win_offset(1) & d.time <= rdp_time(iStim)+win_offset(2);
  
        % Stimulus position
        pos                     = split(ttype{iStim-1},'_');
        rdp_x                   = str2num(pos{1}(2:end));
        rdp_y                   = str2num(pos{2}(2:end));
        stim_cnt                = stim_cnt+1;
        stim_pos(stim_cnt)      = stim_id(rdp_y == y_position, rdp_x == x_position);
        
        % Spike count
        spike_time              = getTrialData(d.time, sIdx, idx.spikes);
        spike_id                = getTrialData(d.value, sIdx, idx.spikes);
        spike_id                = floor(spike_id); % channel-wise for now
        
        % Unit loop
        for iUnit = 1:length(units)
            bl_spikes          	= length(bl_spike_time(bl_spike_id == units(iUnit)));
            nSpikes(iUnit,stim_cnt) = length(spike_time(spike_id == units(iUnit)));
        end
    end
end

%% PLOT: RF Mapping - raw and smoothed

mat_sum = nan(size(stim_id));

for iPlot = 1:2
    f                               = figure('units','normalized','position',[0 0 1 1]);
    units(units == 0)               = [];
    
    for iChan = units
        ax                          = subplot(4,8,iChan);
        
        if sum(units == iChan) > 0            
            for iPos = 1:max(stim_pos)
                mat_sum(stim_id == iPos) = sum(nSpikes(iChan, stim_pos == iPos));
            end
            
            if iPlot == 1
                hm                      = imagesc(mat_sum);
            else
                hm                      = imagesc(imgaussfilt(mat_sum,.75));
            end
            
            ax.XTick                = [1 5 9];
            ax.XTickLabel           = {num2str(x_position(1)) num2str(x_position(5)) num2str(x_position(9))};
            ax.YTick                = [1 5 9];
            ax.YTickLabel           = {num2str(y_position(1)) num2str(y_position(5)) num2str(y_position(9))};
        else
            imagesc(zeros(size(stim_id)))
            axis off
        end
        
        colormap(jet(256))
        ax.Title.String             = ['Chan: ' num2str(iChan)];
        ax.FontSize                 = 14;
        
    end
    
    if iPlot == 1
        print(f, [dest_dir fname '_rf_raw.png'], '-r300', '-dpng');
    else
        print(f, [dest_dir fname '_rf_smoothed.png'], '-r300', '-dpng');
    end
end

%% PLOT: Linear arrangement -> RF Mapping/Tuning

f                           = figure('units','normalized','position',[0 0 .3 3]);
height                      = linspace(.95,.01, 32);
mat_sum                     = nan(size(stim_id));
units(units == 0)           = [];

for iPlot = 1:2%3
    cnt = 0;
    for iChan = 1:32
        if iPlot == 1
            ax              = axes('Position', [.1 height(iChan) .05 .025]); hold on; axis off
        elseif iPlot == 2
            ax              = axes('Position', [.175 height(iChan) .05 .025]); hold on; axis off
        elseif iPlot == 3
            ax              = polaraxes('Position', [.25 height(iChan) .05 .025]); hold on
        end
        
        if sum(units == iChan) > 0
            cnt             = cnt+1;
            
            for iPos = 1:max(stim_pos)
                mat_sum(stim_id == iPos) = sum(nSpikes(iChan, stim_pos == iPos));
            end
            
            %%% Raw
            if iPlot == 1
                hm                 	= imagesc(flipud(mat_sum));
                colormap(jet(256))
                %%% Smoothed
            elseif iPlot == 2
                hm              	= imagesc(flipud(imgaussfilt(mat_sum,.75)));
                colormap(jet(256))
                %%% Tuning
            else
                for iDir = 1:length(dirs)
                    indx         	= dir == dirs(iDir) & chan == iChan;
                    sum_spks(iDir)	= mean(nSpikes_tun(indx));
                end
                
                pl                 	= polarplot(deg2rad(dirs),sum_spks);
                pl.LineWidth       	= 2;
                pl.Color          	= [0 0 0];
                ax.ThetaTick      	= [];
                ax.RTick          	= [];
                ax.ThetaZeroLocation      	= 'top';
                ax.ThetaDir               	= 'clockwise'; % 90 degrees at the right
            end
        else
            imagesc(zeros(9))
        end
    end
end

print(f, [dest_dir fname '_rf_array.png'], '-r300', '-dpng');
