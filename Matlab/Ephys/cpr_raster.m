%% IMPORT
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/

pth                         = '/Users/fschneider/Desktop/EPHYS/pilot_data/';
fname_cpr                   = 'fxs-CPR_20230707-cla-016-01+01.mwk2';
fname_tun                  	= 'fxs-RF_20230707-cla-016-01+02.mwk2';
fname_splt                  = split(fname_cpr,'_');
date                        = fname_splt{2}(1:8);
rec_pos                     = fname_splt{2}(18:22);
dest_dir                    = ['/Users/fschneider/Desktop/EPHYS/RF/' date '/'];

mkdir(dest_dir)

var_import = {
    'INFO_', ...
    'TRIAL_', ...
    'IO_spikes',...
    'IO_joystickStrength',...
    '#stimDisplay'};

cpr                         = CPR_import_mwk2([pth fname_cpr], var_import, false);
tun                      	= CPR_import_mwk2([pth fname_tun], var_import, false);

%% SELECT, PROCESS & PLOT DATA

%%%% Filter Joystick respons ein same direction
%%%% Filter for hits

close all
tstep               = .001;	% Resolution for SDF [s]
nbins               = 50;
PD_offset           = 30;
 
f                   = figure; hold on
ff = figure; hold on

% For each channel
for iChan = 27
    cpr_spikes                  = extractSpikes_CPR(cpr, iChan);
    tun_spikes                  = extractSpikes_tuning(tun, iChan);
    
    snr                         = unique(cpr_spikes.coh_val);
    dirs                        = unique(tun_spikes.dir);
    col                         = hsv(length(snr));

    % Average FR per direction
    for iDir = 1:length(dirs)
        dIdx                    = tun_spikes.dir == dirs(iDir);
        avg_fr(iDir)            = mean(tun_spikes.spikes_hz(dIdx));
    end
    
    % Skip channel if FR < 5Hz
    if sum(avg_fr < 5) == length(avg_fr)
        disp(['Skipped channel ' num2str(iChan)])
        continue
    end
    
    %%% PLOT TRIAL ONSET RESPONSE    
    for iExp = 1:2   
        clear spike_times
        
        if iExp == 1
            spike_times	= cpr_spikes.spikes_start;
            exp_str = ['Chan' num2str(iChan) ' : CPR'];
            time         	= tstep-.25:tstep:.25;       % Time vector
        else
            spike_times	= tun_spikes.spikes_start;
            exp_str = ['Chan' num2str(iChan) ' : TUN'];
            time          	= tstep:tstep:.5;
        end
        
        figure(1)
        ax                  = subplot(3,2,iExp);hold on
        [all_spks, sdf]    	= FR_estimation(spike_times, time);
        
        ax.Title.String     = [exp_str ' TrialON'];
        ax.YLim             = [0 length(spike_times)];
        ax.FontSize         = 14;
        ax.XTick            = [-.2 0 .2 .4];
        ax.YLabel.String  	= 'Trial';
        
        if iExp == 1
            ax.XLim       	= [-.25 .25];
        else
            ax.XLim        	= [0 .5];
        end
        
        ax                  = subplot(3,2,iExp+2);hold on
        h                   = histogram(all_spks,nbins);
        h.FaceColor         = 'k';
        
        mVal_hist           = max(h.Values)+round(max(h.Values)*.1);
        ax.YLim             = [0 mVal_hist];
        ax.XTick            = [-.2 0 .2 .4];
        ax.XLabel.String  	= 'Time [s]';
        ax.YLabel.String  	= 'Spks/Bin';
        ax.FontSize         = 14;
        
        if iExp == 1
            ax.XLim       	= [-.25 .25];
        else
            ax.XLim        	= [0 .5];
        end
        
        ax                  = subplot(3,2,iExp+4);hold on
        pl                  = plot(mean(sdf), 'Color', 'k', 'LineWidth', 1.5);
        mVal_sdf            = max(mean(sdf)) + round(max(mean(sdf))*.1);
        ax.XLim             = [0 500];
        ax.YLim             = [0 mVal_sdf];
        ax.XLabel.String  	= 'Time [s]';
        ax.YLabel.String  	= 'FR [Hz]';
        ax.FontSize         = 14;
        
        if iExp == 1
            ax.XTick      	= [50 250 450];
            ax.XTickLabel  	= {'-.2', '0', '0.2'};
        else
            ax.XTick      	= [1 200 400];
            ax.XTickLabel  	= {'0', '0.2', '0.4'};
        end
    end
    
    print(f, [dest_dir date '_chan' num2str(iChan) '_TrlON'], '-r500', '-dpng');
    
    %%% PLOT TUNING RASTER
    time                    = tstep-.05:tstep:.25;      % Time vector    
    for iExp = 1:2
        f                   = figure;

        for iDir = 1:2
            
            figure(f)
            ax                      = subplot(3,2,iDir);hold on
            
            if iExp == 1 %%% Tuning
                exp_str          	= 'Tun';

                if iDir == 1
                    pd_idx        	= tun_spikes.dir == dirs(avg_fr == max(avg_fr));
                    sdir            = dirs(avg_fr == max(avg_fr));
                else
                    pd_idx       	= tun_spikes.dir == dirs(avg_fr == min(avg_fr));
                    sdir            = dirs(avg_fr == min(avg_fr));
                end
                spike_times         = tun_spikes.spikes_aligned(pd_idx);
                pt                  = patch([.0 .2 .2 .0], [0 0 length(spike_times) length(spike_times)], [.5 .5 .5], 'FaceAlpha',.2, 'EdgeColor','none');

            elseif iExp == 2 %%% CPR
                exp_str          	= 'CPR';
                dir_idx             = [];

                if iDir == 1
                    maxFR           = dirs(avg_fr == max(avg_fr));          % Direction of interest
                    sdir            = maxFR;
                    trigger_min  	= maxFR - PD_offset + 360;              % Window around direction
                    trigger_max   	= maxFR + PD_offset + 360;
                else
                    minFR           = dirs(avg_fr == min(avg_fr));
                    sdir            = minFR;
                    trigger_min  	= minFR - PD_offset + 360;
                    trigger_max   	= minFR + PD_offset + 360;
                end
                              
                for i = 1:3
                    dir_idx(i,:)    = (cpr_spikes.dir_val+360) >= trigger_min*i & (cpr_spikes.dir_val+360) <= trigger_max*i;
                end
                
                pd_idx              = logical(sum(dir_idx)) & cpr_spikes.coh_val == .9;
                spike_times         = cpr_spikes.spikes_aligned(pd_idx);
                pt                  = patch([.0 .25 .25 .0], [0 0 length(spike_times) length(spike_times)], [.5 .5 .5], 'FaceAlpha',.2, 'EdgeColor','none');
            end        

            [all_spks, sdf]         = FR_estimation(spike_times, time);
           
            if iDir == 1
                ax.Title.String     = ['max(FR): ' num2str(sdir) 'deg'];
            else
                ax.Title.String     = ['min(FR): ' num2str(sdir) 'deg'];
            end
            
            ax.YLim             = [0 length(spike_times)];
            ax.XLim             = [-.05 .25];
            ax.FontSize         = 14;
            ax.XTick            = [0 .1 .2];
            ax.YLabel.String  	= 'State';
            
            ax                  = subplot(3,2,iDir+2);hold on
            h                   = histogram(all_spks,nbins);
            h.FaceColor         = 'k';
            
            if iDir == 1
                mVal_hist     	= max(h.Values)+round(max(h.Values)*.1);
            end
            
            ax.XLim             = [-.05 .25];
            ax.YLim             = [0 mVal_hist];
            ax.XTick            = [0 .1 .2];
            ax.XLabel.String  	= 'Time [s]';
            ax.YLabel.String  	= 'Spks/Bin';
            ax.FontSize         = 14;
            
            ax                  = subplot(3,2,iDir+4);hold on
            pl                  = plot(mean(sdf), 'Color', 'k', 'LineWidth', 1.5);
            
            if iDir == 1
                mVal_sdf       	= max(mean(sdf)) + round(max(mean(sdf))*.1);
            end
            
            ax.XLim             = [0 300];
            ax.YLim             = [0 mVal_sdf];
            ax.XTick            = [50 150 250];
            ax.XTickLabel       = {'0', '0.1', '0.2'};
            ax.XLabel.String  	= 'Time [s]';
            ax.YLabel.String  	= 'FR [Hz]';
            ax.FontSize         = 14;
            
            
            %%%% ECC - FR relationship
            if iExp == 2
                figure(ff)
                ax = subplot(2,1,iDir);hold on
                for iCoh = 1:length(snr)
                    cohIdx = cpr_spikes.coh_val ==  snr(iCoh);
                    sc = scatter(cellfun(@length,cpr_spikes.spikes_stEnd(sum(dir_idx) & cohIdx)),cellfun(@mean,cpr_spikes.js_ecc(sum(dir_idx) & cohIdx)));
                    sc.MarkerFaceColor = col(iCoh,:);
                    sc.MarkerEdgeColor = 'none';
                    sc.MarkerFaceAlpha = .7;
                    ls = lsline;
                end
                
                %             sc = scatter(cellfun(@length,cpr_spikes.spikes_stEnd(pd_idx)),cellfun(@mean,cpr_spikes.js_ecc(pd_idx)));
                %             sc.MarkerFaceColor = col(2,:);
                %             sc.MarkerEdgeColor = 'none';
                %             sc.MarkerFaceAlpha = .4;
                %             ls = lsline;
                
                set(ls(1),'color', col(1,:))
                set(ls(2),'color', col(2,:))
                set(ls(3),'color', col(3,:))
                
                ax.FontSize = 14;
                ax.XLabel.String = 'Spike count';
                ax.YLabel.String = 'Avg joystick eccentricity';
                
                if iDir == 1
                    ax.Title.String     = ['max(FR): ' num2str(sdir) 'deg'];
                else
                    ax.Title.String     = ['min(FR): ' num2str(sdir) 'deg'];
                end
            end
        end
        print(f, [dest_dir date '_chan' num2str(iChan) '_' exp_str '_raster'], '-r500', '-dpng');
    end
    
    %%% COMPARE TUNING TO CPR
    f                           = figure;
    ax                          = polaraxes;
    pp                          = polarplot(deg2rad(dirs),avg_fr);
    pp.Color                    = [0 0 0];
    pp.LineWidth                = 3;
    ax.ThetaZeroLocation      	= 'top';
    ax.ThetaDir               	= 'clockwise'; % 90 degrees at the right
    ax.FontSize                 = 14;
    ax.Title.String             = ['Chan: ' num2str(iChan)];
    hold on
    
    for iCoh = 1:length(snr)
        cohIdx                  = cpr_spikes.coh_val == snr(iCoh);
        spks                    = cpr_spikes.spikes_hz(cohIdx);
        dir_rad                 = deg2rad(cpr_spikes.dir_val(cohIdx));
        
        ps(iCoh)                    = polarscatter(dir_rad,spks);
        ps(iCoh).MarkerFaceColor    = col(iCoh,:);
        ps(iCoh).MarkerFaceAlpha    = .5;
        ps(iCoh).MarkerEdgeColor    = 'none';
    end
    
    legend([ps pp],{'CPR 20%','CPR 40%','CPR 90%','Tun 90%'})
    
    print(f, [dest_dir date '_chan' num2str(iChan) '_FR_comparison'], '-r500', '-dpng');
%     close all
end
%% FUNCIONS

function out = extractSpikes_CPR(in, chan)

idx.tStart                  = in.event == 'TRIAL_start';
idx.tStart                  = in.event == 'TRIAL_start';
idx.tEnd                    = in.event == 'TRIAL_end';
idx.outcome                 = in.event == 'TRIAL_outcome';
idx.type                    = in.event == 'TRIAL_type';
idx.spikes                  = in.event == 'IO_spikes';
idx.js_ecc                  = in.event == 'IO_joystickStrength';

trl.start                   = in.time(idx.tStart);
trl.end                     = in.time(idx.tEnd);

win_offset                  = [50e3 250e3];

nSpikes                     = [];
cnt                         = 0;
dir_val                 	= [];
coh_val                   	= [];

for iTrl = 1:length(trl.start)
    
    % Trial index
    trlIdx                  = [];
    trlIdx                  = in.time >= trl.start(iTrl) & in.time <= trl.end(iTrl);
    
    start_idx               = in.time >= trl.start(iTrl)-win_offset(2) & in.time <= trl.start(iTrl)+win_offset(2);
    spike_time_start        = getTrialData(in.time, start_idx, idx.spikes);
    spike_id_start          = floor(getTrialData(in.value, start_idx, idx.spikes));
    out.spikes_start{iTrl}  = (spike_time_start(spike_id_start == chan) - trl.start(iTrl)) ./ 1e3;
    
    % Trial outcome
    outcome              	= getTrialData(in.value, trlIdx, idx.outcome);
    ttype                   = getTrialData(in.value, trlIdx, idx.type);
    rdp_time                = getTrialData(in.time, trlIdx, idx.type);
    
    % Stimulus position loop
    for iStim = 1:length(ttype)
        
        % Stimulus index - offset by certain lag
        sIdx                = [];
        sIdx            	= in.time >= rdp_time(iStim)+win_offset(1) & in.time <= rdp_time(iStim)+win_offset(2);
        dur                 = (rdp_time(iStim)+win_offset(2) - rdp_time(iStim)+win_offset(1)) / 1e6;
        
        % Spike count
        spike_time        	= getTrialData(in.time, sIdx, idx.spikes);
        spike_id          	= getTrialData(in.value, sIdx, idx.spikes);
        spike_id            = floor(spike_id); % channel-wise for now
        
        tIdx            	= in.time >= rdp_time(iStim)-win_offset(1) & in.time <= rdp_time(iStim)+win_offset(2);
        spike_time_trans    = getTrialData(in.time, tIdx, idx.spikes);
        spike_id_trans      = floor(getTrialData(in.value, tIdx, idx.spikes));
     
                
        if iStim ~= length(ttype)
            evt             = rdp_time(iStim+1);
        else
            evt             = trl.end(iTrl);
        end
        
        stEndIdx            = in.time >= evt-1000e3 & in.time <= evt;
        js_ecc_stEnd       	= getTrialData(in.value, stEndIdx, idx.js_ecc);
        spike_time_stEnd   	= getTrialData(in.time, stEndIdx, idx.spikes);
        spike_id_stEnd     	= floor(getTrialData(in.value, stEndIdx, idx.spikes));
        
        cnt                     = cnt+1;
        tmp                     = split(ttype{iStim},'_');
        out.spikes_sum(cnt)    	= sum(spike_id == chan);
        out.spikes_hz(cnt)     	= sum(spike_id == chan)  / dur;
        out.spikes_aligned{cnt} = (spike_time_trans(spike_id_trans == chan) - rdp_time(iStim)) ./ 1e3;
        out.spikes_stEnd{cnt}   = (spike_time_stEnd(spike_id_stEnd == chan) - evt) ./ 1e3;
        out.js_ecc{cnt}         = js_ecc_stEnd;
        out.dir_val(cnt)      	= str2num(tmp{2}(4:end));
        out.coh_val(cnt)      	= str2num(tmp{1}(4:end));
    end
end
end

function out = extractSpikes_tuning(in, chan)

idx.task                    = in.event == 'INFO_task';
idx.tStart                  = in.event == 'TRIAL_start';
idx.tEnd                    = in.event == 'TRIAL_end';
idx.outcome                 = in.event == 'TRIAL_outcome';
idx.type                    = in.event == 'TRIAL_type';
idx.spikes                  = in.event == 'IO_spikes';

trl.start                   = in.time(idx.tStart);
trl.end                     = in.time(idx.tEnd);

tmp_task                    = in.value(idx.task);
tmp_task_ts                 = in.time(idx.task);

for iTrl = 1:length(trl.start)
    trl.task{iTrl}          = tmp_task{find(tmp_task_ts > trl.start(iTrl),1,'first')};
end

win_offset                  = [50e3 250e3];
cnt                         = 0;
dir                         = [];
nSpikes                     = [];

for iTrl = 1:length(trl.start)
    
    if ~strcmp(trl.task{iTrl}, 'RF_tuning')
        continue
    end
    
    % Trial index
    trlIdx                  = [];
    trlIdx                  = in.time >= trl.start(iTrl) & in.time <= trl.end(iTrl);
    
    start_idx               = in.time >= (trl.start(iTrl)) & in.time <= (trl.start(iTrl)+500e3);
    spike_time_start        = getTrialData(in.time, start_idx, idx.spikes);
    spike_id_start          = floor(getTrialData(in.value, start_idx, idx.spikes));
    out.spikes_start{iTrl}  = (spike_time_start(spike_id_start == chan) - trl.start(iTrl)) ./ 1e3;
    
    % Trial outcome
    outcome              	= getTrialData(in.value, trlIdx, idx.outcome);
    ttype                   = getTrialData(in.value, trlIdx, idx.type);
    rdp_time                = getTrialData(in.time, trlIdx, idx.type);
    
    if strcmp(outcome, 'fixation break')
        nStim            	= length(ttype) - 1; % Exclude last stimulus
    else
        nStim             	= length(ttype);
    end
    
    % Stimulus position loop
    for iStim = 1:nStim
        
        % Stimulus index - offset by certain lag
        sIdx                = [];
        tIdx                = [];
        sIdx                = in.time >= rdp_time(iStim)+win_offset(1) & in.time <= rdp_time(iStim)+win_offset(2);
        tIdx                = in.time >= rdp_time(iStim)-win_offset(1) & in.time <= rdp_time(iStim)+win_offset(2);
        dur                 = (rdp_time(iStim)+win_offset(2) - rdp_time(iStim)+win_offset(1)) / 1e6;
        
        % Spike count
        spike_time        	= getTrialData(in.time, sIdx, idx.spikes);
        spike_id          	= floor(getTrialData(in.value, sIdx, idx.spikes));
        
        spike_time_trans  	= getTrialData(in.time, tIdx, idx.spikes);
        spike_id_trans     	= floor(getTrialData(in.value, tIdx, idx.spikes));
        
        cnt                     = cnt+1;
        out.spikes_sum(cnt)     = sum(spike_id == chan);
        out.spikes_hz(cnt)      = sum(spike_id == chan) ./ dur;
        out.spikes_aligned{cnt} = (spike_time_trans(spike_id_trans == chan) - rdp_time(iStim)) ./ 1e3;
        out.dir(cnt)            = str2num(ttype{iStim}(4:end));
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

function [all, sdf] = FR_estimation(spike_times, time)

sdf                 = [];
all                 = [];

for iTrial = 1:length(spike_times)
    
    spks            = spike_times{iTrial} ./1e3;  	% Get all spikes of respective trial
    all             = [all spks];                   % Concatenate spikes of all trials
    xspikes         = repmat(spks,3,1);             % Replicate array
    yspikes      	= nan(size(xspikes));           % NaN array
    
    if ~isempty(yspikes)
        yspikes(1,:) = iTrial-1;                   	% Y-offset for raster plot
        yspikes(2,:) = iTrial;
    end
    
    % Plot trial raster
    pl           	= plot(xspikes, yspikes, 'Color', 'k', 'LineWidth',1.25);
    
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