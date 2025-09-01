% This function contains a summary of all important details for each
% experiment.
function par = getRecParam(animalID)

par = struct();
par.pth.pl2         = ['/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_electrophysiology/' animalID '/pl2/'];       % Source directory - Plexon [pl2]
par.pth.h5          = ['/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_electrophysiology/' animalID '/h5/'];        % Source directory - Behavior [h5]
par.pth.spikes      = ['/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_electrophysiology/' animalID '/sorted/'];    % Source directory - Spikes [mat/plx]
par.pth.dest_dir  	= ['/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_electrophysiology/' animalID '/summary/'];   % Destination directory for output [png/mat/pdf]
c                   = 1; % Session counter

switch animalID
    case 'Nilan'
        %%% REC_NIL_005_20250221 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_006_20250228 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_007_20250305 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_008_20250306 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_009_20250307 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sid                        	= ['rec' num2str(c)];           % Session ID
        par.(sid).date           	= '20250307';                 	% Date of recording
        par.(sid).fname_pl2        	= '20250710_nil_CPR_block1_phy4_rec009_fxs.pl2';    % File name of raw Plexon data [pl2]
        par.(sid).fname_h5        	= '20250307_nil_CPR_block1_phy4_fxs.h5';            % File name of behavior/spikes [h5]
        par.(sid).exp              	= {'RFmapping','CPRsolo'};      % Experiment with dyad ID, e.g. {'RFmapping','CPRsolo','CPRdyadic_nilann','CPRdyadic_nilfxs'}
        par.(sid).nProbe           	= 1;                            % Number of V-Probes
        par.(sid).probe_id         	= {'red_dot',''};               % Probe ID
        par.(sid).tip_guide_tube	= [1.0];                        % Distance [mm] of guide tube below dura [probe1 probe 2]
        par.(sid).tip_electrode   	= [6.6];                        % Distance [mm] of electrode tip below guide tube tip [probe1 probe 2] -> aligned recording position
        par.(sid).probe_coord      	= {[0.5 -2.5],[]};              % Grid coordinates of probes {[x1 y1] [x2 y2]}
        
        par                         = add_information(par, sid);  	% Add channel-wise information for each probe - see function below
        c                       	= c+1;                          % Add to session counter

        %%% REC_NIL_033_20250627 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sid                        	= ['rec' num2str(c)];           % Session ID
        par.(sid).date           	= '20250627';                  	% Date of recording
        par.(sid).fname_pl2        	= '20250627_nil_CPR_block1_phy4_rec033_ann.pl2';    % File name of raw Plexon data [pl2]
        par.(sid).fname_h5        	= '20250627_nil_CPR_block1_phy4_ann.h5';            % File name of behavior/spikes [h5]
        par.(sid).exp              	= {'RFmapping','CPRsolo'};      % Experiment with dyad ID, e.g. {'RFmapping','CPRsolo','CPRdyadic_nilann','CPRdyadic_nilfxs'}
        par.(sid).nProbe           	= 2;                            % Number of V-Probes
        par.(sid).probe_id         	= {'green_cross','blue dot'};  	% Probe ID
        par.(sid).tip_guide_tube	= [3.0 2.5];                    % Distance [mm] of guide tube below dura [probe1 probe 2]
        par.(sid).tip_electrode   	= [5.9 6.1];                   	% Distance [mm] of electrode tip below guide tube tip [probe1 probe 2] -> aligned recording position
        par.(sid).probe_coord      	= {[5.5 -2.5],[-5.5 1.5]};    	% Grid coordinates of probes {[x1 y1] [x2 y2]}
        
        par                         = add_information(par, sid);  	% Add channel-wise information for each probe - see function below
        c                       	= c+1;                          % Add to session counter
        
        %%% REC_NIL_034_20250703 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sid                        	= ['rec' num2str(c)];           % Session ID
        par.(sid).date           	= '20250703';                  	% Date of recording
        par.(sid).fname_pl2        	= '20250703_nil_CPR_block1_phy4_rec034_fxs.pl2';    % File name of raw Plexon data [pl2]
        par.(sid).fname_h5        	= '20250703_nil_CPR_block1_phy4_fxs.h5';            % File name of behavior/spikes [h5]
        par.(sid).exp              	= {'RFmapping','CPRsolo'};      % Experiment with dyad ID, e.g. {'RFmapping','CPRsolo','CPRdyadic_nilann','CPRdyadic_nilfxs'}
        par.(sid).nProbe           	= 1;                            % Number of V-Probes
        par.(sid).probe_id         	= {'green_cross',''};           % Probe ID
        par.(sid).tip_guide_tube	= [3.0];                        % Distance [mm] of guide tube below dura [probe1 probe 2]
        par.(sid).tip_electrode   	= [7.0];                        % Distance [mm] of electrode tip below guide tube tip [probe1 probe 2] -> aligned recording position
        par.(sid).probe_coord      	= {[2.5 -2.5],[]};              % Grid coordinates of probes {[x1 y1] [x2 y2]}
        
        par                         = add_information(par, sid);  	% Add channel-wise information for each probe - see function below
        c                       	= c+1;                          % Add to session counter
        
        %%% REC_NIL_035_20250710 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sid                        	= ['rec' num2str(c)];           % Session ID
        par.(sid).date           	= '20250710';                  	% Date of recording
        par.(sid).fname_pl2        	= '20250710_nil_CPR_block1_phy4_rec035_ann.pl2';    % File name of raw Plexon data [pl2]
        par.(sid).fname_h5        	= '20250710_nil_CPR_block1_phy4_ann.h5';            % File name of behavior/spikes [h5]
        par.(sid).exp              	= {'RFmapping','CPRsolo'};      % Experiment with dyad ID, e.g. {'RFmapping','CPRsolo','CPRdyadic_nilann','CPRdyadic_nilfxs'}
        par.(sid).nProbe           	= 1;                            % Number of V-Probes
        par.(sid).probe_id         	= {'green_cross',''};           % Probe ID
        par.(sid).tip_guide_tube	= [4.0];                        % Distance [mm] of guide tube below dura [probe1 probe 2]
        par.(sid).tip_electrode   	= [5.0];                        % Distance [mm] of electrode tip below guide tube tip [probe1 probe 2] -> aligned recording position
        par.(sid).probe_coord      	= {[4.5 -3.5],[]};              % Grid coordinates of probes {[x1 y1] [x2 y2]}
        
        par                         = add_information(par, sid);  	% Add channel-wise information for each probe - see function below
        c                       	= c+1;                          % Add to session counter
        
        %%% REC_NIL_036_20250711 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        sid                        	= ['rec' num2str(c)];           % Session ID
        par.(sid).date           	= '20250711';                  	% Date of recording
        par.(sid).fname_pl2        	= '20250711_nil_CPR_block1_phy4_rec036_ann.pl2';    % File name of raw Plexon data [pl2]
        par.(sid).fname_h5        	= '20250711_nil_CPR_block1_phy4_ann.h5';            % File name of behavior/spikes [h5]
        par.(sid).exp              	= {'RFmapping','CPRsolo'};      % Experiment with dyad ID, e.g. {'RFmapping','CPRsolo','CPRdyadic_nilann','CPRdyadic_nilfxs'}
        par.(sid).nProbe           	= 2;                            % Number of V-Probes
        par.(sid).probe_id         	= {'green_cross','blue_dot'};   % Probe ID
        par.(sid).tip_guide_tube	= [3.5 3.5];                    % Distance [mm] of guide tube below dura [probe1 probe 2]
        par.(sid).tip_electrode   	= [2.6 4.0];                    % Distance [mm] of electrode tip below guide tube tip [probe1 probe 2] -> aligned recording position
        par.(sid).probe_coord      	= {[5.5 -2.5],[-4.5 0.5]};      % Grid coordinates of probes {[x1 y1] [x2 y2]}
        
        par                         = add_information(par, sid);  	% Add channel-wise information for each probe - see function below
        c                       	= c+1;                          % Add to session counter
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 'Dunkin'
        
        
end
end

% Not sure this is needed anymore. We might be able to automise
% this based on the output from the spike sorter
%         cluster                                     = {                     % fill in [unit label]
%             {'MU'};... % Channel 1
%             {'MU'};... % Channel 2
%             {'MU','MU'};... % Channel 3
%             {'MU'};... % Channel 4
%             {'SU'};... % Channel 5
%             {''};... % Channel 6
%             {'SU'};... % Channel 7
%             {''};... % Channel 8
%             {'MU'};... % Channel 9
%             {'MU', 'SU'};... % Channel 10
%             {'MU','MU'};... % Channel 11
%             {'MU'};... % Channel 12
%             {'MU'};... % Channel 13
%             {'SU'};... % Channel 14
%             {'SU'};... % Channel 15
%             {'MU'};... % Channel 16
%             {'MU'};... % Channel 17
%             {'MU','MU'};... % Channel 18
%             {'SU'};... % Channel 19
%             {'MU','MU'};... % Channel 20
%             {''};... % Channel 21
%             {''};... % Channel 22
%             {'SU','SU','MU'};... % Channel 23
%             {'SU'};... % Channel 24
%             {'SU','MU'};... % Channel 25
%             {''};... % Channel 26
%             {'MU'};... % Channel 27
%             {'MU'};... % Channel 28
%             {''};... % Channel 29
%             {''};... % Channel 30
%             {''};... % Channel 31
%             {'MU'};... % Channel 32
%             };

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% HELPER FUNCTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function par = add_information(par, sid)

%%% Constant: V-Probe specs %%%
nChan                               = 32;	% Number of channels of V-Probe
offset                              = .3; 	% Distance tip to 1st contact [mm]
spacing                             = .1; 	% Inter-electrode spacing [mm]

%%% Add specific information for each V-probe %%
for iProbe = 1:par.(sid).nProbe
    pid                             = ['probe' num2str(iProbe)];        % Probe number
    par.(sid).(pid).probeID         = par.(sid).probe_id{iProbe};       % Probe ID
    par.(sid).(pid).probe_coord     = par.(sid).probe_coord{iProbe};    % Grid coordniates of probe
    
    % Channel numbers [superficial deep]
    if iProbe == 1
        par.(sid).(pid).nChan   	= 1:nChan;
    elseif iProbe == 2
        par.(sid).(pid).nChan    	= (1:nChan) + nChan;
    end
    
    % Calculate contact position
    contact_position                = fliplr(offset:spacing:(spacing*(nChan-1)+offset)); % Contact distance to tip of electrode [mm]
    contact_depth                   = (par.(sid).tip_guide_tube(iProbe) + par.(sid).tip_electrode(iProbe)) - contact_position; % Contact depth from dura [mm]
    par.(sid).(pid).array_depth     = contact_depth; % [superficial deep]
    
    for iChan = par.(sid).(pid).nChan
        cid                             = ['ch' num2str(iChan)]; % Channel ID
        par.(sid).(pid).(cid).chan_num  = iChan;
        
        % Assign recording location
        if iChan < 33
            par.(sid).(pid).(cid).coord	= [par.(sid).probe_coord{iProbe} contact_depth(iChan)]; % Center of grid at level of dura [0x 0y 0z]mm
            par.(sid).(pid).(cid).field	= 'MT/MST';
        else
            par.(sid).(pid).(cid).coord	= [par.(sid).probe_coord{iProbe} contact_depth(iChan-nChan)]; % Center of grid at level of dura [0x 0y 0z]mm
            par.(sid).(pid).(cid).field	= 'IPS';
        end
    end
end
end