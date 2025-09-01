
function par = getRecParam(animalID)

par = struct();
par.spath  	= ['Y:\EPHYS\RAWDATA\NHP\Neuralynx\FigureGround\' animalID '\'];    % Source directory [pl2/log-files]
par.dpath  	= ['C:\Rec\' animalID '\'];                                         % Data directory [mat/png-files]
par.dpath  	= ['D:\Rec\' animalID '\'];                                     % Data directory [mat/png-files]
c          	= 1;                                                                % Counter

switch animalID
    
    case 'Nilan'
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Constant %%%
        nChan                                       = 32;                   % Number of channels of V-Probe
        offset                                      = .3;                   % Distance tip to 1st contact [mm]
        spacing                                     = .1;                   % Inter-electrode spacing [mm]
        
        %%% REC_NIL_005_20250221 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_006_20250228 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_007_20250305 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_008_20250306 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_009_20250307 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        par.date                                    = '20250307';                                   % Date of recording
        par.fname_pl2                               = 'fxs-CPR20250307_nil_009-01+01_sorted.pl2';   % File name of sorted Plexon file
        par.fname_h5                                = '20250307_nil_CPRsolo_block1_phy4.h5';        % File name of combined data HDF5 file  
  
        par.(['rec' num2str(c)]).probe_id           = 'red_dot';        	% Probe ID
        par.(['rec' num2str(c)]).tip_guide_tube   	= 1.0;                  % fill in [mm below dura]
        par.(['rec' num2str(c)]).tip_electrode     	= 6.6;                  % fill in [mm below guide tube tip]
        par.(['rec' num2str(c)]).nChan              = 1:32;               	% [superficial deep]
        
        %%% Calculate electrode position %%%
        contact_position                          	= fliplr(offset:spacing:(spacing*(nChan-1)+offset)); % Contact distance to tip of electrode [mm]
        contact_depth                               = tip_guide_tube + (tip_electrode-contact_position); % Contact depth from dura [mm]
        par.(['rec' num2str(c)]).array_depth        = contact_depth;      	% [superficial deep]

        cluster                                     = {                     % fill in [unit label]
            {'MU'};... % Channel 1
            {'MU'};... % Channel 2
            {'MU','MU'};... % Channel 3
            {'MU'};... % Channel 4
            {'SU'};... % Channel 5
            {''};... % Channel 6
            {'SU'};... % Channel 7
            {''};... % Channel 8
            {'MU'};... % Channel 9
            {'MU', 'SU'};... % Channel 10
            {'MU','MU'};... % Channel 11
            {'MU'};... % Channel 12
            {'MU'};... % Channel 13
            {'SU'};... % Channel 14
            {'SU'};... % Channel 15
            {'MU'};... % Channel 16
            {'MU'};... % Channel 17
            {'MU','MU'};... % Channel 18
            {'SU'};... % Channel 19
            {'MU','MU'};... % Channel 20
            {''};... % Channel 21
            {''};... % Channel 22
            {'SU','SU','MU'};... % Channel 23
            {'SU'};... % Channel 24
            {'SU','MU'};... % Channel 25
            {''};... % Channel 26
            {'MU'};... % Channel 27
            {'MU'};... % Channel 28
            {''};... % Channel 29
            {''};... % Channel 30
            {''};... % Channel 31
            {'MU'};... % Channel 32
            };
 
        %%% Add channel-wise information %%%
        for iChan = 1:32
            par.(['rec' num2str(c)]).(['ch' num2str(iChan)]).coord   = [0.5 -2.5 contact_depth(iChan)]; % Center of grid at level of dura [0x 0y 0z]mm
            par.(['rec' num2str(c)]).(['ch' num2str(iChan)]).field   = 'MT/MST';
            par.(['rec' num2str(c)]).(['ch' num2str(iChan)]).clus    = cluster{iChan};
        end
        c                                           = c+1;                  % Count sessions
                
        %%% REC_NIL_010_20250312 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_011_20250313 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_012_20250314 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_013_20250320 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% REC_NIL_014_20250321 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

                
    case 'Dunkin'
        
        
end