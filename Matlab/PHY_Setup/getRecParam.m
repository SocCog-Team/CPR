
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
        
        %%% Variable %%%%
        tip_guide_tube                              = 1.5;                  % fill in
        tip_electrode                               = 7.7;                  % fill in 
        cluster                                     = {                     % fill in
            {'SU', 'MU'};... % Channel 1
            {'fill_me'};.... % Channel 2
            {'fill_me'};.... % Channel 3
            {'fill_me'};.... % Channel 4
            {'fill_me'};.... % Channel 5
            {'fill_me'};.... % Channel 6
            {'fill_me'};.... % Channel 7
            {'fill_me'};.... % Channel 8
            {'fill_me'};.... % Channel 9
            {'fill_me'};.... % Channel 10
            {'fill_me'};.... % Channel 11
            {'fill_me'};.... % Channel 12
            {'fill_me'};.... % Channel 13
            {'fill_me'};.... % Channel 14
            {'fill_me'};.... % Channel 15
            {'fill_me'};.... % Channel 16
            {'fill_me'};.... % Channel 17
            {'fill_me'};.... % Channel 18
            {'fill_me'};.... % Channel 19
            {'fill_me'};.... % Channel 20
            {'fill_me'};.... % Channel 20
            {'fill_me'};.... % Channel 21
            {'fill_me'};.... % Channel 22
            {'fill_me'};.... % Channel 23
            {'fill_me'};.... % Channel 24
            {'fill_me'};.... % Channel 25
            {'fill_me'};.... % Channel 26
            {'fill_me'};.... % Channel 27
            {'fill_me'};.... % Channel 28
            {'fill_me'};.... % Channel 29
            {'fill_me'};.... % Channel 30
            {'fill_me'};.... % Channel 31
            {'fill_me'};.... % Channel 32
            };
        
        %%% Calculate electrode position %%%
        contact_position                          	= fliplr(offset:spacing:(spacing*(nChan-1)+offset)); % Contact distance to tip of electrode [mm]
        contact_depth                               = tip_guide_tube + (tip_electrode-contact_position); % Contact depth from dura [mm]
        
        par.(['rec' num2str(c)]).date          = '20250313';        	% Date of recording
        par.(['rec' num2str(c)]).nChan         = 1:32;               	% [superficial deep]
        par.(['rec' num2str(c)]).array_depth   = contact_depth;      	% [superficial deep]
        
        %%% Add channel-wise information %%%
        for iChan = 1:32
            par.(['rec' num2str(c)]).(['ch' num2str(iChan)]).coord   = [-1.5 1.5 contact_depth(iChan)]; % Center of grid at level of dura [0x 0y 0z]mm
            par.(['rec' num2str(c)]).(['ch' num2str(iChan)]).field   = 'MT/MST';
            par.(['rec' num2str(c)]).(['ch' num2str(iChan)]).clus    = cluster{iChan};
        end
        c                                           = c+1;                  % Count sessions
        
    case 'Dunkin'
        
        
end