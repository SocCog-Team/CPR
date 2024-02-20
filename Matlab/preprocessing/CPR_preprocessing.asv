%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CPR DATA PREPROCESSING %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all

% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/preprocessing
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

% Copy subject summary from server
source_pth              = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/';
local_pth               = '/Volumes/T7_Shield/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';
[status,msg,msgID]      = copyfile([source_pth fname], local_pth);

% Import subject summary table
x                       = readtable([local_pth fname]);
sbj_lst                 = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

% Set params
skip_processed_files    = false;
copy_if_present         = false;

% Preprocess data
for iSubj = 29:length(sbj_lst)
    
    disp(['Processing subject: ' sbj_lst{iSubj}])
    
    % Create local directories
    data_pth                = [local_pth sbj_lst{iSubj} '/raw/'];
    summ_pth                = [local_pth sbj_lst{iSubj} '/summary/'];
    [status,msg,msgID]      = mkdir(data_pth);
    [status,msg,msgID]      = mkdir(summ_pth);
            
    % Change into data source directory
    cd([source_pth sbj_lst{iSubj} '/raw/'])

    % Extract .mwk2 file names from source directory
    mwk2_files              = dir('*.mwk2');
    
    for iFile = 1:length(mwk2_files)
        
        % Skip irrelevant files
        if contains(mwk2_files(iFile).name, 'CPRtrain') || ...
                (contains(mwk2_files(iFile).name, 'RT') && length(mwk2_files(iFile).name) > 24)
            continue
        end
        
        % Copy file to local directory
        if copy_if_present || ~isfile([data_pth mwk2_files(iFile).name])
            [status,msg,msgID]	= copyfile([source_pth sbj_lst{iSubj} '/raw/' mwk2_files(iFile).name], data_pth);
        end
        
        % Skip files if already processed
        if skip_processed_files
            fid = strsplit(mwk2_files(iFile).name,'_');
            
            if isfile([summ_pth strjoin(fid(2:3),'_') '.mat'])
                continue
            elseif isfile([summ_pth strjoin(fid(2:3),'_') 'fit.mat'])
                continue
            elseif isfile([summ_pth strjoin(fid(1:4),'_') '_tbl.mat'])
                continue
            end
        end
        
%         if contains(mwk2_files(iFile).name,'RT')
%             % Extract reaction time profile
%             [rt.dat,rt.trg_dir] = RT_polar(mwk2_files(iFile).name, data_pth, summ_pth);
%             % Save variables
%             save([summ_pth sbj_lst{iSubj} '_RT.mat'], 'rt', '-v7.3')
%         end
%         
%         if contains(mwk2_files(iFile).name,'SNR')
%             % Fit psychometric curves
%             [psy_func]   	= CPR_psychometric_function(data_pth, mwk2_files(iFile).name, summ_pth);
%             % Save variables
%             save([summ_pth sbj_lst{iSubj} '_SNRfit.mat'], 'psy_func', '-v7.3')
%         end
        
        if contains(mwk2_files(iFile).name,'CPRsolo') && length(mwk2_files(iFile).name) >= 44
            % Import CPRsolo data
            [t,d,s]         = CPR_psych_import(data_pth, mwk2_files(iFile).name);
        end
         
        if contains(mwk2_files(iFile).name,'CPRagent') && length(mwk2_files(iFile).name) >= 44
            % Import CPRagent data
            [t,d,s]         = CPR_psych_import(data_pth, mwk2_files(iFile).name);
        end
    end
    close all
end

%% Import dyadic data

xd                      = readtable([local_pth fname],'Sheet','Dyads');

for iDyad = 19:71
    % Specify directories
    data_pth            = [local_pth 'Dyad' num2str(iDyad) '/raw/'];
    summ_pth            = [local_pth 'Dyad' num2str(iDyad) '/summary/'];
    [status,msg,msgID] 	= mkdir(data_pth);
    [status,msg,msgID] 	= mkdir(summ_pth);
    
    % Extract .mwk2 file names from source directory
    cd([source_pth 'Dyad' num2str(iDyad) '/raw/'])
    mwk2_files              = dir('*.mwk2');

    for iFile = 1:length(mwk2_files)
        % Copy file to local directory
        if copy_if_present || ~isfile([data_pth mwk2_files(iFile).name])
            [status,msg,msgID]	= copyfile([source_pth 'Dyad' num2str(iDyad) '/raw/' mwk2_files(iFile).name], data_pth);
        end 
        
        if contains(mwk2_files(iFile).name,'CPRdyadic')
            % Import CPRdyadic data
            [t,d,s]         = CPR_psych_import(data_pth, mwk2_files(iFile).name);
        end
    end
end
