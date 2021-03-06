clear all 
close all

% Add relevant directories
addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

% Import subject summary table
pth                     = '/Users/fschneider/Documents/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';
x                       = readtable([pth fname]);
sbj_lst                 = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];
dte_format              = 'yyyymmdd';
clear_dir_flag          = false;

% Clear all summary directories
if clear_dir_flag == true
    for iSubj = 1:length(sbj_lst)
        summ_pth            = [pth sbj_lst{iSubj} '/summary/'];
        
        if isfolder(summ_pth)
            rmdir([pth sbj_lst{iSubj} '/summary/'],'s')
            mkdir([pth sbj_lst{iSubj} '/summary/'])
        end
    end
end

% Preprocess data
for iSubj = 1:length(sbj_lst)
        
    clear fname_RT fname_SNR fname_SOLO1 fname_SOLO2 fname_AGNT1 fname_AGNT2 fname_SOLO1r fname_SOLO2r
    disp(['Processing subject: ' sbj_lst{iSubj}])
    
    % Specify directories
    data_pth                = [pth sbj_lst{iSubj} '/raw/'];
    summ_pth                = [pth sbj_lst{iSubj} '/summary/'];
    cd(data_pth)
    
    % Pre-tests
    if ~isnat(x.DateTraining(iSubj))
        fname_RT            = [datestr(x.DateTraining(iSubj),dte_format) '_' sbj_lst{iSubj} '_RT_fxs.mwk2'];
        fname_SNR           = [datestr(x.DateTraining(iSubj),dte_format) '_' sbj_lst{iSubj} '_SNR_fxs.mwk2'];
            
        % Extract reaction time profile
        [rt, trg_dir]       = RT_polar(fname_RT, data_pth, summ_pth);
            
        % Fit psychometric curves
        CPR_psychometric_function(data_pth, fname_SNR, summ_pth)
    end
    
    % Import CPRsolo data
    if ~isnat(x.DateSolo1(iSubj))
        try
            [t1,d1,s1]   	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo1(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block1_psycho4_fxs.mwk2']);
            [t2,d2,s2]    	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo1(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block2_psycho3_fxs.mwk2']);
        catch
            [t1,d1,s1]    	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo1(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block1_psycho3_fxs.mwk2']);
            [t2,d2,s2]     	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo1(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block2_psycho4_fxs.mwk2']);
        end
    end
    
%     % Import CPRagent data
%     if ~isnat(x.DateAgent(iSubj))
%         try
%             [t1,d1,s1]   	= CPR_psych_import([sbj_lst{iSubj} 'Agnt'], data_pth, [datestr(x.DateAgent(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRagent_block1_psycho4_fxs.mwk2']);
%             [t2,d2,s2]   	= CPR_psych_import([sbj_lst{iSubj} 'Agnt'], data_pth, [datestr(x.DateAgent(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRagent_block2_psycho3_fxs.mwk2']);
%         catch
%             [t1,d1,s1]    	= CPR_psych_import([sbj_lst{iSubj} 'Agnt'], data_pth, [datestr(x.DateAgent(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRagent_block1_psycho3_fxs.mwk2']);
%             [t2,d2,s2]     	= CPR_psych_import([sbj_lst{iSubj} 'Agnt'], data_pth, [datestr(x.DateAgent(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRagent_block2_psycho4_fxs.mwk2']);
%         end
%     end
    
    % Import CPRsolo_repeated data
    if ~isnat(x.DateSolo2(iSubj))
        try
            [t1,d1,s1]      = CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo2(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block1_psycho4_fxs.mwk2']);
            [t2,d2,s2]   	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo2(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block2_psycho3_fxs.mwk2']);
        catch
            [t1,d1,s1]    	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo2(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block1_psycho3_fxs.mwk2']);
            [t2,d2,s2]   	= CPR_psych_import(sbj_lst{iSubj}, data_pth, [datestr(x.DateSolo2(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block2_psycho4_fxs.mwk2']);
        end
    end
 
%     % CPR-4AFC psychometric curve comparison
%     if exist(fname_SNR,'file') == 2 && exist([summ_pth DATA_TABLE])
%     end
end

%% Import dyadic data

% xd                      = readtable([pth fname],'Sheet','Dyads');
%     
% for iDyad = 14:17%xd.Dyad'
%     % Specify directories
%     data_pth            = [pth 'Dyad' num2str(iDyad) '/raw/'];
%     summ_pth            = [pth 'Dyad' num2str(iDyad) '/summary/'];
%     cd(data_pth)
%     
%     % Wipe directory
%     if isfolder(summ_pth)
%         rmdir(summ_pth,'s')
%         mkdir(summ_pth)
%     else
%         mkdir(summ_pth)
%     end
%     
%     % Specify file names
%     fname_block1        = [datestr(xd.Date(iDyad),dte_format) '_' xd.Block1_PSY4{iDyad} xd.Block1_PSY3{iDyad} '_CPRdyadic_block1_fxs.mwk2'];
%     fname_block2        = [datestr(xd.Date(iDyad),dte_format) '_' xd.Block2_PSY4{iDyad} xd.Block2_PSY3{iDyad} '_CPRdyadic_block2_fxs.mwk2'];
%     
%     % Import data: Subject naming convention: Primary (PSY4) - Secondary (PSY3)
%     [t_b1,d_b1,s_b1]    = CPR_psych_import([xd.Block1_PSY4{iDyad} xd.Block1_PSY3{iDyad}], data_pth, fname_block1);
%     [t_b2,d_b2,s_b2]   	= CPR_psych_import([xd.Block2_PSY4{iDyad} xd.Block2_PSY3{iDyad}], data_pth, fname_block2);
% 
% %     % Merge data table for each subject
% %     subj1_merged        = [t_b1{1}; t_b2{1}];
% %     subj2_merged        = [t_b1{2}; t_b2{1}];
%     
% end
