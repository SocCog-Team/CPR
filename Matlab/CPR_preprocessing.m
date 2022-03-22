clear all 
close all

addpath /Users/fschneider/ownCloud/Shared/MWorks_MatLab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/WIP/
addpath /Users/fschneider/Documents/GitHub/CPR/Matlab/Helper_functions/
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

pth                     = '/Users/fschneider/Documents/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';
x                       = readtable([pth fname]);
sbj_lst                 = x.Abbreviation;
dte_format              = 'yyyymmdd';

for iSubj = 1:length(sbj_lst)
    summ_pth            = [pth sbj_lst{iSubj} '/summary/'];
    
    if isfolder(summ_pth)
        rmdir([pth sbj_lst{iSubj} '/summary/'],'s')
        mkdir([pth sbj_lst{iSubj} '/summary/'])
    end
end

for iSubj = 1:length(sbj_lst)
    
    % Specify directories
    data_pth            = [pth sbj_lst{iSubj} '/raw/'];
    summ_pth            = [pth sbj_lst{iSubj} '/summary/'];
    cd(data_pth)
    
    % Specify file names
    fname_RT            = [datestr(x.DateTraining(iSubj),dte_format) '_' sbj_lst{iSubj} '_RT_fxs.mwk2'];
    fname_SNR           = [datestr(x.DateTraining(iSubj),dte_format) '_' sbj_lst{iSubj} '_SNR_fxs.mwk2'];
    fname_SOLO1         = [datestr(x.DateSolo1(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block1_fxs.mwk2'];
    fname_SOLO2         = [datestr(x.DateSolo1(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block2_fxs.mwk2'];
    fname_AGNT1         = [datestr(x.DateAgent(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRagent_block1_fxs.mwk2'];
    fname_AGNT2       	= [datestr(x.DateAgent(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRagent_block2_fxs.mwk2'];
    
    if ~isnat(x.DateSolo2(iSubj))
        fname_SOLO1r	= [datestr(x.DateSolo2(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block1_fxs.mwk2'];
        fname_SOLO2r	= [datestr(x.DateSolo2(iSubj),dte_format) '_' sbj_lst{iSubj} '_CPRsolo_block2_fxs.mwk2'];
    end
    
    % Extract reaction time profile
    if exist(fname_RT,'file') == 2
        RT_polar(fname_RT, data_pth, summ_pth)
    end
    
    % Fit psychometric curves
    if exist(fname_SNR,'file') == 2
        CPR_psychometric_function(data_pth, fname_SNR, summ_pth)
    end
    
    % Import CPRsolo data
    if exist(fname_SOLO1,'file') == 2 && exist(fname_SOLO2,'file') == 2
        [t,d,s] 	= CPR_psych_import(sbj_lst{iSubj}, data_pth, {fname_SOLO1,fname_SOLO2});
    end
    
    % Import CPRsolo_repeated data
    if exist(fname_SOLO1r,'file') == 2 && exist(fname_SOLO2r,'file') == 2
        [t,d,s]     = CPR_psych_import(sbj_lst{iSubj}, data_pth, {fname_SOLO1r,fname_SOLO2r});
    end
    
    % Import CPRagent data
    if exist(fname_AGNT1,'file') == 2 && exist(fname_AGNT2,'file') == 2
        [t,d,s]     = CPR_psych_import(sbj_lst{iSubj}, data_pth, {fname_AGNT1,fname_AGNT2});
    end
    
%     % CPR-4AFC psychometric curve comparison
%     if exist(fname_SNR,'file') == 2 && exist([summ_pth DATA_TABLE])
%     end
end

%% Import dyadic data

xd                      = readtable([pth fname],'Sheet','Dyads');
    
for iDyad = xd.Dyad
    % Specify directories
    data_pth            = [pth 'Dyad' num2str(iDyad) '/raw/'];
    summ_pth            = [pth 'Dyad' num2str(iDyad) '/summary/'];
    cd(data_pth)
    
    % Specify file names
    fname_block1        = [datestr(xd.Date(iDyad),dte_format) '_' xd.Block1_PSY4{iDyad} xd.Block1_PSY3{iDyad} '_CPRdyadic_block1_fxs.mwk2'];
    fname_block2        = [datestr(xd.Date(iDyad),dte_format) '_' xd.Block2_PSY4{iDyad} xd.Block2_PSY3{iDyad} '_CPRdyadic_block1_fxs.mwk2'];
    
    % Import data: Subject naming convention: Primary (PSY4) - Secondary (PSY3)
    [t_b1,d_b1,s_b1]    = CPR_psych_import([xd.Block1_PSY4{iDyad} xd.Block1_PSY3{iDyad}], data_pth, fname_block1);
    [t_b2,d_b2,s_b2]   	= CPR_psych_import([xd.Block2_PSY4{iDyad} xd.Block2_PSY3{iDyad}], data_pth, fname_block2);

    % Merge data table for each subject
    subj1_merged        = [t_b1{1}; t_b2{1}];
    subj2_merged        = [t_b1{2}; t_b2{1}];
    
end
