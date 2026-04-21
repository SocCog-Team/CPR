local_pth               = '/Volumes/T7_Shield/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';

x                       = readtable([local_pth fname]);
id                      = cellfun(@lower, x.Abbreviation, 'UniformOutput', false);

for iSubj = 1:length(id)
    cd([local_pth id{iSubj} '/summary/'])
    load([id{iSubj} '_RT.mat'])
    
    RT_median(iSubj)  	= median(rt.dat);
    RT_iqr(iSubj)     	= iqr(rt.dat);
end

info_tbl                = [table(id,'VariableNames',{'ID'})  x(:,3:5) x(:,7) table(RT_median',RT_iqr','VariableNames',{'RT_median', 'RT_IQR'})];
dest_dir                = '/Users/fschneider/ownCloud/Wibral_lab/';

writetable(info_tbl,[dest_dir 'subject_info.csv']);
writetable(info_tbl,[dest_dir 'subject_info.xlsx']);