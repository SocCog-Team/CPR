spth = '/Volumes/T7_Shield/CPR_psychophysics/';
dpth = '/Volumes/DPZ/KognitiveNeurowissenschaften/CNL/DATA/fxs/CPR_psychophysics/';

cd(spth);
folders = dir;

% Import subject summary table
fname                   = 'Subjects_summary.xlsx';
x                       = readtable([dpth fname]);
sbj_lst                 = x.Abbreviation;
sbj_lst(cellfun(@isempty,sbj_lst)) = [];

for i = 1:length(folders)
    
    % Only incude folders with three or six characters
    if (length(folders(i).name) == 3 || length(folders(i).name) == 6) && isfolder([spth folders(i).name '/summary/'])
        % Check if part of subject list or dyad
        if sum(cellfun(@(x) contains(x,folders(i).name),sbj_lst)) == 1 || contains(folders(i).name,'Dyad')
        else
            continue
        end
        
        disp(['Copy summary: ' folders(i).name])
        new_dir = [dpth folders(i).name  '/summary/'];
        copyfile([spth folders(i).name '/summary/'], new_dir)
    end
    disp('Done!')
    
%     disp(['Copy H5: ' folders(i).name])
%     if length(folders(i).name) == 3 || length(folders(i).name) == 6
%         new_dir = [dpth folders(i).name  '/raw/'];
%         copyfile([spth folders(i).name '/raw/'], new_dir)
%     end
%     disp('Done!')
end