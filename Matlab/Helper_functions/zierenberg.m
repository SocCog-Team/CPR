%% Writing cycle-wise stimulus- and joystick direction

cd /Users/fschneider/Desktop/Bastian/summary_folder

mat_lst = dir('*.mat');
dest_dir = '/Users/fschneider/Desktop/zierenberg/';

for iFile = 1:length(mat_lst)
    
    if contains(mat_lst(iFile).name, 'CPRsolo')
        clear out
        load([mat_lst(iFile).folder '/' mat_lst(iFile).name])
                    
        csvwrite([dest_dir out.info.exp '_' out.info.subject_id '_' out.info.date '_' out.info.setup '_' out.info.block '_rdp_coh.csv'],out.coherence)
        csvwrite([dest_dir out.info.exp '_' out.info.subject_id '_' out.info.date '_' out.info.setup '_' out.info.block '_rdp_coh_num.csv'],out.block_num)
        
        for iCyc = 1:length(out.raw.rdp_dir) 
            str = [dest_dir out.info.exp '_' out.info.subject_id '_' out.info.date '_' out.info.setup '_' out.info.block '_cycle' num2str(iCyc)];
            csvwrite([str '_rdp_dir.csv'],out.raw.rdp_dir{iCyc}')
            csvwrite([str '_js_dir.csv'],out.raw.js_dir{iCyc}')
            csvwrite([str '_js_tlt.csv'],out.raw.js_ecc{iCyc}')           
        end      
    end
end