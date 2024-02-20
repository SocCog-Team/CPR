function out = CPR_target_reward(t)

out                 = [];

for i = 1:size(t.trg_ts,1)
    for iTrg = 1:length(t.trg_ts{i})
        
        [~,idx]     = min(abs(t.frme_ts{i} - t.trg_ts{i}(iTrg)));
        str         = t.js_str{i}(idx);
        acc         = abs(1 - abs(t.rdp_dir(i) - t.js_dir{i}(idx)) / 180);
  
        trg(1)      = str * acc;
        trg(2)      = t.ss_coh(i);
        trg(3)      = t.trg_hit{i}(iTrg);
        out         = [out; trg];
        
    end
end

figure
idx                 = logical(out(:,3));
vs                  = violinplot(out(idx,1),out(idx,2));

end