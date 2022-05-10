% Target score
s_cnt                	= s_cnt+1;
n                       = sum(cellfun(@numel,t.trg_ts(logical(t.trg_shown))));
mat                     = t.trg_score;
t_cnt                	= 0;
score                   = [];

for i = 1:size(mat,1)
    if t.trg_shown(i) == 0
        continue
    end
    
    for j = 1:size(mat{i})
        if ~isnan(mat{i}(j))
            t_cnt           = t_cnt+1;
            score(t_cnt)    = mat{i}(j);
        end
    end
end

df                      = diff(score);
df(df < 0)              = [];
total                   = cumsum(df);
perc(s_cnt)         	= total(end)/n;