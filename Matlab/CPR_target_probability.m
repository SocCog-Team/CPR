close all

f               = figure; hold on
p_targ          = [.01 .005 .001];
n_smpl          = 1000;
pow             = 1:n_smpl;
cl              = [.3 .3 .3];
p_interest      = .5;

sl              = line([125 125],[0 1], 'Color',[.75 .75 .75], 'LineStyle', ':', 'LineWidth',2);
st              = text(135,.99,'Min state duration', 'Color',[.75 .75 .75]);

for i = 1:size(p_targ,2)
        
    p_NO_target   	= ((1-p_targ(i)) .^ (pow-1));
    p_target      	= p_NO_target * p_targ(i); % P(X=k)= (1-p)^(k-1) * p
    sum_p_target   	= cumsum(p_target);
    [val, idx]   	= min(abs(sum_p_target-p_interest));
    
    p(i)            = plot(sum_p_target, 'LineWidth',2);
    lv              = line([idx idx],[0 sum_p_target(idx)],'Color', cl, 'LineWidth',1);
    lh              = line([0 idx],[sum_p_target(idx) sum_p_target(idx)],'Color',cl, 'LineWidth',1);
    t               = text(idx+5,p_interest-.03,num2str(idx*10));
end

ax                  = gca;
ax.FontSize         = 20;
ax.XTickLabel       = cellfun(@str2num, ax.XTickLabel) / 100;
ax.YLabel.String    = 'Target Probability';
ax.XLabel.String    = 'Time [s]';
ax.Title.String     = 'Probability of Target Appearance';
p_targ              = p_targ * 100;
lg                  = legend([p(1) p(2) p(3)],...
    {[num2str(p_targ(1)) '%'] ...
    [num2str(p_targ(2)) '%'] ...
    [num2str(p_targ(3)) '%']},...
    'Location','southeast');
lg.Title.String     = {'Sample-wise'; 'target probability'};
lg.FontSize         = 12;