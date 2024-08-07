% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Session ID %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')

pth             = '/Volumes/T7_Shield/CPR_psychophysics/';      % Local hard drive
cnt             = 0;
lb_fs           = 8;

for iSubj = 1:length(solo_cr)
    % Go to directory
    cd([pth solo_cr{iSubj}.id '/summary/'])
    
    % Check files
    mat_files 	= dir('*.mat');                         % Get all .mat files in directory
    
    % Extract date and id
    for iFile = 1:length(mat_files)
        if contains(mat_files(iFile).name,'CPR')
            tmp             = split(mat_files(iFile).name,'_');
            
            if strcmp(tmp{4}, 'block1') && ~strcmp(tmp{2}, 'agnt')
                cnt             = cnt+1;
                sub_id{cnt}     = lower(tmp{2});
                exp_id{cnt}     = tmp{3};
                exp_date{cnt}   = tmp{1};
            end
        end
    end
end

for iDyad = 19:71
    % Go to directory
    cd([pth '/Dyad' num2str(iDyad) '/summary/'])

    % Check files
    mat_files 	= dir('*.mat');                  
    
    % Extract date and id
    for iFile = 1:length(mat_files)
        if contains(mat_files(iFile).name,'CPR')
            tmp             = split(mat_files(iFile).name,'_');
            if strcmp(tmp{4}, 'block1')
                cnt             = cnt+1;
                sub_id{cnt}     = tmp{2};
                exp_id{cnt}     = tmp{3};
                exp_date{cnt}   = tmp{1};
            end
        end
    end
end

spool = unique(sub_id);

for iSubj = 1:length(spool)
    sIdx = cellfun(@(x) strcmp(x,spool{iSubj}),sub_id);
    
    dte = exp_date(sIdx);
    exp = exp_id(sIdx);
    
    [dte_sorted, idx] = sort(dte);
    exp_sorted = exp(idx);
    
    for iExp = 1:length(exp_sorted)
        if strcmp(exp_sorted{iExp}, 'CPRsolo')
            flag(iSubj,iExp) = 1;
        elseif strcmp(exp_sorted{iExp}, 'CPRdyadic')
            flag(iSubj,iExp) = 2;
        elseif strcmp(exp_sorted{iExp}, 'CPRagent')
            flag(iSubj,iExp) = 3;
        end
    end
end

[~,flag_sorted] = sort(sum(flag ~= 0,2));

f                   = figure('units','centimeters','position',[0 0 7 7]); hold on
% col                 = cbrewer('qual', 'Set3', 3, 'PCHIP');
col                 = cool(3);
im                  = imagesc(flag(flag_sorted,:));
ax                  = gca;
ax.XLabel.String  	= 'Session number';
ax.YLabel.String 	= 'Subject';
ax.FontSize       	= lb_fs;
ax.XLim          	= [.5 9.5];
ax.YLim            	= [.5 38.5];

colormap([0 0 0; col])

h(1)                = plot(NaN,NaN,'.','Color', col(1,:));
h(2)                = plot(NaN,NaN,'.','Color', col(2,:));
h(3)                = plot(NaN,NaN,'.','Color', col(3,:));
[lg,icons]          = legend(h,'Solo', 'Dyadic_Human', 'Dyadic_Computer', 'Interpreter', 'none');
lg.Location         = 'southeast';
lg.Box              = 'on';
icons               = findobj(icons,'Type','line');
icons               = findobj(icons,'Marker','none','-xor');
set(icons,'MarkerSize',30);

axis equal
axis tight

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1x', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1x', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Subject-wise score over time %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
local_pth = '/Volumes/T7_Shield/CPR_psychophysics/';
x = readtable([pth 'Subjects_summary.xlsx']); % Spreadsheet
sbj_lst = x.Abbreviation; % Subject ID list
sbj_lst(cellfun(@isempty,sbj_lst)) = [];
[trg_score, exp, dte]   = getFinalScore(local_pth, sbj_lst);
dest_dir = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/';
save([dest_dir 'subject_scores.mat'], 'trg_score', 'exp', 'dte', '-v7.3');
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/subject_scores.mat')

f                      	= figure('units','centimeters','position',[0 0 7 7]); hold on
alph                    = .3;

for iSub = 1:size(trg_score,1)
    % sort to date
    for i = 1:size(dte,2)
        if isempty(dte{iSub,i})
            continue
        end
        
        if dte{iSub,i} == "202220818"
          dte{iSub,i} =  "20220818";
        end
        
        t(i) = datetime(dte{iSub,i},'InputFormat','yyyyMMdd');
    end

    [v,idx]             = sort(t);
    dat              	= cell2mat(trg_score(iSub,idx));
    plt              	= plot(1:length(dat),dat,'Color',[.5 .5 .5 alph/2], 'LineWidth', 1); uistack(plt,'bottom')
    sc                  = scatter(1:length(dat),dat,'Marker','x','MarkerFaceColor',[.5 .5 .5],'MarkerFaceAlpha',alph,'MarkerEdgeColor',[.5 .5 .5],'MarkerEdgeAlpha',alph); uistack(plt,'bottom');

    hold on
    P               	= polyfit(1:length(dat),dat,1);
    slpe(iSubj)       	= P(1);
    yfit            	= P(1)*(1:length(dat))+P(2);  % P(1) is the slope and P(2) is the intercept
    pf              	= plot(1:length(dat),yfit,'k', 'LineWidth', 1);
    uistack(pf,'top') 

end

ax3                     = gca;
ax3.XLabel.String       = 'Exp. block';
ax3.YLabel.String       = 'Final cum. score';
ax3.FontSize            = lb_fs;
ax3.XLim                = [1 18];
ax3.YLim                = [50 275];

axis square
axis tight

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1c', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1c', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Response lag %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')

f                      	= figure('units','centimeters','position',[0 0 18 18]); hold on
dim                     = [.1 .175];
height                  = [0 linspace(.8,.07,5)];
colmn                   = linspace(.075, .8,8);
r = 1;

snr = unique(solo_cr{1}.coh);
cnt = 0;

axis off

for iSubj = 1:length(solo_cr)
    cnt = cnt+1;
    avg_xc = [];
   for iCoh = 1:length(snr) 
       
       if iSubj == 1
           lab{iCoh} = num2str(round(snr(iCoh),2)*100);
           
       end
       cIdx = (solo_cr{iSubj}.coh == snr(iCoh))';
       avg_xc(iCoh,:) = median(solo_cr{iSubj}.sxc(cIdx,151:end));
   end
   
   c = mod(iSubj,8);
   flag = 0;
   if c == 0
       c = 8;
       flag = 1;
   end
   
   if c == 1
       r = r+1;
   end
   
   ax	= axes('Position', [colmn(c) height(r) dim]); hold on
   imagesc(avg_xc ./ max(max(avg_xc)))
   
   ax.YTick = [1:7];
   ax.XTick = [75 150];
   ax.XTickLabels = cellfun(@(x) round(str2num(x).*8.333),ax.XTickLabels,'UniformOutput',false);
   ax.YTickLabels = lab;
   ax.YLabel.String = 'Coherence';
   ax.XLabel.String = 'Lag [ms]';
   ax.FontSize = lb_fs;
   ax.XTickLabelRotation = 30;
   colormap(gray(256));
   axis tight
   
   if mod(iSubj,8) ~= 1
       ax.YAxis.Visible = 'off';
   end
   
   if iSubj < 33
       ax.XAxis.Visible = 'off';
   end
end

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1l', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1l', '-r500', '-dsvg', '-painters');


f               = figure('units','centimeters','position',[0 0 3 3]);
cb              = colorbar;
cb.Ticks        = [0 1];
cb.Position     = [.4 .1 .2 .8];
cb.Label.FontSize = 9;
cb.Label.String = {'Cross-Correlation'  '[norm]'};
cl              = linspace(0,.99,256)';

axis off
colormap(gray(256))

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1_col_bar', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Extract subject-wise RT data %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
local_pth               = '/Volumes/T7_Shield/CPR_psychophysics/';
fname                   = 'Subjects_summary.xlsx';
x                       = readtable([local_pth fname]);
id                      = cellfun(@lower, x.Abbreviation, 'UniformOutput', false);

% Get Median
for iSubj = 1:length(id)
    cd([local_pth id{iSubj} '/summary/'])
    load([id{iSubj} '_RT.mat'])
    mRT(iSubj)          = median(rt.dat);
end

% Sort
[~,idx]                 = sort(mRT);

% Get full data
c = 0;
for iSubj = idx
    cd([local_pth id{iSubj} '/summary/'])
    load([id{iSubj} '_RT.mat'])
    
    c = c+1;
    RT{c}           = rt.dat';
    dir{c}          = rt.trg_dir';
    num{c}          = repmat(c,[length(RT{c}) 1]);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Avg Subject-wise RT %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f                      	= figure('units','centimeters','position',[0 0 17 7]);
lb_fs                 	= 8;
lg_fs                 	= 8;
vl                      = violinplot(cell2mat(RT'),cell2mat(num'));

for i = 1:length(vl)
    vl(i).ViolinColor{1}        = [.5 .5 .5];
    vl(i).ViolinAlpha{1}       	= .1;
    vl(i).ViolinPlot.EdgeAlpha  = 0;
    vl(i).BoxColor              = [0 0 0];
    vl(i).BoxWidth              = .075;
    vl(i).BoxPlot.FaceAlpha     = 1;
    vl(i).BoxPlot.EdgeAlpha     = 1;
end

ax0                     = gca;
ax0.YLim                = [150 800];
ax0.XTick               = [];
ax0.YLabel.String       = 'Reaction time [ms]';
ax0.XLabel.String       = 'Subjects';
ax0.FontSize            = lb_fs;

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1a', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1a', '-r500', '-dsvg', '-painters');

f                      	= figure('units','centimeters','position',[0 0 5 5]);
vec                     = 0:30:360;

for iSubj = 1:length(RT)
    for iDir = 1:size(vec,2)-1
        didx            = dir{iSubj} > vec(iDir) & dir{iSubj} <= vec(iDir+1);
        mrt(iDir)       = median(RT{iSubj}(didx));
        tdir(iDir)      = mean([vec(iDir) vec(iDir+1)]);
    end
    
    pl                  = polarplot(deg2rad([tdir tdir(1)]), [mrt mrt(1)]);
    pl.LineWidth        = 2;
    pl.Color            = [0 0 0 .3];
    
    hold on
end

ax1                     = gca;
ax1.RLim                = [0 400];
ax1.FontSize            = lb_fs;

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1b', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1b', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [trg_score, exp, dte] = getFinalScore(pth, sbj_lst)

for iSubj = 1:length(sbj_lst)
    
    disp(['Score extraction subject: ' sbj_lst{iSubj}])
    data_pth                      	= [pth sbj_lst{iSubj} '/summary/'];
    cd(data_pth)
    
    mat_files                      	= dir('*.mat');
    c                               = 0;
    for iFile = 1:length(mat_files)
        clear tmp_tbl
        
        if contains(mat_files(iFile).name,'CPR') && ~contains(mat_files(iFile).name,'agnt')
            tmp_tbl              	= load(mat_files(iFile).name);            
            nT                      = sum(cellfun(@numel,tmp_tbl.t.trg_ts(logical(tmp_tbl.t.trg_shown))));
            if nT >= 600
                c                  	= c+1;
                trg_score{iSubj,c} 	= getTargetScore(tmp_tbl.t);
                exp{iSubj,c}       	= tmp_tbl.t.exp(1);
                dte{iSubj,c}      	= tmp_tbl.t.date(1);
            end
        end
    end
    
    dcnt = c;
    
    for iDyad = 19:71
        cd(['/Volumes/T7_Shield/CPR_psychophysics/Dyad' num2str(iDyad) '/summary/'])
        
        dfiles                          = dir('*.mat');
        
        if isempty(dfiles)
            continue
        end
        
        dnames                          = {dfiles.name};
        fsplit                          = cellfun(@(x) strsplit(x,'_'),dnames,'UniformOutput',false);
        subj                            = cellfun(@(x) x(2), fsplit);
        idx                             = cellfun(@(x) strcmp(x,lower(sbj_lst{iSubj})), subj);
        
        if sum(idx) == 0
            continue
        else
            files                      	= dnames(idx);
            b1                          = load(files{1});
            nT1                      	= sum(cellfun(@numel,b1.t.trg_ts(logical(b1.t.trg_shown))));
            if nT1 >= 600
                dcnt                  	= dcnt+1;
                trg_score{iSubj,dcnt} 	= getTargetScore(b1.t);
                exp{iSubj,dcnt}        	= b1.t.exp(1);
                dte{iSubj,dcnt}        	= b1.t.date(1);
            end
            
            b2                          = load(files{2});
            nT2                      	= sum(cellfun(@numel,b2.t.trg_ts(logical(b2.t.trg_shown))));
            if nT2 >= 600
                dcnt                   	= dcnt+1;
                trg_score{iSubj,dcnt}  	= getTargetScore(b2.t);
                exp{iSubj,dcnt}        	= b2.t.exp(1);
                dte{iSubj,dcnt}        	= b2.t.date(1);
            end
        end
    end
end
end

function [trg_score] = getTargetScore(in)
c                           = 0;

% Target response params
for iState = 1:size(in.trg_ts,1)
    if in.trg_shown(iState) == true
        for iTrg = 1:length(in.trg_ts{iState})
            c                   = c +1;
            score_cum(c)        = in.trg_score{iState}(iTrg);
            score_coh(c)        = in.rdp_coh(iState);
            
            if in.trg_score{iState}(iTrg) == 0
                score_hi(c)  	= false;
            else
                score_hi(c)  	= in.trg_hit{iState}(iTrg);
            end
        end
    end
end
trg_score                   = sum(score_cum(logical(score_hi)));
end
