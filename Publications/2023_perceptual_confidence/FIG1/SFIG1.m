% Add relevant directories
addpath /Users/fschneider/Documents/MATLAB/CircStat2012a/
addpath /Users/fschneider/Documents/GitHub/Violinplot-Matlab
addpath /Users/fschneider/Documents/MATLAB/cbrewer/

close all
clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get subject-wise RT %%%
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
%%% Subject-wise RT %%%
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
    pl.Color            = [.5 .5 .5 .3];
    
    hold on
end

ax1                     = gca;
ax1.RLim                = [0 400];
ax1.FontSize            = lb_fs;

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1b', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1b', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNR %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f                      	= figure('units','centimeters','position',[0 0 12 7]); hold on
lw                      = 2;
col                     = [.4 .4 .4];
alph = .4;

for iSubj = 1:length(id)
    try
        load([local_pth id{iSubj} '/summary/' id{iSubj} '_SNRfit.mat'])
    catch
        continue
    end
    
    psy_func.model_snr(1) = .0001;
    psy_func.snr(1)       = .01;

    pl                      = plot((psy_func.model_snr),psy_func.model,'-','color',[col alph],'linewidth',lw);
    sc                      = scatter(psy_func.snr,psy_func.hir);
    sc.Marker               = '.';
    sc.MarkerEdgeColor      = col;
    sc.MarkerEdgeAlpha      = alph;
    sc.MarkerFaceAlpha      = alph;
    sc.SizeData             = 100;
end

for i = 1:length(psy_func.snr)
    lab{i} = num2str(round(psy_func.snr,2)*100);
end

ax1 = gca;
ax1.XScale                  = 'log';
ax1.YLim                    = [0 1];
ax1.XLim                    = [.01 1];
ax1.FontSize                = lb_fs;
ax1.XTick                   = psy_func.snr;
ax1.XTickLabel              = lab;
ax1.YTick                   = [0 .25 .5 .75 1];
ax1.XLabel.String           = 'Coherence [%]';
ax1.YLabel.String           = '% correct';
ax1.XTickLabelRotation      = 30
ln                          = line(ax1.XLim,[.25 .25],'Color', 'k', 'LineStyle', '--');

for i = 1:length(ax1.XTickLabel)
    if i == 1
        ax1.XTickLabel{i}   = 0;
    else
        ax1.XTickLabel{i} 	= round(psy_func.snr(i),2);
    end
end

grid on
grid minor 
grid minor

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1d', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1d', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Subject-wise score over time %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% [trg_score, exp, dte]   = getFinalScore(local_pth, id);
% dest_dir = '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/';
% save([dest_dir 'subject_scores.mat'], 'trg_score', 'exp', 'dte', '-v7.3');
load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/subject_scores.mat')

f                      	= figure('units','centimeters','position',[0 0 7 7]); hold on
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

    P               	= polyfit(1:length(dat),dat,1);
    slpe(iSubj)       	= P(1);
    yfit            	= P(1)*(1:length(dat))+P(2);  % P(1) is the slope and P(2) is the intercept
    pf              	= plot(1:length(dat),yfit,'k');
end

ax3                     = gca;
ax3.XLabel.String       = 'Exp. block';
ax3.YLabel.String       = 'Final cum. score';
ax3.FontSize            = lb_fs;
ax3.XLim                = [1 19];
ax3.YLim                = [50 275];

print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1c', '-r500', '-dpng');
print(f, '/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/FIG1/SFIG1/SFIG1c', '-r500', '-dsvg', '-painters');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

load('/Users/fschneider/Documents/GitHub/CPR/Publications/2023_perceptual_confidence/var_plot/solo_correlation.mat')

f                      	= figure('units','centimeters','position',[0 0 15 15]); hold on

snr = unique(solo_cr{1}.coh);
cnt = 0;
for iSubj = 1:length(solo_cr)
    cnt = cnt+1;
    avg_xc = [];
   for iCoh = 1:length(snr) 
       
       if iSubj == 1
           lab{iCoh} = num2str(round(snr(iCoh),2));
           
       end
       cIdx = (solo_cr{iSubj}.coh == snr(iCoh))';
       avg_xc(iCoh,:) = median(solo_cr{iSubj}.sxc(cIdx,151:end));
   end
   
   ax = subplot(5,8,cnt);
   imagesc(avg_xc ./ max(max(avg_xc)))
   
   ax.YTick = [1:7];
   ax.XTick = [75 150];
   ax.XTickLabels = cellfun(@(x) round(str2num(x).*8.333),ax.XTickLabels,'UniformOutput',false);
   ax.YTickLabels =lab;
   ax.FontSize = lb_fs;
   colormap(gray(256));
   
   if mod(iSubj,8) ~= 1
       ax.YAxis.Visible = 'off';
   end
   
   if iSubj < 33
       ax.XAxis.Visible = 'off';
   end
end

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
            c                       = c+1;
            tmp_tbl              	= load(mat_files(iFile).name);
            trg_score{iSubj,c}    	= getTargetScore(tmp_tbl.t);
            exp{iSubj,c}           	= tmp_tbl.t.exp(1);
            dte{iSubj,c}         	= tmp_tbl.t.date(1);
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
        idx                             = cellfun(@(x) strcmp(x,sbj_lst{iSubj}), subj);
        
        if sum(idx) == 0
            continue
        else
            files                      	= dnames(idx);
            b1                          = load(files{1});
            dcnt                       	= dcnt+1;
            trg_score{iSubj,dcnt}    	= getTargetScore(b1.t);
            exp{iSubj,dcnt}          	= b1.t.exp(1);
            dte{iSubj,dcnt}          	= b1.t.date(1);
            
            b2                          = load(files{2});
            dcnt                       	= dcnt+1;
            trg_score{iSubj,dcnt}    	= getTargetScore(b2.t);
            exp{iSubj,dcnt}         	= b2.t.exp(1);
            dte{iSubj,dcnt}          	= b2.t.date(1);
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
