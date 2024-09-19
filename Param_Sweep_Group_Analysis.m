%% Param Sweep Group Analysis
addpath(genpath('/Users/nkb/Documents/NCAN/code/BLAES_stimSweep'))
addpath(genpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools'))
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
bci2ktools(BCI2KPath);
%% Load Data and Add necessary field for group analysis
regionTable = readtable(fullfile("/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param","BrainRegionDictionary.xlsx"));
boxpath = ('/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param');
groupPath = fullfile(boxpath,'group');
subject_info = readtable(fullfile(boxpath,'Subject_Locations.xlsx'));
labs = {'stim_side', 'pair', 'anode', 'cathode', 'anode_region', 'cathode_region'};
files = dir(fullfile(groupPath,'*_cohens.mat'));
subjects = cell(length(files),1);
brains = struct;
for i=1:length(files)
    subject = split(files(i).name,"_cohens");
    subject = subject{1};
    if contains(subject,'UIC')
        UtahFlag = 1;
    else
        UtahFlag = 0;
    end
    subjects{i} =subject;
    
    stimDat = subject_info(ismember(subject_info.Subject,subject),:);
    if i<2
       
        
        temp = load(fullfile(groupPath,files(1).name));
        temp = temp.outstruct;
        tt = repmat({subject},length(temp),1);
        if UtahFlag
            localChans = {temp.channel};
        else
            localChans = remakeChannels({temp.channel},{temp.shank});
        end
        [temp(:).subject] = deal(tt{:});
        [temp(:).rec_chan] = deal(localChans{:});
        chanIdx = [temp.loc] == stimDat.loc1;
        stimInfo = cell(length(temp),6); % {side, pair, anode, cathode, anode region, cathode region}
        [pair1_chan, pair1_polarity] = parse_stimpair(stimDat.Pair1);
        [pair2_chan, pair2_polarity] = parse_stimpair(stimDat.Pair2);
        allChans = {pair1_chan{:}; pair2_chan{:}}';
        allChans = {allChans{:}};
        allRegion = cell(size(allChans));
        regions = {temp.region};

        for j=1:length(allChans)
            logIdx = ismember({temp.rec_chan},allChans{j});
            region = unique(regions(logIdx));
            allRegion{j} = region{1};

        end
        insertDat = repmat({stimDat.side1{1}, stimDat.Pair1{1}, pair1_chan{1}, pair1_chan{2}, allRegion{1},allRegion{2}},sum(chanIdx),1);
        stimInfo(chanIdx,:) = insertDat;
        insertDat = repmat({stimDat.side2{1}, stimDat.Pair2{1}, pair2_chan{1}, pair2_chan{2}, allRegion{3},allRegion{4}},sum(~chanIdx),1);
        stimInfo(~chanIdx,:) = insertDat;
        for j=1:length(labs)
            [temp(:).(labs{j})] = deal(stimInfo{:,j});
        end
        corrdat = temp;
        
    else
        temp = load(fullfile(groupPath,files(i).name));
        temp = temp.outstruct;
        tt = repmat({subject},length(temp),1);
        [temp(:).subject] = deal(tt{:});
        if UtahFlag
            localChans = {temp.channel};
        else
            localChans = remakeChannels({temp.channel},{temp.shank});
        end

        [temp(:).rec_chan] = deal(localChans{:});
        chanIdx = [temp.loc] == stimDat.loc1;
        stimInfo = cell(length(temp),6); % {side, pair, anode, cathode, anode region, cathode region}
        [pair1_chan, pair1_polarity] = parse_stimpair(stimDat.Pair1);
        [pair2_chan, pair2_polarity] = parse_stimpair(stimDat.Pair2);
        allChans = {pair1_chan{:}; pair2_chan{:}}';
        allChans = {allChans{:}};
        allRegion = cell(size(allChans));
        regions = {temp.region};

        for j=1:length(allChans)
            logIdx = ismember({temp.rec_chan},allChans{j});
            region = unique(regions(logIdx));
            allRegion{j} = region{1};

        end
        insertDat = repmat({stimDat.side1{1}, stimDat.Pair1{1}, pair1_chan{1}, pair1_chan{2}, allRegion{1},allRegion{2}},sum(chanIdx),1);
        stimInfo(chanIdx,:) = insertDat;
        insertDat = repmat({stimDat.side2{1}, stimDat.Pair2{1}, pair2_chan{1}, pair2_chan{2}, allRegion{3},allRegion{4}},sum(~chanIdx),1);
        stimInfo(~chanIdx,:) = insertDat;
        for j=1:length(labs)
            [temp(:).(labs{j})] = deal(stimInfo{:,j});
        end
        corrdat = [corrdat temp];
    end
    tempBrain = load(fullfile(boxpath,subject,sprintf('%s_MNI_new.mat',subject)));
    tempBrain = verifyElectrodes(tempBrain,unique({temp.channel}));
    brains(i).brain = tempBrain;
    brains(i).subject = subject;
end
clear temp tt
recSide = {corrdat.shank}; stimSide = {corrdat.stim_side};
ipsi = checkIpsiShank(recSide,stimSide)';
[corrdat(:).rec_side] = deal(ipsi{:});
stimregion = identifyStimRegion({corrdat.anode_region},{corrdat.cathode_region});
[corrdat(:).stimregion] = deal(stimregion{:});
all_rois = {corrdat.region};
pattern = '^ctx-(rh|lh)-';
all_rois = cellfun(@(x) regexprep(x, pattern, ''), all_rois, 'UniformOutput', false);
pattern = '^(Right|Left)-';
all_rois = cellfun(@(x) regexprep(x, pattern, ''), all_rois, 'UniformOutput', false);
[corrdat(:).region_class] = deal(all_rois{:});
simplified_regions = poolBrainRegions(all_rois,regionTable,boxpath);
[corrdat(:).pooled_region] = deal(simplified_regions{:});
gammaDist = squeeze(mean(cat(3,corrdat.gamma_change),1));
thetaDist = squeeze(mean(cat(3,corrdat.theta_change),1));
gammaMeans = squeeze(median((mean(cat(3,corrdat.gamma_change),2)),1));
thetaMeans = squeeze(median((mean(cat(3,corrdat.theta_change),2)),1));
[x,y]=histcounts(categorical({corrdat.subject}));
for j=1:length(corrdat)
    sub = corrdat(j).subject;
    l = ismember(y,sub);
    numCompares = x(l);
    theta_dist = thetaDist(:,j);
    theta_p = signrank(theta_dist-1);
    gamma_dist = gammaDist(:,j);
    gamma_p = signrank(gamma_dist-1);
    corrdat(j).gamma_m = gammaMeans(j);
    corrdat(j).gamma_p = gamma_p * numCompares;
    corrdat(j).theta_m = thetaMeans(j);
    corrdat(j).theta_p = theta_p * numCompares;
    corrdat(j).post_p = corrdat(j).post_p * numCompares;
    corrdat(j).stim_p = corrdat(j).stim_p * numCompares;

end
%% Plot All Trajectories
trajFig = figure;
ax = gca;
ax = plotAllCoverageMNI(ax,brains,0);

%% stim vs post-stim coherence
stim_d = [corrdat.stim_d];
post_d = [corrdat.post_d];
cmap = [[1 1 1]; [1 0 0]; [0 1 0]; [0 0 1]];
alpha_sig = 0.05;
agg = [[corrdat.stim_p];[corrdat.post_p]]';
both_sig = all(agg < alpha_sig, 2);
stim_sig = agg(:,1) < alpha_sig ;
gamma_sig = [corrdat.gamma_p] < alpha_sig;
theta_sig = [corrdat.theta_p] < alpha_sig;
stim_nonsig = agg(:,2) >= alpha_sig;
stim_res = stim_sig & stim_nonsig;
post_sig = agg(:,2) < alpha_sig;
post_nonsig = agg(:,1) >= alpha_sig;
post_res = post_sig & post_nonsig;
sig_map = both_sig + 2*stim_res + 3*post_res+1;
colors = cmap(sig_map,:);
% inverse_map = (both_sig + stim_res + post_res -1 ) *-1;
inverse_map = logical(sig_map -1);

figure
hold on
scatter(stim_d(~inverse_map),post_d(~inverse_map),10,[0 0 0]  )
scatter(stim_d(both_sig),post_d(both_sig),30,[1 0 0])
scatter(stim_d(stim_res),post_d(stim_res),30,[0 0 1])
scatter(stim_d(post_res),post_d(post_res),30,[0 1 0])
hold off

xlabel('Stim vs Baseline (d)')
ylabel('Post-stim vs Baseline (d)')
legend({'ns','both sig', 'stim sig only','post sig only'})
title(sprintf('During Stim vs Post Stim Increase in Coherence n=%.d',length(subjects)))
[rho,pval] = corr(agg(:,1),agg(:,2),"Type","Spearman");
result = sprintf("R = %.4f\np = %.6f",rho,pval);
text(0,0.6,result)

%% Ipsi vs Contra Results in stim_d vs post_d
sideIdx = ismember({corrdat.rec_side},'ipsi');
sideIdx = sideIdx(both_sig);
stim_d_s = stim_d(both_sig);
post_d_s = post_d(both_sig);
figure('Name','Ipsi v Contra All Responses')
subplot(2,1,1)
hold on
scatter(stim_d_s(sideIdx),post_d_s(sideIdx),20,'blue')
title('All Significant Responses Ipsilateral to Stim Site')
xlabel('Stim vs Baseline (d)')
ylabel('Post-stim vs Baseline (d)')
[rho2,pval2] = corr(stim_d_s(sideIdx)',post_d_s(sideIdx)',"Type","Spearman");
result = sprintf("R = %.4f\np = %.4f",rho2,pval2);
xlim([0 1])
text(0,0.6,result)
hold off
subplot(2,1,2)
title('All Significant Responses Contralateral to Stim Site')
xlabel('Stim vs Baseline (d)')
ylabel('Post-stim vs Baseline (d)')
hold on
scatter(stim_d_s(~sideIdx),post_d_s(~sideIdx),20,'red')
[rho3,pval3] = corr(stim_d_s(~sideIdx)',post_d_s(~sideIdx)',"Type","Spearman");
result = sprintf("R = %.4f\np = %.4f",rho3,pval3);
text(0,0.6,result)
xlim([0 1])
% legend({'Ipsi','Contra'})
hold off
%% region plot
close all
both_sig = all(agg < alpha_sig, 2);
ROIs = {corrdat(both_sig).region};
pattern = '^ctx-(rh|lh)-';
ROIs = cellfun(@(x) regexprep(x, pattern, ''), ROIs, 'UniformOutput', false);
pattern = '^(Right|Left)-';
ROIs = cellfun(@(x) regexprep(x, pattern, ''), ROIs, 'UniformOutput', false);

all_rois = {corrdat.region};
pattern = '^ctx-(rh|lh)-';
all_rois = cellfun(@(x) regexprep(x, pattern, ''), all_rois, 'UniformOutput', false);
pattern = '^(Right|Left)-';
all_rois = cellfun(@(x) regexprep(x, pattern, ''), all_rois, 'UniformOutput', false);
all_dat = [corrdat.post_d];



colors = [[0 0 0]; [1,0,0]];
roi_cat = categorical(ROIs);
post_stim = [corrdat(both_sig).post_d];
thresh = 0.3; threshloc = post_stim > thresh;
thresh_stim = post_stim(threshloc);
thresh_roi = categorical(ROIs(threshloc));
figure
subplot(2,1,1)
% boxplot(post_stim,roi_cat);
swarmchart(roi_cat,post_stim,[],colors(threshloc+1,:));
ylabel("Post Stimulation Coherence (d)")
% swarmchart(post_stim,roi_cat,[],colors(threshloc+1,:));
title('All Significant Responses')
subplot(2,1,2)
swarmchart(thresh_roi,thresh_stim);
ylabel("Post Stimulation Coherence (d)")
title('Regions with d > 0.3')
% boxplot(thresh_stim,thresh_roi);
figure
histogram(categorical(ROIs(threshloc)),'DisplayOrder','ascend');
title("Brain Regions with Large Responses (d > 0.3)")
figure
histogram(roi_cat,'DisplayOrder','ascend')
title("All Brain Regions with d > 0.3")
%% Proportion of Significant Responses
close all
both_sig = all(agg < 0.05, 2);
ROIs = {corrdat(both_sig).region};
pattern = '^ctx-(rh|lh)-';
ROIs = cellfun(@(x) regexprep(x, pattern, ''), ROIs, 'UniformOutput', false);
pattern = '^(Right|Left)-';
ROIs = cellfun(@(x) regexprep(x, pattern, ''), ROIs, 'UniformOutput', false);
roi_cat = categorical(ROIs);
post_stim = [corrdat(both_sig).post_d];
thresh = 0.3; threshloc = post_stim > thresh;
thresh_stim = post_stim(threshloc);
thresh_roi = categorical(ROIs(threshloc));
allname = {corrdat.region};
pattern = '^ctx-(rh|lh)-';
allname = cellfun(@(x) regexprep(x, pattern, ''), allname, 'UniformOutput', false);
pattern = '^(Right|Left)-';
allname = cellfun(@(x) regexprep(x, pattern, ''), allname, 'UniformOutput', false);
[allct, allname] = histcounts(categorical(allname));

figure %
bar(allname,allct)
[ct, name] = histcounts(roi_cat);
uniques = unique(name);
loc = ismember(allname,name);
allname = allname(loc);
[allname, idx1] = sort(allname);
[name, idx2] = sort(name);
allct_target = allct(loc);
allct_target = allct_target(idx1);
ct = ct(idx2);

figure %
bar(allname,allct_target);
hold on
bar(name,ct);
hold off
title('Number of Responses within regions')
legend({'All Responses','Significant Responses'})
ratio = ct ./ allct_target;
[ratio, id] = sort(ratio);
sortname = name(id);

figure %
bar(sortname,ratio)
title('Proportion of Resonant Responses to all Responses Within a Brain Region')

loc4 = ismember(sortname,thresh_roi);
figure %
bar(sortname(loc4),ratio(loc4))
title('Proportion of Large Resonant Responses to all Responses Within a Brain Region')
%% Multiple Regression of Factors
Y = [corrdat.post_d]';
features = {"freq" "current" "pw" "loc"};% "stim_d" "charge" "chargePerPhase"};
X = zeros(length(Y),length(features)+1);
regs = struct;
for i=1:length(features)
    regs(i).feature = features{i};
    feat = [corrdat.(features{i})]';
    XX = [feat, ones(numel(feat),1)];
    % [regs(i).b,regs(i).bint,~,~,regs(i).stats] = regress(Y,XX);
    [regs(i).b,~,~,~,regs(i).stats] = regress(Y,XX);
    X(:,i) = categorical(feat);
end
X(:,end) = ones(length(Y),1);
[b,bint,r,rint,stats] = regress(Y,X);

%% Effect of Parameters Across Different Regions via PCA, Clumped Regions
conditions = {'stimregion','current','pw','freq'};
stim_d = [corrdat.stim_d];
post_d = [corrdat.post_d];
agg = [stim_d; post_d]';
both_sig = all(agg < alpha_sig, 2);
[Acoeff,Ascore,Alatent,Atsquared,Aexplained,Amu] = pca(agg);
Aproj = Acoeff(:,1);
Aproj_all = agg * Aproj;
figs = struct;
for i=1:numel(conditions)
figs(i).figure = plotPCAConditionAcrossRegions(Aproj_all,{corrdat.pooled_region},{corrdat.(conditions{i})},unique(regionTable.region),post_sig);
figs(i).condition = conditions{i};
end

%%
conditions = {'stimregion','pw','freq','current'};
stim_d = [corrdat.stim_d];
post_d = [corrdat.post_d];
agg = [stim_d; post_d]';
both_sig = all(agg < alpha_sig, 2);
[Acoeff,Ascore,Alatent,Atsquared,Aexplained,Amu] = pca(agg);
Aproj = Acoeff(:,1);
Aproj_all = agg * Aproj;
figs = struct;
for i=1:numel(conditions)
figs(i).figure = plotPCAConditionAcrossRegions(Aproj_all,{corrdat.region_class},{corrdat.(conditions{i})},targetSubregions,sig_ipsi,'red');
figs(i).condition = conditions{i};
end
%%
targetRegions = {'thalamus', 'hippocampal complex','amygdala','cingulate','temporal','ofc'};
targetSubregions = unique(regionTable.subregion(ismember(regionTable.region,targetRegions)));
for i=1:numel(conditions)
% figs(i).figure = plotPCAConditionAcrossRegions(Aproj_all,{corrdat.pooled_region},{corrdat.(conditions{i})},targetRegions,post_sig,'blue');
cat_hist(post_d,{corrdat.region_class},{corrdat.(conditions{i})},targetSubregions,post_sig);
figs(i).condition = conditions{i};
end
%%
% close all
conditions = {'stimregion','pw','freq'};
ipsi2stim = ismember({corrdat.rec_side},'ipsi');
% ipsi2stim = ismember({corrdat.rec_side},'contra');
% titlestring = "All Contralateral Channels to Stimulation, Pooled to Common Functional Region";

sig_ipsi = [post_sig ipsi2stim'];
sig_ipsi = all(sig_ipsi == 1, 2);
sig_ipsi = ipsi2stim;
titlestring = "All Ipsilateral Channels to Stimulation, Pooled to Common Functional Region";
% titlestring = "Significant Ipsilateraclose l Channels to Stimulation, Pooled to Common Functional Region";
for i=1:length(conditions)
% [~, tcl]=plotConditionAcrossRegions(corrdat,{corrdat.region_class},targetSubregions,conditions{i},sig_ipsi);
% [~, tcl]=plotConditionAcrossRegions(corrdat,{corrdat.pooled_region},targetRegions,conditions{i},sig_ipsi);
[~, tcl]=plotConditionAcrossRegions(corrdat,{corrdat.pooled_region},targetRegions,conditions{i},sig_ipsi);

title(tcl,sprintf("%s\nvaried %s",titlestring,conditions{i}))
% conditions = {'stimregion','current','pw','freq'};
end
%%
for i=1:length(conditions)
% plotConditionAcrossRegions(corrdat,{corrdat.pooled_region},unique(regionTable.region),conditions{i},sig_ipsi);
plotConditionAcrossRegions(corrdat,{corrdat.pooled_region},unique(regionTable.region),conditions{i},0);


end
%%
function fig = plotPCAConditionAcrossRegions(PCAdata,RegionData, ...
    condData,ROIs,indexer,color)

if sum(indexer) > 0
    PCAdata = PCAdata(indexer);
    condData = condData(indexer);
    RegionData = RegionData(indexer);
end
% figure('Name',figName,'Visible','on');
fig = figure;
rows = ceil(sqrt(length(ROIs)));
cols = ceil(length(ROIs)/rows);
tcl = tiledlayout(rows,cols,"TileSpacing","compact");
for i=1:numel(ROIs)
    locs = ismember(RegionData,ROIs{i});
    data = PCAdata(locs); 
    conds = {condData{locs}};
    numflags = cellfun(@isnumeric, conds);
    ax=nexttile(tcl);
   

    if any(numflags)
        conds = cell2mat(conds);
        scatter(ax,data,conds,color);
    else
        scatter(ax,data,categorical(conds),color);
    end
    
    title(ax,ROIs{i})
end
linkaxes(tcl.Children,'x')
end
function [fig, tcl] = plotConditionAcrossRegions(dataStruct,RegionData,ROIs,condition,indexer)
colors = [0 0,1.0000;1.0000,0,0;0,1.0000,0;0,0,0.1724;1.0000,0.1034,0.7241;1.0000,    0.8276         0;    0    0.3448         0;    0.5172    0.5172    1.0000;    0.6207    0.3103    0.2759;    0    1.0000    0.7586;    0    0.5172    0.5862;    0         0    0.4828;    0.5862    0.8276    0.3103;    0.9655    0.6207    0.8621;    0.8276    0.0690    1.0000;    0.4828    0.1034    0.4138;    0.9655    0.0690    0.3793;    1.0000    0.7586    0.5172;    0.1379    0.1379    0.0345;    0.5517    0.6552    0.4828;    0.9655    0.5172    0.0345;    0.5172    0.4483         0;    0.4483    0.9655    1.0000;    0.6207    0.7586    1.0000];
x = [dataStruct.stim_d];
y = [dataStruct.post_d];
condData = {dataStruct.(condition)};
numflag = 0;
if any(cellfun(@isnumeric, condData))
    numflag = 1;
    condData = cell2mat(condData);
end
allconds = unique(condData);
N = length(allconds);
cmap = colors(1:N,:);


if sum(indexer) > 0
    x= x(indexer);
    y= y(indexer);
    condData = condData(indexer);
    RegionData = RegionData(indexer);
end
% figure('Name',figName,'Visible','on');
fig = figure;
rows = ceil(sqrt(length(ROIs)));
cols = ceil(length(ROIs)/rows);
tcl = tiledlayout(rows,cols,"TileSpacing","compact");
for i=1:numel(ROIs)
    locs = ismember(RegionData,ROIs{i});
    conds = condData(locs);
    ax=nexttile(tcl);
    hold on
    for j=1:N
        if numflag
        conLoc = ismember(conds,allconds(j));
        else
        conLoc = ismember(conds,allconds{j});
        end
        xx = x(conLoc);
        yy = y(conLoc);
        % color = repmat(cmap(j,:)x,length(yy),1);
        scatter(ax,xx,yy,20,cmap(j,:),'filled','MarkerFaceAlpha',0.3);
    end
    hold off

    title(ax,ROIs{i})
end
if numflag
    leg = cellfun(@num2str,num2cell(allconds),'UniformOutput', false);
else
    leg = allconds;
end
hL = legend(leg);
fontsize(hL,20,'points')
hL.Layout.Tile='South';
linkaxes(tcl.Children,'xyz')
xlabel(tcl,'Stim vs Baseline Coherence (d)')
ylabel(tcl,'Post stim vs Baseline Coherence (d)')
xlim([-0.1 1])
% ylim([-0.1 1])
end
function [channel, polarity] = parse_stimpair(stimpair)
% takes stimpair string of format AR10(-)_AR9(+) 
% ie shank, Number, (polarity) and returns the anode and cathode names and
% polarities in that order
str = split(stimpair, '_');
channel = cell(2,1);
polarity = channel;
for i=1:2
pattern = '([A-Za-z]+\d+)\(([+-])\)';
tokens = regexp(str{i}, pattern, 'tokens');

if ~isempty(tokens)
    channel{i} = tokens{1}{1};
    polarity{i} = tokens{1}{2};
else
    channel{i} = [];
    polarity{i} = [];
end
end
[polarity, idx] = sort(polarity);
channel = channel(idx);
end

function chans = remakeChannels(channels,shanks)
chans = cell(length(channels),1);
for i=1:length(shanks)
    strs = split(shanks{i},'_');
    patt = '\d+';
    chan = channels{i};
    num = regexp(chan,patt);

    chans{i} = strcat(strs{1},strs{2},chan(num));
end
end

function out = checkIpsiShank(shank, stimSide)
    % Split the 'shank' cell array and extract the last part of each entry
    sides = cellfun(@(x) x{end}, cellfun(@(x) split(x, '_'), shank, 'UniformOutput', false), 'UniformOutput', false);
    
    % Compare each extracted side with the corresponding entry in 'stimSide'
    isIpsi = cellfun(@strcmp, sides, stimSide);
    
    % Allocate the output based on the comparison results
    out = cell(size(shank));
    out(isIpsi) = {'ipsi'};
    out(~isIpsi) = {'contra'};
end

function out = identifyStimRegion(anode,cathode)

out = cell(size(anode));
for i = 1:length(anode)
    if strcmpi(anode{i},cathode{i})
        o = split(anode{i},'-');
        out{i} = strjoin({o{1},'amyg'},'-');
    else
        t = strjoin({cathode{i},anode{i}},'-');
        strs = split(t,'-');
        strs = unique(strs);
        strs = sort(strs);
        o = strjoin({strs{1:end-1},'amyg'},'-');
        out{i} = o;
    end
end
end

function fig = cat_hist(data,RegionData, ...
    condData,ROIs,indexer)
colors = [1 0 0; 0 0 1; 0 1 0; 1 1 0; 0 1 1; 1 0 1];
if sum(indexer) > 0
    data = data(indexer);
    condData = condData(indexer);
    RegionData = RegionData(indexer);
end
% figure('Name',figName,'Visible','on');
fig = figure;
rows = ceil(sqrt(length(ROIs)));
cols = ceil(length(ROIs)/rows);
tcl = tiledlayout(rows,cols,"TileSpacing","compact");
for i=1:numel(ROIs)
    locs = ismember(RegionData,ROIs{i});
    dat = data(locs); 
    conds = {condData{locs}}';
    numflags = cellfun(@isnumeric, conds, 'UniformOutput',false);
    numflags = cell2mat(numflags);
    ax=nexttile(tcl);
    
    
    if any(numflags)
        conds = cell2mat(conds);
        condLabs = unique(conds);
        condLabs = num2cell(condLabs);
        condStr = cellfun(@num2str,condLabs,'UniformOutput',false);
    else
        condLabs = unique(conds);
        condStr = condLabs;
    end
    
    hold on
    for j=1:length(condLabs)
        
        tvar = condStr{j};
        if numflags
            tvar = str2double(tvar);
        end
        loc = ismember(conds,tvar);
        pdat = dat(loc);
        histogram(ax,pdat,'FaceColor',colors(j,:),'FaceAlpha',0.75,'BinWidth',0.05)
    end
    hold off
    title(ax,ROIs{i})
end
leg = unique(condData);
linkaxes(tcl.Children,'x')
hL = legend(leg);
fontsize(hL,20,'points')
hL.Layout.Tile='South';
ax.YScale ='log';



end

function brain = verifyElectrodes(brain,electrodes)
electrodeNames=brain.electrodes.Name;
locs = ismember(electrodeNames,electrodes);
brain.electrodes.DefinitionIdentifier = brain.electrodes.DefinitionIdentifier(locs);
brain.electrodes.Annotation = brain.electrodes.Annotation(locs);
brain.electrodes.Label = brain.electrodes.Label(locs);
brain.electrodes.Name = brain.electrodes.Name(locs);
brain.electrodes.Location = brain.electrodes.Location(locs,:);   
end