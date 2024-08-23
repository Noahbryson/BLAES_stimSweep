
loc = '/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param/method_building';
files = dir(fullfile(loc,'*.mat'));
if ~contains("selected_data.mat",{files.name})
    theta_epochs = load(fullfile(loc,files(1).name));
    theta_epochs = theta_epochs.local_data;
    for i=2:length(files)
        temp = load(fullfile(loc,files(i).name));
        theta_epochs = [theta_epochs temp.local_data];
    end
else
    theta_epochs = load(fullfile(loc,"selected_data.mat"));
    load(fullfile(loc,"fs.mat"));
    theta_epochs = theta_epochs.exportStruct;
end
[~,sortIdx] = sort([theta_epochs.charge]);
theta_epochs = theta_epochs(sortIdx);
%% Theta Env
close all
theta_band = [4 10]; % Hz
smoothing_window = 0.1; % in seconds
filterOrder = computeStableFilter(theta_band,'bandpass',fs);
postStimStart = floor(2*fs);
for i=1:length(theta_epochs)
    [theta_epochs(i).theta_power,theta_epochs(i).theta_phase, theta_epochs(i).theta_filt] = timeseriesPower(theta_epochs(i).pre_stim_post',fs,theta_band,filterOrder, ...
        'smooth',smoothing_window,'baselineDuration',1);
    theta_epochs(i).theta_power = zscore(theta_epochs(i).theta_power);
    theta_epochs(i).theta_power = avg_baseline_correct(theta_epochs(i).theta_power,1,fs);
    theta_epochs(i).baseline_corr = channelCoherence(theta_epochs(i).baseline);
    theta_epochs(i).stim_corr = channelCoherence(theta_epochs(i).signals);
    theta_epochs(i).post_corr = channelCoherence(theta_epochs(i).pre_stim_post(:,postStimStart:end));
end
%%
chan = 100;
figure(1)
plot(theta_epochs(chan).theta_filt,'Color',[1 1 1 0.25])
hold on
plot(mean(theta_epochs(chan).theta_filt,2),'Color',[1 0 0])

figure(2)
plot(theta_epochs(chan).theta_power,'Color',[1 1 1 0.25])
hold on
plot(mean(theta_epochs(chan).theta_power,2),'Color',[1 0 0])

figure(3)
plot(theta_epochs(chan).pre_stim_post','Color',[1 1 1 0.25])
hold on
plot(mean(theta_epochs(chan).pre_stim_post),'Color',[1 0 0])
%% Compute Differences in Power


%% Visually Explore Data

channels = unique([theta_epochs.channel_idx]);
channel = channels(103);
datIdx = ismember([theta_epochs.channel_idx],channel);
% datIdx = 2471;
% plot_channel(theta_epochs(datIdx),fs,1)
plot_PSD(theta_epochs(datIdx),fs)
% plot_time_frequency(theta_epochs(datIdx),fs)
% plot_correlation(theta_epochs(datIdx))

%%
function plot_channel(data,fs,showCorrelation)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' raw');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
for i=1:length(data)
    ax = nexttile(tcl);
    y = data(i).pre_stim_post;
    if showCorrelation
        base = mean(data(i).baseline_corr);
        stim = mean(data(i).stim_corr);
        post = mean(data(i).post_corr);
        d = post - base;
        cor_str = sprintf('baseline %.3f, stim %.3f post %.3f',base,stim,post);
        tit = sprintf("%s\nd=%.3f\n%s",data(i).label,d,cor_str);
        title(tit)
    else
        title(data(i).label)
    end
    t = linspace(0,size(y,2)/fs,size(y,2));
    hold on
    plot(ax,t,y','Color',[1 1 1 0.25]);
    plot(ax,t,mean(y),'Color',[1 0 0])
end

linkaxes(tcl.Children);
end

function plot_correlation(data)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' correlation');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
for i=1:length(data)
    ax = nexttile(tcl);
    base = data(i).baseline_corr;
    stim = data(i).stim_corr;
    post = data(i).post_corr;
    % d = post - base;
    % cor_str = sprintf('baseline %.3f, stim %.3f post %.3f', ...
    %     mean(base),mean(stim),mean(post));
    % % tit = sprintf("%s\nd=%.3f\n%s",data(i).label,d,cor_str);
    tit = sprintf("%s\n %.3f %s",data(i).label,data(i).charge,data(i).chargeUnits);
    title(tit)


    hold on
    histogram(base,'FaceColor',[1 0 0],'FaceAlpha',0.6,'BinWidth',0.1)
    histogram(stim,'FaceColor',[0 1 0],'FaceAlpha',0.6,'BinWidth',0.1)
    histogram(post,'FaceColor',[0 0 1],'FaceAlpha',0.6,'BinWidth',0.1)
end
hL = legend({'pre','stim','post'});
fontsize(hL,20,'points')
hL.Layout.Tile='East';

linkaxes(tcl.Children);

end

function plot_time_frequency(data,fs)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' time-freq');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
F = linspace(1,fs,fs);
M = 50;
g = hann(M);
L = 0.5*M;
for i=1:length(data)
    ax = nexttile(tcl);
    y = mean(data(i).pre_stim_post,1);
    % t = linspace(0,length(y)/fs,length(y));
    
    [s,f,t] = spectrogram(y,g,L,F,fs);
    mesh(t,f,abs(s).^2)
    view(2), axis tight
    % h=imagesc(t,f,abs(s).^2);
    % ylim([0 30])
    tit = sprintf("%s\n %.3f %s",data(i).label,data(i).charge,data(i).chargeUnits);
    title(tit)
    ax.YScale = 'log';

end
linkaxes(tcl.Children);
% ylim([0 30])
end
function plot_PSD(data,fs)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' PSD');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
endloc = floor(2*fs);
stimloc = floor(fs);
for i=1:length(data)
    ax = nexttile(tcl);
    hold on
    
    Ndft = 128;
    leg = {'post','stim','stim-interp','pre'};
    % post stim
    y = mean(data(i).pre_stim_post(:,endloc:end),1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    [Pxx,f] = pwelch(y,g,L,[],fs);
    % Pxx = 20*log(Pxx);
    loglog(f,Pxx,'Color',[0 0 1])
    
    % stim
    y = mean(data(i).pre_stim_post(:,stimloc:endloc),1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    [Pxx,f] = pwelch(y,g,L,[],fs);
    % Pxx = 20*log(Pxx);
    loglog(f,Pxx,'Color',[1 0 0])


    % stim - noise corrected
    y = mean(data(i).signals,1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    [Pxx,f] = pwelch(y,g,L,[],fs);
    % Pxx = 20*log(Pxx);
    loglog(f,Pxx,'Color',[1 0 1 0.7])

    % pre stim
    y = mean(data(i).pre_stim_post(:,1:stimloc),1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    [Pxx,f] = pwelch(y,g,L,[],fs);
    % Pxx = 20*log(Pxx);
    loglog(f,Pxx,'Color',[0 1 0])
    tit = sprintf("%s\n %.3f %s",data(i).label,data(i).charge,data(i).chargeUnits);
    title(tit)
    ax.XScale = 'log';
    ax.YScale = 'log';
end
linkaxes(tcl.Children,'x');
hL = legend(leg);
fontsize(hL,20,'points')
hL.Layout.Tile='East';
% ylim([0 30])
end
%%
close all
% trajplot(dat,0,pwd,fs,[],"theta-power",'field','theta_power','avgField',1)
% trajplot(dat,0,pwd,fs,[],"raw",'field','pre_stim_post','avgField',1)
trajplot(theta_epochs,0,pwd,fs,[],"theta-filt",'field','theta_filt','avgField',1)
%%
close all
idx = [theta_epochs.current]==1;
dat1 = theta_epochs(idx);
trajplot(dat1,0,pwd,fs,[],"theta-power",'field','theta_power','avgField',1)

trajplot(dat1,0,pwd,fs,[],"theta-filt",'field','theta_filt','avgField',1)
%% Filt vs Interp
chan = 9;
test = getLowPassData(theta_epochs(chan).raw_sig',60,4,fs);
test = getHighPassData(test,2,8,fs);
close all
fig = figure(1);
subplot(2,1,1)
plot(test,'Color',[1 1 1 0.15])
hold on
plot(mean(test,2),'Color',[1 0 0])
subplot(2,1,2)
plot(theta_epochs(chan).signals','Color',[1 1 1 0.15])
hold on
plot(mean(theta_epochs(chan).signals),'Color',[1 0 0])