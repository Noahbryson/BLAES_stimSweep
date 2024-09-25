user = expanduser('~'); % Get local path for interoperability on different machines, function in my tools dir. 
if ispc
    boxpath = fullfile(user,"Box/Brunner Lab"); % Path to data
     BCI2KPath = "C:\BCI2000\BCI2000";
else
    boxpath =  fullfile(user,'Library/CloudStorage/Box-Box/Brunner Lab'); % Path to data
    BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
end
datapath = fullfile(boxpath,"/DATA/BLAES/BLAES_param");
addpath(genpath(fullfile(user,'Documents/NCAN/code/BLAES_stimSweep')));
addpath(genpath(fullfile(user,'Documents/NCAN/code/MATLAB_tools')));
rawfilesavepath = fullfile(boxpath,'/Posters/SfN2024/Raw Files');
loc = fullfile(datapath,'method_building');
files = dir(fullfile(loc,'*.mat'));
if ~contains("selected_data.mat",{files.name})
    data = load(fullfile(loc,files(1).name));
    data = data.local_data;
    for i=2:length(files)
        temp = load(fullfile(loc,files(i).name));
        data = [data temp.local_data];
    end
else
    data = load(fullfile(loc,"selected_data.mat"));
    load(fullfile(loc,"fs.mat"));
    data = data.exportStruct;
end
[~,sortIdx] = sort([data.charge]);
data = data(sortIdx);


%%
theta_band = [4 10]; % Hz
smoothing_window = 0.1; % in seconds
filterOrder = computeStableFilter(theta_band,'bandpass',fs);
postStimStart = floor(2*fs);
sample_num=100;
% [data(sample_num).theta_power,data(sample_num).theta_phase, data(sample_num).theta_filt] = timeseriesPower(data(sample_num).pre_stim_post',fs,theta_band,filterOrder, ...
%     'smooth',smoothing_window,'baselineDuration',1);
% samp_theta_power   = zscore(data(sample_num).theta_power);
% samp_theta_power   = avg_baseline_correct(data(sample_num).theta_power,1,fs);
samp_baseline_corr = channelCoherence(data(sample_num).baseline);
samp_stim_corr     = channelCoherence(data(sample_num).signals);
samp_post_corr     = channelCoherence(data(sample_num).pre_stim_post(:,postStimStart:end));

%%
chan = 2;
figure(1)
plot(data(chan).theta_filt,'Color',[1 1 1 0.25])
hold on
plot(mean(data(chan).theta_filt,2),'Color',[1 0 0])

figure(2)
plot(data(chan).theta_power,'Color',[1 1 1 0.25])
hold on
plot(mean(data(chan).theta_power,2),'Color',[1 0 0])

figure(3)
plot(data(chan).pre_stim_post','Color',[1 1 1 0.25])
hold on
plot(mean(data(chan).pre_stim_post),'Color',[1 0 0])

%% Visually Explore Data
close all
dark = 1;
channels = unique([data.channel_idx]);
channel = channels(3);
datIdx = ismember([data.channel_idx],channel);
% datIdx = 2471;
% fig5=plot_theta(data(datIdx),fs,1);
% fig1=plot_channel(data(datIdx),fs,0,0);
% fig2=plot_PSD(data(datIdx),fs,0,1);
% figTaper=plot_PSD(data(datIdx),fs,0,1);
% figGam = plot_gamma(data(datIdx),fs,1);

% fig3=plot_time_frequency(data(datIdx),fs,'FFT');
fig4=plot_correlation(data(datIdx));

%%
function fig=plot_theta(data,fs,dark)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' raw');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
for i=1:length(data)
    ax = nexttile(tcl);
    y = data(i).theta_power;
    
    title(data(i).label)
    t = linspace(0,size(y,1)/fs,size(y,1));
    hold on
    lim = max(abs(y(:)));
    if dark
    plot(ax,t,y,'Color',[1 1 1 0.25]);
    else
    plot(ax,t,y,'Color',[0 0 0 0.25]);
    end
    plot(ax,t,mean(y,2),'Color',[1 0 0],'LineWidth',3)
    ylim([-lim lim])
end

linkaxes(tcl.Children,'x');
end

function fig=plot_gamma(data,fs,dark)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' raw');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
for i=1:length(data)
    ax = nexttile(tcl);
    y = data(i).gamma_power;
    
    title(data(i).label)
    t = linspace(0,size(y,1)/fs,size(y,1));
    hold on
    lim = max(abs(y(:)));
    if dark
    plot(ax,t,y,'Color',[1 1 1 0.25]);
    else
    plot(ax,t,y,'Color',[0 0 0 0.25]);
    end
    plot(ax,t,mean(y,2),'Color',[1 0 0],'LineWidth',3)
    ylim([-lim lim])
end

linkaxes(tcl.Children,'x');
end

function fig=plot_channel(data,fs,showCorrelation,dark)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' raw');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
for i=1:length(data)
    ax = nexttile(tcl);
    set(ax,'ButtonDownFcn',@fig_from_subplot)
    y = data(i).full_trial;
    if showCorrelation
        base = mean(data(i).baseline_corr);
        stim = mean(data(i).stim_corr);
        post = mean(data(i).post_corr);
        d1 = post - base;
        d2 = stim - base;
        cor_str = sprintf('baseline %.3f, stim %.3f post %.3f',base,stim,post);
        tit = sprintf("%s\nstim d=%.3f post d=%.3f\n%s",data(i).label,d2,d1,cor_str);
        title(tit)
    else
        title(data(i).label)
    end
    t = linspace(0,size(y,2)/fs,size(y,2));
    hold on
    lim = max(abs(y(:)));
    if dark
    plot(ax,t,y','Color',[1 1 1 0.25]);
    else
    plot(ax,t,y','Color',[0 0 0 0.25]);
    end
    plot(ax,t,mean(y),'Color',[1 0 0],'LineWidth',1.5,'LineStyle','--')
    ylim([-lim lim])
end

linkaxes(tcl.Children,'x');
end

function fig=plot_correlation(data)

rows = ceil(sqrt(length(data)));
cols = ceil(length(data)/rows);
chan = unique({data.channel});
figName = strcat(chan,' correlation');
fig=figure('Name',figName{1},'Visible','on');
set(fig,'Position',[100 100 1980 1020])
tcl = tiledlayout(rows,cols);
for i=1:length(data)
    ax = nexttile(tcl);
    set(ax,'ButtonDownFcn',@fig_from_subplot)

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

function fig=plot_time_frequency(data,fs,method)

if method == "wavelet"
    waveFlag = 1;
else
    waveFlag = 0;
end

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
function fig=plot_PSD(data,fs,multitaper,ratio)
frequencies = linspace(1,200,200);
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
    set(ax,'ButtonDownFcn',@fig_from_subplot)

    Ndft = 128;
    if ratio
        leg = {'post','stim','stim-interp'};
    else
    
    leg = {'post','stim','stim-interp','pre'};
    end 
    % pre stim
    y = mean(data(i).pre_stim_post(:,1:stimloc),1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    if multitaper
        [Pxx_pre,f] = pmtm(y,3,M,frequencies,fs);
    else
        [Pxx_pre,f] = pwelch(y,g,L,frequencies,fs);
    end    % Pxx = 20*log(Pxx);
    if ~ratio
    loglog(f,Pxx_pre,'Color',[0 1 0])
    hold on
    end


    
    % post stim
    y = mean(data(i).pre_stim_post(:,endloc:end),1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);

    if multitaper
        [Pxx,f] = pmtm(y,3,M,frequencies,fs);
    else
        [Pxx,f] = pwelch(y,g,L,frequencies,fs);
    end
    % Pxx = 20*log(Pxx);
    if ratio
        semilogx(f,Pxx./Pxx_pre,'Color',[0 0 1])
        hold on
    else
        loglog(f,Pxx,'Color',[0 0 1])
    end

    % semilogx(f,Pxx,'Color',[0 0 1])

    % stim
    y = mean(data(i).pre_stim_post(:,stimloc:endloc),1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    if multitaper
        [Pxx,f] = pmtm(y,3,M,frequencies,fs);
    else
        [Pxx,f] = pwelch(y,g,L,frequencies,fs);
    end
    % Pxx = 20*log(Pxx);
    if ratio
        semilogx(f,Pxx./Pxx_pre,'Color',[1 0 0])
    else
        loglog(f,Pxx,'Color',[1 0 0])
    end



    % stim - noise corrected
    y = mean(data(i).signals,1);
    M = length(y);
    g = hann(floor(M*0.5));
    L = floor(0.25*M);
    if multitaper
        [Pxx,f] = pmtm(y,3,M,frequencies,fs);
    else
        [Pxx,f] = pwelch(y,g,L,frequencies,fs);
    end    % Pxx = 20*log(Pxx);
    if ratio
        semilogx(f,Pxx./Pxx_pre,'Color',[1 0 1 0.7])
    else
        loglog(f,Pxx,'Color',[1 0 1 0.7])
    end
    
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
trajplot(data,0,pwd,fs,[],"theta-filt",'field','theta_filt','avgField',1)
%%
close all
idx = [data.current]==1;
dat1 = data(idx);
trajplot(dat1,0,pwd,fs,[],"theta-power",'field','theta_power','avgField',1)

trajplot(dat1,0,pwd,fs,[],"theta-filt",'field','theta_filt','avgField',1)
%% Filt vs Interp
chan = 9;
test = getLowPassData(data(chan).raw_sig',60,4,fs);
test = getHighPassData(test,2,8,fs);
close all
fig = figure(1);
subplot(2,1,1)
plot(test,'Color',[1 1 1 0.15])
hold on
plot(mean(test,2),'Color',[1 0 0])
subplot(2,1,2)
plot(data(chan).signals','Color',[1 1 1 0.15])
hold on
plot(mean(data(chan).signals),'Color',[1 0 0])