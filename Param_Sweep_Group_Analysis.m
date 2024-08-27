%% Param Sweep Group Analysis
addpath(genpath('/Users/nkb/Documents/NCAN/code/BLAES_stimSweep'))
addpath(genpath('/Users/nkb/Documents/NCAN/code/MATLAB_tools'))
BCI2KPath = '/Users/nkb/Documents/NCAN/BCI2000tools';
bci2ktools(BCI2KPath);
boxpath = ('/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param');
groupPath = fullfile(boxpath,'group');
files = dir(fullfile(groupPath,'*_cohens.mat'));
corrdat = load(fullfile(groupPath,files(1).name));
corrdat = corrdat.outstruct;
for i=2:length(files)
    temp = load(fullfile(groupPath,files(i).name));
    temp = temp.outstruct;
    corrdat = [corrdat temp];
end
clear temp
%% TODO: attach brain region to each row in datastruct, then investigate what points relate to region.
% implicate a network in the phenomena i am observing!
%%
stim_d = [corrdat.stim_d];
post_d = [corrdat.post_d];
cmap = [[1 1 1]; [1 0 0]; [0 1 0]; [0 0 1]];
p = 0.05;
agg = [[corrdat.stim_p];[corrdat.post_p]]';
both_sig = all(agg < p, 2);
stim_sig = agg(:,1) < p ;
stim_nonsig = agg(:,2) >= p;
stim_res = stim_sig & stim_nonsig;
post_sig = agg(:,2) < p;
post_nonsig = agg(:,1) >= p;
post_res = post_sig & post_nonsig;
sig_map = both_sig + 2*stim_res + 3*post_res+1;
colors = cmap(sig_map,:);
% inverse_map = (both_sig + stim_res + post_res -1 ) *-1;
inverse_map = logical(sig_map -1);
figure
hold on
scatter(stim_d(~inverse_map),post_d(~inverse_map),10,[1 1 1]  )
scatter(stim_d(both_sig),post_d(both_sig),20,[255 179 186]/256)
scatter(stim_d(stim_res),post_d(stim_res),20,[186 255 201]/256)
scatter(stim_d(post_res),post_d(post_res),20,[186 225 255]/256)
hold off

xlabel('Stim vs Baseline (d)')
ylabel('Post-stim vs Baseline (d)')
legend({'ns','both sig', 'stim sig only','post sig only'})