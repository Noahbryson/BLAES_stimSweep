boxpath = '/Users/nkb/Library/CloudStorage/Box-Box/Brunner Lab/DATA/BLAES/BLAES_param';
saveroot = "/Users/nkb/Documents/NCAN/patients";
locs = dir(fullfile(boxpath,"UIC*"));
subjects = {locs.name};
for i=1:length(subjects)
    subpath = fullfile(saveroot,subjects{i});
    if ~exist(subpath,"dir")
        mkdir(subpath)
    end
    savepath = fullfile(subpath,"electrodes");
    if ~exist(savepath,"dir")
        mkdir(savepath)
    end
    cmap_fp = fullfile(boxpath,subjects{i},"ChannelMap.mat");
    bcidir = dir(fullfile(boxpath,subjects{i},'*.dat'));
    bcifp = fullfile(bcidir.folder,bcidir.name);
    utahElectrodes2dat(cmap_fp,bcifp,savepath);
end
