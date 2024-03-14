%% Video Manipulation
% takes a video and adds a patch in the lower left corner every ten seconds
% for photodiode triggering. 
root = 'C:\Paradigms';
stimuliDir = fullfile(root,'\tasks\BLAES\BLAES_param_sweep\stimuli'); % path to folder containing stimuli
checkDir(stimuliDir);
videos = dir(strcat(stimuliDir,'\','*.mp4'));
video_num = 4;
parmName_dec = videos(video_num).name(1:20);
vidPath = fullfile(videos(video_num).folder,videos(video_num).name);

savePath = fullfile('stimuli',sprintf('test_%s',videos(video_num).name));
%%
video = VideoReader(vidPath);
%%
duration = video.Duration;
fps = video.FrameRate;
numFrame = video.NumFrames;
factors = factor(numFrame);
bigFactor = factors(end)*factors(end-1);
split = numFrame/bigFactor;
h = video.Height;
w = video.Width;
%%
writer = VideoWriter(savePath, 'MPEG-4');
writer.FrameRate=fps;
open(writer);
framecount = 0;
for i=1:bigFactor
    current = read(video,[framecount+1, framecount+split]);
    
    if mod(i,3) == 0 
        IMG = addPhotoPatch(current,255,h,w,endFrame-start);
    else
        IMG = addPhotoPatch(current,1,h,w,endFrame-start);
    end
    framecount = framecount+split;
    writeVideo(writer,IMG);
    fprintf('%d frame written, %d minutes of video\n',framecount,framecount/fps/60)
end

%%
close(writer);


%%
function IMG = addPhotoPatch(IMG,color,height,width,numFrames)
hspan = 0.9*height:1:height;
wspan = 1:1:0.1*width;
IMG(hspan,wspan,:,:) = color;
end