



function [path,idx,onset] = sliceVideo(videoDir,v,onset,duration,idx)


startFrame = onset * v.FrameRate;
endFrame = (onset+duration)*v.FrameRate;
newDir = fullfile(videoDir,v.Name(1:end-4));
if ~isfolder(newDir)
    mkdir(newDir)
end

v1 = read(v,[startFrame,endFrame]);

end