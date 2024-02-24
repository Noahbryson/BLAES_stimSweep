function X = checkDir(dir)

if ~isfolder(dir)
    disp('No Directory Found');
    X=0;
else
    disp('Directory Found')
    X=1;
end
end