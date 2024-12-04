function trajplot(epoch_datastruct,saveFigs,dir,fs,colors)
colors = [0 0,1.0000;1.0000,0,0;0,1.0000,0;0,0,0.1724;1.0000,0.1034,0.7241;1.0000,    0.8276         0;    0    0.3448         0;    0.5172    0.5172    1.0000;    0.6207    0.3103    0.2759;    0    1.0000    0.7586;    0    0.5172    0.5862;    0         0    0.4828;    0.5862    0.8276    0.3103;    0.9655    0.6207    0.8621;    0.8276    0.0690    1.0000;    0.4828    0.1034    0.4138;    0.9655    0.0690    0.3793;    1.0000    0.7586    0.5172;    0.1379    0.1379    0.0345;    0.5517    0.6552    0.4828;    0.9655    0.5172    0.0345;    0.5172    0.4483         0;    0.4483    0.9655    1.0000;    0.6207    0.7586    1.0000];
if ~exist(dir,'dir')
    mkdir(dir);
end
trajectories = unique({epoch_datastruct.shank});
traj_idx = {epoch_datastruct.shank};
newLoop = 1;
for tr=1:length(trajectories)
    locs = strcmp(traj_idx,trajectories{tr});
    fig=figure('Visible','on');
    set(fig,"PaperSize",[8 11]);
    fig.PaperPosition = [0 0 8 11];
    set(gcf,'Position',[100 100 3000 2000])
    tcl = tiledlayout(4,4);
    % traj_channels = unique({traj_data.channel});
    traj_channels = unique({epoch_datastruct(locs).channel});
    for l=1:length(traj_channels)
        chan = traj_channels{l};
        local_data = epoch_datastruct(strcmp({epoch_datastruct.channel},chan));
        ax=nexttile(tcl);
        set(gca,'ButtonDownFcn',@fig_from_subplot)
        hold on
        for j=1:length(local_data)
            y = local_data(j).avg + j;
            stdev = local_data(j).std;
            t = linspace(-local_data(j).timeOffset,local_data(j).timeOffset + length(y)/fs/2,length(y));
            plot(ax,t,y,'LineWidth',2,'DisplayName',local_data(j).label,'Color',colors(j,:));
        end
        hold off
        tit = title(chan,'Interpreter','none');
        fontsize(tit,24,'points')
    end
    hL = legend({local_data.label});
    fontsize(hL,20,'points')
    hL.Layout.Tile='East';
    if saveFigs
        exportgraphics(fig,fullfile(dir,sprintf("%s traj.pdf",trajectories{tr})),"Append",false);

        if newLoop
            exportgraphics(fig,fullfile(dir,'trajs.pdf'),"Append",false);
            newLoop = 0;
        else
            exportgraphics(fig,fullfile(dir,'trajs.pdf'),"Append",true);
        end
    end

end
linkaxes(tcl.Children, 'x');
end