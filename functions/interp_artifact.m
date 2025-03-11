function signalOut = interp_artifact(signal,spikelocs,fs,interp_duration)
signalOut = signal;
interp_width = floor((interp_duration/1000)*fs);

for i=1:length(spikelocs)
    idx = spikelocs(i);
    onsetVal = idx-floor(interp_width/2);
    offsetVal = onsetVal + interp_width;
    % if offsetVal < 0 || onsetVal < 0
    %     disp('issue')
    % end
    if offsetVal <= length(signal) && onsetVal > 0
        replace = linspace(signal(onsetVal),signal(offsetVal),offsetVal-onsetVal+1);
        signalOut(onsetVal:offsetVal) = replace;

    end
end

% figure(1)
% aa1=subplot(2,1,1);
% plot(signal);
% hold on
% scatter(spikelocs,signal(spikelocs));
% hold off
% aa2=subplot(2,1,2);
% plot(signalOut,'LineWidth',2);
% hold on
% scatter(spikelocs,signalOut(spikelocs));
% hold off



end