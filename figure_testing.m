%% Text mapping figure test
figure(1)
subplot(3,10,1)
[t,y] = generate_theta_burst_waveform(50,1);

plot(t,y);
ylim([-0.5,2.5])
xlim([0,0.2])
subplot(3,10,30)