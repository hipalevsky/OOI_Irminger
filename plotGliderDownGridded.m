function plotGliderDownGridded(plotting,inddown_good,fignum)

%Inputs are the glider data structure and the indices of down profiles with
%real data
     
figure(fignum); clf
for i = 1:length(inddown_good)
    subplot(2,ceil(length(inddown_good)/2),i)
    plot(plotting.O2conc(:,inddown_good(i)),plotting.depth_grid,'b'); hold on;
    plot(plotting.O2conc(:,inddown_good(i)-1),plotting.depth_grid,'k'); hold on;
    plot(plotting.O2conc(:,inddown_good(i)+1),plotting.depth_grid,'k'); hold on;
    set(gca,'YDir','reverse'); axis([270 290 0 1000]); ylabel('Depth (m)'); xlabel('O_2 concentration (raw)');
end

end