%% Visualize individual profiles to compare up and down profiles

daystart = [datenum(2014,9,15) datenum(2015,1,15) datenum(2015,4,1)]; %input the day to pull profiles from for comparison
timewindow = 2; %window to pull profiles from
    numcomp = length(daystart);

figure(1); clf
for j = 1:numcomp
    indup = find(Yr1_wfpgrid.time_start > daystart(j) & Yr1_wfpgrid.time_start < daystart(j) + timewindow & Yr1_wfpgrid.updown' > 0);
    inddown = find(Yr1_wfpgrid.time_start > daystart(j) & Yr1_wfpgrid.time_start < daystart(j) + timewindow & Yr1_wfpgrid.updown' < 0);

    subplot(numcomp,4,1 + 4*(j-1))
    for i = 1:length(indup)
        plot(Yr1_wfpgrid.O2conc(:,indup(i)),Yr1_wfpgrid.depth_grid,'b.','markersize',10,'linewidth',2); hold on;
    end
    for i = 1:length(inddown)
        plot(Yr1_wfpgrid.O2conc(:,inddown(i)),Yr1_wfpgrid.depth_grid,'r.','markersize',10,'linewidth',2); hold on;
    end
    set(gca,'YDir','reverse'); axis([240 290 100 2650]); ylabel('Depth (m)'); xlabel('O_2 concentration (raw)'); title(datestr(daystart(j)));
    
    subplot(numcomp,4,2 + 4*(j-1))
    for i = 1:length(indup)
        plot(diff(Yr1_wfpgrid.O2conc(:,indup(i))),Yr1_wfpgrid.depth_grid(2:end),'b.','markersize',10,'linewidth',2); hold on;
    end
    for i = 1:length(inddown)
        plot(diff(Yr1_wfpgrid.O2conc(:,inddown(i))),Yr1_wfpgrid.depth_grid(2:end),'r.','markersize',10,'linewidth',2); hold on;
    end
    set(gca,'YDir','reverse'); axis([-5 5 100 2650]); ylabel('Depth (m)'); xlabel('O_2 concentration gradient'); title(datestr(daystart(j)));

    subplot(numcomp,4,3 + 4*(j-1))
    for i = 1:length(indup)
        plot(Yr1_wfpgrid.T(:,indup(i)),Yr1_wfpgrid.depth_grid,'b.','markersize',10,'linewidth',2); hold on;
    end
    for i = 1:length(inddown)
        plot(Yr1_wfpgrid.T(:,inddown(i)),Yr1_wfpgrid.depth_grid,'r.','markersize',10,'linewidth',2); hold on;
    end
    set(gca,'YDir','reverse'); axis([1 5 100 2650]); ylabel('Depth (m)'); xlabel('T (deg C)'); title(datestr(daystart(j))); legend('location','northwest','Up','Down');

    subplot(numcomp,4,4 + 4*(j-1))
    for i = 1:length(indup)
        plot(Yr1_wfpgrid.S(:,indup(i)),Yr1_wfpgrid.depth_grid,'b.','markersize',10,'linewidth',2); hold on;
    end
    for i = 1:length(inddown)
        plot(Yr1_wfpgrid.S(:,inddown(i)),Yr1_wfpgrid.depth_grid,'r.','markersize',10,'linewidth',2); hold on;
    end
    set(gca,'YDir','reverse'); axis([34.84 34.92 100 2650]); ylabel('Depth (m)'); xlabel('S'); title(datestr(daystart(j)));

end

%% Plot histograms of entire record comparing up and down plots

indup = find(Yr2_wfpgrid.updown > 0);
inddown = find(Yr2_wfpgrid.updown < 0);

O2grad_up = -5*diff(Yr2_wfpgrid.O2conc(:,indup));
O2grad_down = 5*diff(Yr2_wfpgrid.O2conc(:,inddown));
O2conc_up = Yr2_wfpgrid.O2conc(:,indup);
O2conc_down = Yr2_wfpgrid.O2conc(:,inddown);
        tol = 10;
    upcut = find(abs(O2grad_up) < tol); 
    downcut = find(abs(O2grad_down) < tol);
        minval = 240; maxval = 290;
    upcutconc = find(O2conc_up > minval & O2conc_up < maxval); 
    downcutconc = find(O2conc_down > minval & O2conc_down < maxval); 
    
figure(2); clf
    subplot(221)
histogram(O2conc_up(upcutconc)); hold on;
histogram(O2conc_down(downcutconc)); hold on;
title('O_2 concentration, \muM (raw)'); legend('Up','Down')
    subplot(222)
histogram(O2grad_up(upcut)); hold on;
histogram(O2grad_down(downcut)); hold on;
title('O_2 gradient (dO_2/dz), \muM/m'); legend('Up','Down')
    subplot(223)
histogram(Yr2_wfpgrid.T(:,indup)); hold on;
histogram(Yr2_wfpgrid.T(:,inddown)); hold on;
title('T, deg C'); legend('Up','Down')
    subplot(224)
histogram(Yr2_wfpgrid.S(:,indup)); hold on;
histogram(Yr2_wfpgrid.S(:,inddown)); hold on;
xlim([34.8 34.95])
title('S'); legend('Up','Down')

figure(3); clf
bw = 1;
    subplot(211)
histogram(Yr2_wfpgrid.O2conc(1:261,indup),'BinWidth',bw); hold on;
histogram(Yr2_wfpgrid.O2conc(1:261,inddown),'BinWidth',bw); hold on;
xlim([240 290]); title('150-1450 meters, O2 concentration'); legend('Up','Down')
    subplot(212)
histogram(Yr2_wfpgrid.O2conc(262:end,indup),'BinWidth',bw); hold on;
histogram(Yr2_wfpgrid.O2conc(262:end,inddown),'BinWidth',bw); hold on;
xlim([240 290]); title('1455 - 2600 meters, O2 concentration'); legend('Up','Down')

