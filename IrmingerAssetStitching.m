%% Stitch together assets for Years 1 and 2
% Points for comparison:
    % mixed layer (which includes fixed depth assets most of year)
    % somewhere in the 150-1000 m range where glider + WFP overlap
% Options for initial calibration:
    % 1. Do separate Winkler gain corrections for each asset and then
    % compare
    % 2. Do single Winkler gain correction on best-calibrated asset and
    % then carry that correction to all other assets at beginning of year
    % 3. Some sort of optimization that combines both of these (?)
    
%% Pin everything to stable deep ocean measurement in deep WFP data - plotting on isotherms
% Use rate of change on deep isotherms to correct for drift in WFP O2 and salinity
figure(1); clf
thermplot = Yr1_wfpgrid_therm.therm_grid(21:4:47); %select stable deep isotherms
thermstr = {'2.1 deg C', '2.3 deg C', '2.5 deg C', '2.7 deg C', '2.9 deg C', '3.1 deg C', '3.3 deg C'};
clear C h; bot = nicecolor('rrywwwwww'); top = nicecolor('rrykkkkkk');
C = [linspace(top(1),bot(1),length(thermplot))' linspace(top(2),bot(2),length(thermplot))' linspace(top(3),bot(3),length(thermplot))'];
yearspan = {'2014-2015','2015-2016'};
M = 15;

for i = 1:2
    if i == 1
        plotting = Yr1_wfpgrid_therm;
        plotting.time_start = Yr1_wfpgrid_therm.time_start(Yr1_wfpgrid_therm.ind_pair);
        O2gaincorr = plotting.O2conc*gain_wfp(1); %note need a better format for saving gain
        %Initialize arrays to hold slopes
            O2slope = NaN*ones(length(plotting.therm_grid),2);
            O2slope_err = NaN*ones(length(plotting.therm_grid),2);
            Tslope = NaN*ones(length(plotting.therm_grid),2);
            Tslope_err = NaN*ones(length(plotting.therm_grid),2);
            Sslope = NaN*ones(length(plotting.therm_grid),2);
            Sslope_err = NaN*ones(length(plotting.therm_grid),2);
    elseif i == 2
        plotting = Yr2_wfpgrid_therm;
        plotting.time_start = Yr2_wfpgrid_therm.time_start(Yr2_wfpgrid_therm.ind_pair);
        O2gaincorr = plotting.O2conc*gain_wfp(2); %note need a better format for saving gain
    end
    
%Set up to calculate varying drift rate over time
tstep = 90; %number of days in each time period for derivative plotting
tmin = min(plotting.time_start); tmax = max(plotting.time_start);
t = [floor(tmin)+1:ceil(tmax)];
tgrid = [floor(tmin)+1:tstep:tmax]-floor(tmin); tgrid = [tgrid max(t)-min(t)+1];
tgrid2 = [1 60 max(t)-min(t)+1];

%%% Look for trends over time in individual isotherms
for k = 1:length(plotting.therm_grid)
    %ind = find(~isnan(O2gaincorr(k,:)));
    ind = find(O2gaincorr(k,:) > nanmean(O2gaincorr(k,:)) - 2*nanstd(O2gaincorr(k,:)) & O2gaincorr(k,:) < nanmean(O2gaincorr(k,:)) + 2*nanstd(O2gaincorr(k,:)));
    [P,S] = polyfit(plotting.time_start(ind), O2gaincorr(k,ind)',1);
    O2slope(k,i) = P(1); O2int(k) = P(2);
    if sum(size(S.R)) == 4
        err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
        O2slope_err(k,i) = err(1); O2int_err(k) = err(2);
        % Calculate exponential fit
        if k >=12 & k<=82
            f = fit(plotting.time_start(ind) - min(plotting.time_start), O2gaincorr(k,ind)','exp2');
            fitvalues = coeffvalues(f);
            coefa(k,i) = fitvalues(1); coefb(k,i) = fitvalues(2); coefc(k,i) = fitvalues(3); coefd(k,i) = fitvalues(4);
        end
    end
    for j = 1:length(tgrid2)-1
        tind = find(plotting.time_start >= tgrid2(j) + floor(tmin) & plotting.time_start < tgrid2(j+1) + floor(tmin));
        intind = intersect(tind,ind);
        numpts_pieceslope(k,j,i) = length(intind);
        P = polyfit(plotting.time_start(intind), O2gaincorr(k,intind)',1);
        O2pieceslope(k,j,i) = P(1);
    end
    [P,S] = polyfit(plotting.time_start(ind), plotting.T(k,ind)',1);
    Tslope(k,i) = P(1); Tint(k) = P(2);
    if sum(size(S.R)) == 4
        err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
        Tslope_err(k,i) = err(1); Tint_err(k) = err(2);
    end
    [P,S] = polyfit(plotting.time_start(ind), plotting.S(k,ind)',1);
    Sslope(k,i) = P(1); Sint(k) = P(2);
    if sum(size(S.R)) == 4
        err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
        Sslope_err(k,i) = err(1); Sint_err(k) = err(2);
    end
end
%%% Plot trends at specific isotherms
subplot(3,2,i)
    plottime = [1:ceil(length(plotting.time_start)*1)];
for j = 1:length(thermplot)
    idtherm = find(plotting.therm_grid == thermplot(j));
    if length(idtherm) == 1
        h(j) = plot(plotting.time_start, O2gaincorr(idtherm,:),'.','color',C(j,:),'markersize',M); hold on;
%         plot(plotting.time_start, O2int(idtherm) + plotting.time_start*O2slope(idtherm,i),'color',C(j,:),'linewidth',1); hold on;
%         plot(plotting.time_start(plottime), coefa(idtherm,i)*exp((plotting.time_start(plottime) - min(plotting.time_start))*coefb(idtherm,i)) + ...
%             coefc(idtherm,i)*exp((plotting.time_start(plottime) - min(plotting.time_start))*coefd(idtherm,i)),'color',C(j,:),'linewidth',1); hold on;
    end
end
datetick('x',3); ylim([274 304]); %legend(h, thermstr)
title(['O_2 concentration from WFP isotherms below winter ventilation, ' yearspan{i}])
subplot(3,2,i+2)
for j = 1:length(thermplot)
    idtherm = find(plotting.therm_grid == thermplot(j));
    if length(idtherm) == 1
        h(j) = plot(plotting.time_start, plotting.S(idtherm,:),'.','color',C(j,:),'markersize',M); hold on;
        plot(plotting.time_start, Sint(idtherm) + plotting.time_start*Sslope(idtherm,i),'color',C(j,:),'linewidth',1); hold on;
    end
end
datetick('x',3); ylim([34.86 34.93]); %legend(h,thermstr,'location','southwest')
title(['Salinity from WFP isotherms below winter ventilation, ' yearspan{i}])
subplot(3,2,i+4)
for j = 1:length(thermplot)
    idtherm = find(plotting.therm_grid == thermplot(j));
    if length(idtherm) == 1
        h(j) = plot(plotting.time_start, plotting.depth(idtherm,:),'.','color',C(j,:),'markersize',M); hold on;
        %plot(plotting.time_start, Sint(idtherm) + plotting.time_start*Sslope(idtherm,i),'color',C(j,:),'linewidth',1); hold on;
    end
end
datetick('x',3); set(gca,'ydir','reverse'); ylim([1400 2700]); legend(h,thermstr,'location','southeast')
title(['Depth of WFP isotherms below winter ventilation, ' yearspan{i}])
end
%%
figure(2); clf
    Mar = 15; %marker size for plotting
    stable = find(plotting.therm_grid >= min(thermplot) & plotting.therm_grid <= max(thermplot)); %indices of stable deep isotherms
    C = colormap(lbmap(length(tgrid)-1,'RedBlue')); C = colormap(jet);
    C2(1,:) = nicecolor('kwb'); C2(2,:) = nicecolor('kww');
for i = 1:2
subplot(2,2,i)
    ind = find(O2slope(:,i) < 0.1 & O2slope(:,i) > -0.1 & isnan(O2slope_err(:,i)) == 0);
    ind2 = intersect(ind,stable);
    for j = 1:length(ind2)
        %Calculate analytical derivative
        expderiv(j,:) = coefa(ind2(j),i).*coefb(ind2(j),i).*exp((t-tmin).*coefb(ind2(j),i)) + ...
            coefc(ind2(j),i).*coefd(ind2(j),i).*exp((t-tmin).*coefd(ind2(j),i));
        for k = 1:length(tgrid)-1
            expderiv_mean(j,k) = nanmean(expderiv(j,tgrid(k):tgrid(k+1)));
        end
    end
%errorbar(plotting.therm_grid(ind),O2slope(ind,i)*365,O2slope_err(ind,i)*365,'k.'); hold on;
for k = 1:length(tgrid) - 1
    plot(plotting.therm_grid(ind2),expderiv_mean(:,k)*365,'.','color',C(k,:),'markersize',Mar); hold on;
end
for k = 1:length(tgrid2) - 1
    plot(plotting.therm_grid(ind2),O2pieceslope(ind2,k,i)*365,'.','color',C2(k,:),'markersize',Mar); hold on;
end
plot(plotting.therm_grid(ind),O2slope(ind,i)*365,'k.','markersize',Mar-5); hold on;
plot(plotting.therm_grid(ind2),O2slope(ind2,i)*365,'m.','markersize',Mar); hold on;
plot(plotting.therm_grid(ind),zeros(size(ind)),'r--'); hold on;
axis([1.6 4 -30 20]); ylabel('dO_2/dt (\mumol kg^{-1}/yr)'); xlabel('Isotherm') %axis([min(thermplot)-int max(thermplot)+int -20 20]);
title(['Oxygen drift by isotherm for ' yearspan{i}])
text(2, 10, ['mean = ' num2str(nanmean(O2slope(ind2,i)*365),3) ' \mumol kg^{-1}/yr'],'color','m')
wfp_O2drift(i) = nanmean(O2slope(ind2,i));
wfp_O2drift_std(i) = nanstd(O2slope(ind2,i));
wfp_O2drift_num(i) = sum(~isnan(O2slope(ind2,i)));
end
for i = 1:2
subplot(2,2,i+2)
    ind = find(Sslope(:,i) < 0.1 & Sslope(:,i) > -0.1 & isnan(Sslope_err(:,i)) == 0);
    ind2 = intersect(ind,stable);
%errorbar(plotting.therm_grid(ind),Sslope(ind,i),Sslope_err(ind,i),'k.'); hold on;
plot(plotting.therm_grid(ind),Sslope(ind,i)*365,'k.','markersize',Mar-5); hold on;
plot(plotting.therm_grid(ind2),Sslope(ind2,i)*365,'b.','markersize',Mar); hold on;
plot(plotting.therm_grid(ind),zeros(size(ind)),'r--'); hold on;
axis([1.6 4 -0.05 0.035]); ylabel('dS/dt (psu/yr)'); xlabel('Isotherm') %axis([min(thermplot)-int max(thermplot)+int -0.05 0.035]);
title(['Salinity drift by isotherm for ' yearspan{i}])
text(2, 0.02, ['mean = ' num2str(nanmean(Sslope(ind2,i)*365),3) ' psu/yr'],'color','b')
wfp_Sdrift(i) = nanmean(Sslope(ind2,i));
end

%% Find points where WFP and gliders were measuring simultaneously and look at offset to find glider drift rate

toldist = 10; %Glider and WFP must be < toldist apart to make comparison (stick with somewhere in the 5-10 range)
toltime = 0.5; %Glider profile start time must be < toltime from the mean start time of a WFP profile

%Run comparison for GL002 in year 1, GL002 in year 2, and GL003 in year 2
for j = 1:3
if j == 1
    glider = Yr1_GL002_grid;
        %Calculate gain-corrected O2
        glider.O2corr = glider.O2conc.*nanmean(G1);
    wfp = Yr1_wfpgrid; yr = 1;
elseif j > 1
    if j == 2
    glider = Yr2_GL002_grid;
        %Calculate gain-corrected O2
        glider.O2corr = glider.O2conc.*Yr2_gain(1);
    elseif j == 3
    glider = Yr2_GL003_grid;
        %Calculate gain-corrected O2
        glider.O2corr = glider.O2conc.*Yr2_gain(2);
    end
    wfp = Yr2_wfpgrid; yr = 2;       
end
    %WFP essentially stays in place, so just use single location
    latwfp = nanmean(wfp.lat);
    lonwfp = nanmean(wfp.lon);
    %Take mean of two start times because using paired profiles (mostly 0.83 days btwn the two paired profiles)
    wfp.time_start_mean = nanmean([wfp.time_start(wfp.ind_pair) wfp.time_start(wfp.ind_pair + 1)],2);
    %Calculate gain and drift-corrected O2 and drift-corrected S for wfp
    wfp.O2corr = wfp.O2conc*gain_wfp(yr) - ...
            wfp_O2drift(yr)*repmat(wfp.time_start_mean - min(wfp.time_start_mean),1,length(wfp.depth_grid))'; %correct for gain (from Winklers) and then for drift
    wfp.Scorr = wfp.S - wfp_Sdrift(yr)*repmat(wfp.time_start_mean - min(wfp.time_start_mean),1,length(wfp.depth_grid))';


%Pull out intersection of the two depth grids
[depth_grid_overlap,dind_glider,dind_wfp] = intersect(glider.depth_grid, wfp.depth_grid);
%Initialize counter
counter = 0;
clear O2diff Sdiff indcomp

for i = 1:length(glider.profile_ind)
    %Check if it is close enough to the WFP
    dist(i) = distlatlon(latwfp, glider.lat(i), lonwfp, glider.lon(i));
    if dist(i) < toldist
        %Check if there is a WFP profile close enough in time
        [tdiff(i), tind(i)] = min(abs(glider.time_start(i) - wfp.time_start_mean));
        if tdiff(i) < toltime
            %Advance counter and keep track of glider profiles with good comparisons
            counter = counter + 1;
            indcomp(counter) = i;
            %Pull out comparison data at all matching depth intervals
            O2diff(:,counter) = glider.O2corr(dind_glider,i) - wfp.O2corr(dind_wfp,tind(i));
            Sdiff(:,counter) = glider.S(dind_glider,i) - wfp.Scorr(dind_wfp,tind(i));
        end
    end
end

%Visualize comparisons over time
    %Calculate linear fits to data
    x = glider.time_start(indcomp); y1 = nanmean(O2diff)'; y2 = nanmean(Sdiff)';
    ind1 = find(isnan(x) == 0 & isnan(y1) == 0); [P1,S1] = polyfit(x(ind1),y1(ind1),1);
    [yfit1,yerr1] = polyval(P1,x(ind1),S1);
    %Divide salinity into two time sections
    tbreak = datenum(2015,1,15);
    ind2a = find(isnan(x) == 0 & isnan(y2) == 0 & x < tbreak);
    ind2b = find(isnan(x) == 0 & isnan(y2) == 0 & x > tbreak);
    [P2a,S2a] = polyfit(x(ind2a),y2(ind2a),1); [P2b,S2b] = polyfit(x(ind2b),y2(ind2b),1);
    [yfit2a,yerr2a] = polyval(P2a,x(ind2a),S2a);
    [yfit2b,yerr2b] = polyval(P2b,x(ind2b),S2b);
    %Record drift rates for later corrections
    gliderO2drift(j) = P1(1);
    gliderSdrift(j,1) = P2a(1);
    gliderSdrift(j,2) = P2b(1);
    gliderO2uncert(j) = nanmean(yerr1);

figure(2 + j); clf
colormap(cmocean('balance'));
    subplot(411)
imagesc(O2diff); caxis([-30 30]); colorbar
title('O_2 concentration difference, Glider - WFP')
    subplot(412)
imagesc(Sdiff); caxis([-0.15 0.15]); colorbar
title('Salinity difference, Glider - WFP'); %xlabel('Profile number')
    subplot(413)
errorbar(glider.time_start(indcomp),nanmean(O2diff),2*nanstd(O2diff)./sum(~isnan(O2diff)),'k.','markersize',Mar); hold on;
plot(x, (x*P1(1) + P1(2)), 'r-'); hold on;
    plot(x(ind1),yfit1 - yerr1,'r--'); hold on; plot(x(ind1),yfit1 + yerr1,'r--'); hold on;
xlim([min(glider.time_start(indcomp))-5, max(glider.time_start(indcomp))+5]); datetick('x',3,'keeplimits');
title('O_2 concentration difference, Glider - WFP (Mean for full profile)')
    subplot(414)
errorbar(glider.time_start(indcomp),nanmean(Sdiff),2*nanstd(Sdiff)./sum(~isnan(Sdiff)),'k.','markersize',Mar); hold on;
plot(x(ind2a), (x(ind2a)*P2a(1) + P2a(2)), 'r-'); hold on;
    plot(x(ind2a),yfit2a + yerr2a,'r--'); hold on; plot(x(ind2a),yfit2a - yerr2a,'r--'); hold on;
plot(x(ind2b), (x(ind2b)*P2b(1) + P2b(2)), 'r-'); hold on;
    plot(x(ind2b),yfit2b + yerr2b,'r--'); hold on; plot(x(ind2b),yfit2b - yerr2b,'r--'); hold on;
xlim([min(glider.time_start(indcomp))-5, max(glider.time_start(indcomp))+5]); datetick('x',3,'keeplimits');
title('Salinity difference, Glider - WFP (Mean for full profile)')

end

%% Pull together wire-following profiler and deep glider data and plot comparison
%Plot comparison at series of overlapping depths
clear C; C(1,:) = nicecolor('kkw'); C(2,:) = nicecolor('rry');
M = 15;
depthplot = [205 405 605 805];

figure(6); clf
for i = 1:2
    if i == 1
        plotting = Yr1_GL002_grid;
%  Version without having yet made glider drift correction
%         O2gaincorr = plotting.O2conc*nanmean(G1); %note need a better format for saving gain
%         Scorr = plotting.S;
% Version with glider drift correction
        O2gaincorr = plotting.O2conc*nanmean(G1) - ...
            gliderO2drift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
        %Correct salinity in two sections
        ind1 = find(plotting.time_start < datenum(2015,1,15)); ind2 = find(plotting.time_start > datenum(2015,1,15));
        Scorr1 = plotting.S(:,ind1) - gliderSdrift(1,1)*repmat(plotting.time_start(ind1) - min(plotting.time_start(ind1)),1,length(plotting.depth_grid))';
        Scorr2 = plotting.S(:,ind2) - gliderSdrift(1,2)*repmat(plotting.time_start(ind2) - min(plotting.time_start(ind2)),1,length(plotting.depth_grid))'...
            - gliderSdrift(1,1)*(min(plotting.time_start(ind2)) - min(plotting.time_start(ind1)));
        Scorr = [Scorr1 Scorr2];
    elseif i == 2
        plotting = Yr1_wfpgrid;
        plotting.time_start = Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair);
        O2gaincorr = plotting.O2conc*gain_wfp(1) - ...
             wfp_O2drift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
        Scorr = plotting.S - wfp_Sdrift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))';
    end
    for j = 1:length(depthplot)
        iddepth = find(plotting.depth_grid == depthplot(j));
            subplot(length(depthplot),3,j*3)
        plot(plotting.time_start, O2gaincorr(iddepth,:),'.','color',C(i,:),'markersize',M); hold on;
        datetick('x'); ylim([270 315])
        title(['O_2 concentration at ' num2str(depthplot(j)) ' meters'])
        %legend('Open Ocean Glider 002','Wire-following profiler','location','southeast');
        ylabel('\mumol/kg');
             subplot(length(depthplot),3,j*3 - 2)
        plot(plotting.time_start, plotting.T(iddepth,:),'.','color',C(i,:),'markersize',M); hold on;
        datetick('x');
        title(['Temperature at ' num2str(depthplot(j)) ' meters'])
        %legend('Open Ocean Glider 002','Wire-following profiler','location','northeast');
        ylabel('deg C');
              subplot(length(depthplot),3,j*3 - 1)
        plot(plotting.time_start, Scorr(iddepth,:),'.','color',C(i,:),'markersize',M); hold on;
        datetick('x');
        title(['Salinity at ' num2str(depthplot(j)) ' meters'])
        legend('Open Ocean Glider 002','Wire-following profiler','location','northwest');
        %ylabel('deg C');
    end
end

    clear C; C(1,:) = nicecolor('ccbw'); C(2,:) = nicecolor('kkw'); C(3,:) = nicecolor('rry');
figure(7); clf
for i = 1:3
    if i == 1
        plotting = Yr2_GL002_grid;
% Version without having yet made glider drift correction
%         O2gaincorr = plotting.O2conc*Yr2_gain(1); %note need a better format for saving gain
        O2gaincorr = plotting.O2conc*Yr2_gain(1) - ...
            gliderO2drift(2)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
        Scorr = plotting.S;
    elseif i == 2
        plotting = Yr2_GL003_grid;
% Version without having yet made glider drift correction
%         O2gaincorr = plotting.O2conc*Yr2_gain(2); %note need a better format for saving gain
        O2gaincorr = plotting.O2conc*Yr2_gain(2) - ...
            gliderO2drift(3)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
        Scorr = plotting.S;
    elseif i == 3
        plotting = Yr2_wfpgrid;
        plotting.time_start = Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair);
        O2gaincorr = plotting.O2conc*gain_wfp(2) - ...
            wfp_O2drift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
        Scorr = plotting.S - wfp_Sdrift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))';
    end
    for j = 1:length(depthplot)
        iddepth = find(plotting.depth_grid == depthplot(j));
            subplot(length(depthplot),3,j*3)
        plot(plotting.time_start, O2gaincorr(iddepth,:),'.','color',C(i,:),'markersize',M); hold on;
        datetick('x'); xlim([datenum(2015,8,5) datenum(2016,1,1)]);  ylim([270 315])
        title(['O_2 concentration at ' num2str(depthplot(j)) ' meters'])
        %legend('Open Ocean Glider 002','Open Ocean Glider 003','Wire-following profiler','location','southeast');
        ylabel('\mumol/kg');
             subplot(length(depthplot),3,j*3 - 2)
        plot(plotting.time_start, plotting.T(iddepth,:),'.','color',C(i,:),'markersize',M); hold on;
        datetick('x'); xlim([datenum(2015,8,5) datenum(2016,1,1)]);  
        title(['Temperature at ' num2str(depthplot(j)) ' meters'])
        legend('Open Ocean Glider 002','Open Ocean Glider 003','Wire-following profiler','location','northeast');
        ylabel('deg C');
              subplot(length(depthplot),3,j*3 - 1)
        plot(plotting.time_start, plotting.S(iddepth,:),'.','color',C(i,:),'markersize',M); hold on;
        datetick('x'); xlim([datenum(2015,8,5) datenum(2016,1,1)]);  
        title(['Salinity at ' num2str(depthplot(j)) ' meters'])
        %legend('Open Ocean Glider 002','Open Ocean Glider 003','Wire-following profiler','location','northeast');
        ylabel('deg C');
    end
end



%% Make a map of array showing all assets and locations of deployment casts
figure(8); clf
    subplot(2,1,1)
plot(Yr1_GL002_grid.lon,Yr1_GL002_grid.lat,'k','linewidth',1); hold on;
plot([nanmean(Yr1_flmA.lon) nanmean(Yr1_flmB.lon) nanmean(Yr1_sb.lon) nanmean(Yr1_wfpgrid.lon)],...
    [nanmean(Yr1_flmA.lat) nanmean(Yr1_flmB.lat) nanmean(Yr1_sb.lat) nanmean(Yr1_wfpgrid.lat)],'b.','markersize',35);
plot(Yr1_disc.lon360-360,Yr1_disc.lat,'r.','markersize',25);
axis([-40 -39 59.6 60.1]); title('Year 1'); legend('Glider track','Mooring locations','Deployment casts','location','northwest');
    subplot(2,1,2)
plot(Yr2_GL002_grid.lon,Yr2_GL002_grid.lat,'k','linewidth',1); hold on;
plot(Yr2_GL003_grid.lon,Yr2_GL003_grid.lat,'k','linewidth',1); hold on;
plot(Yr2_PG001_grid.lon,Yr2_PG001_grid.lat,'k','linewidth',1); hold on;
plot([nanmean(Yr2_flmA.lon) nanmean(Yr2_flmB.lon) nanmean(Yr2_sb.lon) nanmean(Yr2_wfpgrid.lon)],...
    [nanmean(Yr2_flmA.lat) nanmean(Yr2_flmB.lat) nanmean(Yr2_sb.lat) nanmean(Yr2_wfpgrid.lat)],'b.','markersize',35);
plot(Yr2_disc.lon360-360,Yr2_disc.lat,'r.','markersize',25);
axis([-40 -39 59.6 60.1]); title('Year 2');

%% Calibration of fixed depth assets

time_window = 0.2; %all measurements within X days before or after the cast
    %note that there are no data for the Yr 2 flanking moorings unless this
    %tolerance is widened significantly
tol_dist = 30; %limit of distance between cast and asset in order to use for calibration
tol_surf = 35; %depth limit from which to consider surface data from deployment Winklers

for j = 1:7
if j == 1
    asset = Yr1_flmA; discdata = Yr1_disc;
    time_window = 0.2; tol_dist = 15;
elseif j == 2
    asset = Yr1_flmB; discdata = Yr1_disc;
    time_window = 0.2; tol_dist = 30;
elseif j == 3
    asset = Yr1_sb; discdata = Yr1_disc;
    time_window = 0.2; tol_dist = 15;
elseif j == 4
    asset = Yr1_flmA; discdata = Yr2_disc;
    time_window = 5; tol_dist = 15;
elseif j == 5
    asset = Yr2_flmB; discdata = Yr2_disc;
    time_window = 2; tol_dist = 15;
elseif j == 6
    asset = Yr2_sb; discdata = Yr2_disc;
    time_window = 0.2; tol_dist = 15;
elseif j == 7
    asset = Yr2_rid; discdata = Yr2_disc;
    time_window = 0.2; tol_dist = 15;
end
    %Find surface deployment Winklers near the fixed asset
    surfid = find(discdata.depth < tol_surf);
    for i = 1:length(surfid)
        distance(i) = distlatlon(discdata.lat(surfid(i)), nanmean(asset.lat),...
            discdata.lon360(surfid(i)), nanmean(asset.lon)+360);
    end
    closeid = find(distance < tol_dist);
    %Pull out asset sensor measurements at time of Winklers
    time_Winkler = discdata.day(surfid(closeid)) + discdata.time(surfid(closeid));
    for i = 1:length(closeid)
        id_sensor = find(asset.time_mat > time_Winkler(i) - time_window & asset.time_mat < time_Winkler(i) + time_window);
        gain(i) =  nanmean(discdata.oxy(surfid(closeid(i)))./asset.oxygen(id_sensor));
        gainstd(i) =  nanstd(discdata.oxy(surfid(closeid(i)))./asset.oxygen(id_sensor));
    end
    gain_out(j) = nanmean(gain);
    gain_out_std(j) = nanstd(gain);
    gain_out_numpts(j) = length(gain);
    clear gain distance
end

%% Calculate O2 offset from glider for each asset in Year 1

toldist = 10; %Glider and WFP must be < toldist apart to make comparison (stick with somewhere in the 5-10 range)
toltime = 0.25; %Glider profile start time must be < toltime from the mean start time of a WFP profile

    glider = Yr1_GL002_grid;
        %Use time at end of glider profiles because most are up (so ML measured at end)
        glider.time = glider.time_start + glider.duration;
        %Calculate gain and drift corrected O2 for reference glider
        glider.O2surf = glider.O2conc(3,:)*nanmean(G1) - gliderO2drift(1)*(glider.time_start - min(glider.time_start))';
     %Initialize arrays to hold offset data for each comparison asset
     profilenum = length(glider.profile_ind);
     assetO2 = NaN*ones(profilenum,3); assetO2std = NaN*ones(profilenum,3); assetO2offset = NaN*ones(profilenum,3);
for k = 1:3
    if k == 1
    asset = Yr1_flmA;
        asset.O2corr = asset.oxygen(asset.nightind)*gain_out(1);
    elseif k == 2
    asset = Yr1_flmB;
        asset.O2corr = asset.oxygen(asset.nightind)*gain_out(2);
    elseif k == 3
    asset = Yr1_sb;
        asset.O2corr = asset.oxygen(asset.nightind)*gain_out(3);
    end
    asset.time = asset.time_mat(asset.nightind);
    latasset = nanmean(asset.lat); lonasset = nanmean(asset.lon);

for i = 1:length(glider.profile_ind)
    %Check if glider is close enough to asset to calibrate
    dist(i) = distlatlon(latasset, glider.lat(i), lonasset, glider.lon(i));
    if dist(i) < toldist
        %Find asset measurements close enough in time
        tind = find(abs(asset.time - glider.time(i)) < toltime);
        numt(i) = length(tind);
        %If there are measurements within time range
        if length(tind) > 0
            %Calculate a mean and std of the asset value
            assetO2(i,k) = nanmean(asset.O2corr(tind));
            assetO2std(i,k) = nanstd(asset.O2corr(tind));
            %Calculate offset from glider value
            assetO2offset(i,k) = assetO2(i,k) - glider.O2surf(i);
        end
    end
end

end


figure(9); clf
clear C; C(1,:) = nicecolor('k'); C(2,:) = nicecolor('b'); C(3,:) = nicecolor('rrw');
exp1start = [10 0.5];
for k = 1:3
    ind = find(isnan(assetO2offset(:,k)) == 0);
    [P,S] = polyfit(glider.time(ind),assetO2offset(ind,k),1);
    [yfit,yerr] = polyval(P,glider.time,S);
    hc(k) = plot(glider.time, assetO2offset(:,k), '.', 'color', C(k,:), 'markersize', Mar); hold on;
    plot(glider.time,yfit,'-','color', C(k,:)); hold on;
    text(max(glider.time)-70, 15 - 2.5*k, ['drift = ' num2str(P(1)*365,3) ' \mumol kg^{-1}/yr'],'color',C(k,:))
    Yr1_surfdrift(k) = P(1);
    [Yr1_surfdrift_check(k), Yr1_surfdrift_delta(k)] = polyval(P,min(glider.time) + 365,S);
    [Yr1_surfoffset(k),Yr1_surfoffset_delta(k)] = polyval(P,min(glider.time),S);
    %Calculate and plot exponential fits
%         tmin = min(glider.time(ind)); O2min = min(assetO2offset(ind,k));
%     f = fit(glider.time(ind) - tmin, assetO2offset(ind,k) - O2min, 'exp1','StartPoint',exp1start);
%     coeff = coeffvalues(f);
%     plot(glider.time(ind), coeff(1)*exp((glider.time(ind) - tmin).*coeff(2)) + O2min,'-','color', C(k,:)); hold on;
end
legend(hc,'Flanking Mooring A','Flanking Mooring B','Apex Surface Buoy','location','northwest')
title('Surface mixed layer (top 30 m) assets vs. glider in Year 1, 2014-2015')
ylabel('O_2 concentration offset, Asset - Glider')
datetick('x',3);

%% Calculate O2 offset from glider for each asset in Year 2
    glider = Yr2_GL002_grid;
        %Use time at end of glider profiles because most are up (so ML measured at end)
        glider.time = glider.time_start + glider.duration;
        %Calculate gain and drift corrected O2 for reference glider
        glider.O2surf = glider.O2conc(3,:)*Yr2_gain(1) - gliderO2drift(2)*(glider.time_start - min(glider.time_start))';
     %Initialize arrays to hold offset data for each comparison asset
     profilenum = length(glider.profile_ind);
     assetO2 = NaN*ones(profilenum,3); assetO2std = NaN*ones(profilenum,3); assetO2offset = NaN*ones(profilenum,3);
for k = 1:3
    if k == 1
    asset = Yr2_flmA;
        asset.O2corr = asset.oxygen(asset.nightind)*gain_out(4);
    elseif k == 2
    asset = Yr2_rid;
        asset.O2corr = asset.oxygen(asset.nightind)*gain_out(7);
    elseif k == 3
    asset = Yr2_sb;
        asset.O2corr = asset.oxygen(asset.nightind)*gain_out(6);
    end
    asset.time = asset.time_mat(asset.nightind);
    latasset = nanmean(asset.lat); lonasset = nanmean(asset.lon);

for i = 1:length(glider.profile_ind)
    %Check if glider is close enough to asset to calibrate
    dist(i) = distlatlon(latasset, glider.lat(i), lonasset, glider.lon(i));
    if dist(i) < toldist
        %Find asset measurements close enough in time
        tind = find(abs(asset.time - glider.time(i)) < toltime);
        numt(i) = length(tind);
        %If there are measurements within time range
        if length(tind) > 0
            %Calculate a mean and std of the asset value
            assetO2(i,k) = nanmean(asset.O2corr(tind));
            assetO2std(i,k) = nanstd(asset.O2corr(tind));
            %Calculate offset from glider value
            assetO2offset(i,k) = assetO2(i,k) - glider.O2surf(i);
        end
    end
end

end


figure(10); clf
clear C; C(1,:) = nicecolor('k'); C(2,:) = nicecolor('rmkw'); C(3,:) = nicecolor('rrw');
for k = 1:3
    ind = find(isnan(assetO2offset(:,k)) == 0);
    [P,S] = polyfit(glider.time(ind),assetO2offset(ind,k),1);
    [yfit,yerr] = polyval(P,glider.time,S);
    hc(k) = plot(glider.time, assetO2offset(:,k), '.', 'color', C(k,:), 'markersize', Mar); hold on;
    plot(glider.time,yfit,'-','color', C(k,:)); hold on;
    text(max(glider.time)-40, 19 - 1.5*k, ['drift = ' num2str(P(1)*365,3) ' \mumol kg^{-1}/yr'],'color',C(k,:))
    Yr2_surfdrift(k) = P(1);
    [Yr2_surfdrift_check(k), Yr2_surfdrift_delta(k)] = polyval(P,min(glider.time) + 365,S);
    [Yr2_surfoffset(k),Yr2_surfoffset_delta(k)] = polyval(P,min(glider.time),S);
end
legend(hc,'Flanking Mooring A','Apex Instrument Frame','Apex Surface Buoy','location','northwest')
title('Surface mixed layer (top 30 m) assets vs. glider in Year 2, 2015-2016')
ylabel('O_2 concentration offset, Asset - Glider')
datetick('x',3);

%% Plot time series of ML assets
figure(11); clf
    %subplot(2,1,1)
plot(Yr1_flmA.time_mat(Yr1_flmA.nightind), anomaly(Yr1_flmA.oxygen(Yr1_flmA.nightind)*gain_out(1) -...
    Yr1_surfdrift(1)*(Yr1_flmA.time_mat(Yr1_flmA.nightind) - min(Yr1_flmA.time_mat)) - Yr1_surfoffset(1)), 'k.'); hold on;
plot(Yr1_flmB.time_mat(Yr1_flmB.nightind), anomaly(Yr1_flmB.oxygen(Yr1_flmB.nightind)*gain_out(2) - ...
    Yr1_surfdrift(2)*(Yr1_flmB.time_mat(Yr1_flmB.nightind) - min(Yr1_flmB.time_mat)) - Yr1_surfoffset(2)), 'b.'); hold on;
plot(Yr1_sb.time_mat(Yr1_sb.nightind), anomaly(Yr1_sb.oxygen(Yr1_sb.nightind)*gain_out(3) - ...
    Yr1_surfdrift(3)*(Yr1_sb.time_mat(Yr1_sb.nightind) - min(Yr1_sb.time_mat)) - Yr1_surfoffset(3)), '.','color',nicecolor('rrw')); hold on;
plot(Yr1_GL002_grid.time_start, anomaly(Yr1_GL002_grid.O2conc(3,:)*nanmean(G1) -...
    gliderO2drift(1)*(Yr1_GL002_grid.time_start - min(Yr1_GL002_grid.time_start))'),'c.'); hold on; %25m
%plot(Yr1_sb.time_mat,anomaly(O2sol(Yr1_sb.SSS_dosta, Yr1_sb.SST_dosta)),'g-.','linewidth',3); hold on;
%plot(Yr1_GL002_grid.time_start,anomaly(O2sol(Yr1_GL002_grid.S(3,:), Yr1_GL002_grid.T(3,:))),'g-','linewidth',3); hold on;
% xlim([datenum(2014,9,5) datenum(2015,8,20)]); datetick('x',3,'keeplimits'); ylim([-50 100]); ylabel('\mumol/kg')
% title(['Surface (top 30 m) O_2 anomaly from mean (2014-2015)'])
% legend('Flanking Mooring A','Flanking Mooring B','Apex Surface Buoy','Open Ocean Glider 002','location','northwest');
%     subplot(2,1,2)
plot(Yr2_flmA.time_mat(Yr2_flmA.nightind), anomaly(Yr2_flmA.oxygen(Yr2_flmA.nightind)*gain_out(4) - ...
    Yr2_surfdrift(1)*(Yr2_flmA.time_mat(Yr2_flmA.nightind) - min(Yr2_flmA.time_mat)) - Yr2_surfoffset(1)), 'k.'); hold on;
plot(Yr2_flmB.time_mat(Yr2_flmB.nightind), Yr2_flmB.oxygen(Yr2_flmB.nightind)*gain_out(5), 'b.'); hold on;
plot(Yr2_sb.time_mat(Yr2_sb.nightind), anomaly(Yr2_sb.oxygen(Yr2_sb.nightind)*gain_out(6) - ...
    Yr2_surfdrift(2)*(Yr2_sb.time_mat(Yr2_sb.nightind) - min(Yr2_sb.time_mat)) - Yr2_surfoffset(2)), '.','color',nicecolor('rrw')); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.nightind),anomaly(Yr2_rid.oxygen(Yr2_rid.nightind)*gain_out(7) -...
    Yr2_surfdrift(3)*(Yr2_rid.time_mat(Yr2_rid.nightind) - min(Yr2_rid.time_mat)) - Yr2_surfoffset(3)),'.','color',nicecolor('rmkw')); hold on;
plot(Yr2_GL002_grid.time_start, anomaly(Yr2_GL002_grid.O2conc(3,:)*Yr2_gain(1) -...
    gliderO2drift(2)*(Yr2_GL002_grid.time_start - min(Yr2_GL002_grid.time_start))'),'c.'); hold on; %25m
plot(Yr2_GL003_grid.time_start, anomaly(Yr2_GL003_grid.O2conc(3,:)*Yr2_gain(2) -...
    gliderO2drift(3)*(Yr2_GL003_grid.time_start - min(Yr2_GL003_grid.time_start))'),'.','color',nicecolor('ckw')); hold on;  %25m
plot(Yr2_PG001_grid.time_start, anomaly(Yr2_PG001_grid.O2conc(5,:)*Yr2_gain(3)),'.','color',nicecolor('ck')); hold on;  %25m
%plot(Yr2_rid.time_mat,anomaly(O2sol(Yr2_rid.SSS_dosta, Yr2_rid.SST_dosta)),'g-','linewidth',3); hold on;
% xlim([datenum(2015,8,10) datenum(2016,6,1)]); datetick('x',3,'keeplimits'); ylim([-25 50]); ylabel('\mumol/kg')
% title(['Surface (top 30 m) O_2 anomaly from mean (2015-2016)'])
% legend('Flanking Mooring A','Apex Surface Buoy','Apex Instrument Frame',...
%     'Open Ocean Glider 002','Open Ocean Glider 003','Profiling Glider 001','location','northwest');
xlim([datenum(2014,9,5) datenum(2016,6,1)]); datetick('x',12,'keeplimits'); ylim([-45 70]); ylabel('\mumol/kg')
title(['Surface (top 30 m) O_2 anomaly from mean (2014-2016)'])
%     subplot(2,2,3)
% %No data yet from flmA and flmB because of CTD processing issue
% plot(Yr1_sb.time_mat(Yr1_sb.nightind), Yr1_sb.SST_dosta(Yr1_sb.nightind), 'r.'); hold on;
% plot(Yr1_GL002_grid.time_start, Yr1_GL002_grid.T(3,:),'c.'); hold on; %25m
% datetick('x',3);
%     subplot(2,2,4)
% %No data yet from flmA and flmB because of CTD processing issue
% plot(Yr2_sb.time_mat(Yr2_sb.nightind), Yr2_sb.SST_dosta(Yr2_sb.nightind), 'r.'); hold on;
% plot(Yr2_rid.time_mat(Yr2_rid.nightind),Yr2_rid.SST_dosta(Yr2_rid.nightind),'m.'); hold on;
% plot(Yr2_GL002_grid.time_start, Yr2_GL002_grid.T(3,:),'c.'); hold on; %25m
% plot(Yr2_GL003_grid.time_start, Yr2_GL003_grid.T(3,:),'.','color',nicecolor('ckw')); hold on;  %25m
% plot(Yr2_PG001_grid.time_start, Yr2_PG001_grid.T(5,:),'.','color',nicecolor('ck')); hold on;  %25m
% datetick('x',3);

%% Calculate O2 saturation in surface
A = O2sol(Yr1_sb.SSS_dosta, Yr1_sb.SST_dosta);
B = O2sol(Yr1_GL002_grid.S(3,:), Yr1_GL002_grid.T(3,:));

figure(100); clf
plot(Yr1_sb.time_mat,O2sol(Yr1_sb.SSS_dosta, Yr1_sb.SST_dosta),'g-.','linewidth',3); hold on;
plot(Yr1_GL002_grid.time_start,O2sol(Yr1_GL002_grid.S(3,:), Yr1_GL002_grid.T(3,:)),'g--','linewidth',3); hold on;

plot(Yr2_rid.time_mat,O2sol(Yr2_rid.SSS_dosta, Yr2_rid.SST_dosta),'g--','linewidth',3); hold on;

%% Summarize accuracy and drift rates for all optodes
%Gliders: G1 for year 1 and G2 for year 2
gliders_gain = [nanmean(G1) nanmean(G2)];
gliders_gain_std = [nanstd(G1) nanstd(G2)];
gliders_numcal = [sum(~isnan(G1)) sum(~isnan(G2))];
gliders_sem95 = 2*gliders_gain_std./sqrt(gliders_numcal);
gliders_rawsurfmean = [nanmean(Yr1_GL002_grid.O2conc(3,:)) nanmean(Yr2_GL002_grid.O2conc(3,:)) nanmean(Yr2_GL003_grid.O2conc(3,:)) nanmean(Yr2_PG001_grid.O2conc(5,:))];
gliders_corr = gliders_rawsurfmean.*gliders_gain;
gliders_corr_err = [gliders_rawsurfmean.*(gliders_gain + gliders_sem95) - gliders_rawsurfmean.*(gliders_gain - gliders_sem95)]/2;
gliders_percentaccuracy = gliders_corr_err./gliders_corr*100;

gliderO2drift*365;
gliderO2uncert;

%% WFP
wfp_sem95 = 2*gain_wfpstd./sqrt(gain_wfpnum);
wfp_rawmean = [nanmean(Yr1_wfpgrid.O2conc(21,:)) nanmean(Yr2_wfpgrid.O2conc(21,:))]; %250 m
wfp_corr = wfp_rawmean.*gain_wfp;
wfp_corr_err = [wfp_rawmean.*(gain_wfp + wfp_sem95) - wfp_rawmean.*(gain_wfp - wfp_sem95)]/2;
wfp_percentaccuracy = wfp_corr_err./wfp_corr*100;

wfp_O2drift*365;
wfp_O2drift_std*365;
wfp_O2drift_num;

%% See ML assets above

 %% Looking for stable deep ocean measurement in deep WFP data - plotting on depth surfaces
% figure(5); clf
% depthplot = [1400:150:2400];
% depthstr = {'1400 m', '1550 m', '1700 m', '1850 m', '2000 m', '2150 m', '2300 m'};
% clear C h; top = nicecolor('rrywwwwww'); bot = nicecolor('rrykkkkkk');
% C = [linspace(top(1),bot(1),length(depthplot))' linspace(top(2),bot(2),length(depthplot))' linspace(top(3),bot(3),length(depthplot))'];
% yearspan = {'2014-2015','2015-2016'};
% for i = 1:2
%     if i == 1
%         plotting = Yr1_wfpgrid;
%         plotting.time_start = Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair);
%         O2gaincorr = plotting.O2conc*gain_wfp(1); %note need a better format for saving gain
%         %Initialize arrays to hold slopes
%             O2slope = NaN*ones(length(plotting.depth_grid),2);
%             O2slope_err = NaN*ones(length(plotting.depth_grid),2);
%             Tslope = NaN*ones(length(plotting.depth_grid),2);
%             Tslope_err = NaN*ones(length(plotting.depth_grid),2);
%             Sslope = NaN*ones(length(plotting.depth_grid),2);
%             Sslope_err = NaN*ones(length(plotting.depth_grid),2);
%     elseif i == 2
%         plotting = Yr2_wfpgrid;
%         plotting.time_start = Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair);
%         O2gaincorr = plotting.O2conc*gain_wfp(2); %note need a better format for saving gain
%     end
% %%% Look for trends over time in individual depth intervals
% %(Note that this could also be done on density surfaces)
% for k = 1:length(plotting.depth_grid)
%     ind = find(~isnan(O2gaincorr(k,:)));
%     [P,S] = polyfit(plotting.time_start(ind), O2gaincorr(k,ind)',1);
%     O2slope(k,i) = P(1); O2int(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         O2slope_err(k,i) = err(1); O2int_err(k) = err(2);
%     end
%     [P,S] = polyfit(plotting.time_start(ind), plotting.pdens(k,ind)',1);
%     densslope(k,i) = P(1); densint(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         densslope_err(k,i) = err(1); densint_err(k) = err(2);
%     end
%     [P,S] = polyfit(plotting.time_start(ind), plotting.T(k,ind)',1);
%     Tslope(k,i) = P(1); Tint(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         Tslope_err(k,i) = err(1); Tint_err(k) = err(2);
%     end
%     [P,S] = polyfit(plotting.time_start(ind), plotting.S(k,ind)',1);
%     Sslope(k,i) = P(1); Sint(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         Sslope_err(k,i) = err(1); Sint_err(k) = err(2);
%     end
% end
% %%% Plot trends at specific depths
% subplot(3,2,i)
% for j = 1:length(depthplot)
%     iddepth = find(plotting.depth_grid == depthplot(j));
%     h(j) = plot(plotting.time_start, O2gaincorr(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, O2int(iddepth) + plotting.time_start*O2slope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([274 292]); %legend(h, depthstr)
% title(['O_2 concentration from WFP depths below winter ventilation' yearspan(i)])
% subplot(3,2,i+2)
% for j = 1:length(depthplot)
%     iddepth = find(plotting.depth_grid == depthplot(j));
%     h(j) = plot(plotting.time_start, plotting.T(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, Tint(iddepth) + plotting.time_start*Tslope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([2.3 4]); %legend(h, depthstr)
% title(['Temperature from WFP depths below winter ventilation' yearspan(i)])
% subplot(3,2,i+4)
% for j = 1:length(depthplot)
%     iddepth = find(plotting.depth_grid == depthplot(j));
%     h(j) = plot(plotting.time_start, plotting.S(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, Sint(iddepth) + plotting.time_start*Sslope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([34.86 34.93]); legend(h,depthstr,'location','southwest')
% title(['Salinity from WFP depths below winter ventilation' yearspan(i)])
% % subplot(2,2,i+2)
% % for j = 1:length(depthplot)
% %     iddepth = find(plotting.depth_grid == depthplot(j));
% %     h(j) = plot(plotting.time_start, plotting.pdens(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
% %     plot(plotting.time_start, densint(iddepth) + plotting.time_start*densslope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% % end
% % datetick('x'); ylim([1027.74 1027.86]); %legend(h, depthstr)
% % title(['Density from WFP depths below winter ventilation' yearspan(i)])
% end
% 
% figure(50); clf
% for i = 1:2
% subplot(3,2,i)
%     ind = find(O2slope(:,i) < 0.1 & O2slope(:,i) > -0.1 & isnan(O2slope_err(:,i)) == 0);
% %errorbar(plotting.depth_grid(ind),O2slope(ind,i)*365,O2slope_err(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),O2slope(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),zeros(size(ind)),'r--'); hold on;
% axis([300 2500 -20 20]); ylabel('dO_2/dt (\mumol kg^{-1}/yr)'); xlabel('Depth')
% title(['Oxygen drift by depth for ' yearspan(i)])
% end
% for i = 1:2
% subplot(3,2,i+2)
%     ind = find(Tslope(:,i) < 0.1 & Tslope(:,i) > -0.1 & isnan(Tslope_err(:,i)) == 0);
% %errorbar(plotting.depth_grid(ind),Tslope(ind,i)*365,Tslope_err(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),Tslope(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),zeros(size(ind)),'r--'); hold on;
% axis([300 2500 -0.7 0.35]); ylabel('dT/dt (deg C/yr)'); xlabel('Depth')
% title(['Temperature drift by depth for ' yearspan(i)])
% end
% for i = 1:2
% subplot(3,2,i+4)
%     ind = find(Sslope(:,i) < 0.1 & Sslope(:,i) > -0.1 & isnan(Sslope_err(:,i)) == 0);
% %errorbar(plotting.depth_grid(ind),Sslope(ind,i),Sslope_err(ind,i),'k.'); hold on;
% plot(plotting.depth_grid(ind),Sslope(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),zeros(size(ind)),'r--'); hold on;
% axis([300 2500 -0.07 0.035]); ylabel('dS/dt (psu/yr)'); xlabel('Depth')
% title(['Salinity drift by depth for ' yearspan(i)])
% end


%% Looking for stable deep ocean measurement in deep WFP data - plotting on depth surfaces, now corrected for drift in S and O2
% figure(3); clf
% depthplot = [1400:150:2400];
% depthstr = {'1400 m', '1550 m', '1700 m', '1850 m', '2000 m', '2150 m', '2300 m'};
% clear C h; top = nicecolor('rrywwwwww'); bot = nicecolor('rrykkkkkk');
% C = [linspace(top(1),bot(1),length(depthplot))' linspace(top(2),bot(2),length(depthplot))' linspace(top(3),bot(3),length(depthplot))'];
% yearspan = {'2014-2015','2015-2016'};
% for i = 1:2
%     if i == 1
%         plotting = Yr1_wfpgrid;
%         plotting.time_start = Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair);
%         Scorr = plotting.S - wfp_Sdrift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))';
%         O2gaincorr = plotting.O2conc*gain_wfp(1) - ...
%             wfp_O2drift(1)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
%         %Initialize arrays to hold slopes
%             O2slope = NaN*ones(length(plotting.depth_grid),2);
%             O2slope_err = NaN*ones(length(plotting.depth_grid),2);
%             Tslope = NaN*ones(length(plotting.depth_grid),2);
%             Tslope_err = NaN*ones(length(plotting.depth_grid),2);
%             Sslope = NaN*ones(length(plotting.depth_grid),2);
%             Sslope_err = NaN*ones(length(plotting.depth_grid),2);
%     elseif i == 2
%         plotting = Yr2_wfpgrid;
%         plotting.time_start = Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair);
%         Scorr = plotting.S - wfp_Sdrift(2)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))';
%         O2gaincorr = plotting.O2conc*gain_wfp(2) - ...
%             wfp_O2drift(2)*repmat(plotting.time_start - min(plotting.time_start),1,length(plotting.depth_grid))'; %correct for gain (from Winklers) and then for drift
%     end
% %%% Look for residual trends over time in individual depth intervals
% for k = 1:length(plotting.depth_grid)
%     ind = find(~isnan(O2gaincorr(k,:)));
%     [P,S] = polyfit(plotting.time_start(ind), O2gaincorr(k,ind)',1);
%     O2slope(k,i) = P(1); O2int(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         O2slope_err(k,i) = err(1); O2int_err(k) = err(2);
%     end
%     [P,S] = polyfit(plotting.time_start(ind), plotting.T(k,ind)',1);
%     Tslope(k,i) = P(1); Tint(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         Tslope_err(k,i) = err(1); Tint_err(k) = err(2);
%     end
%     [P,S] = polyfit(plotting.time_start(ind), Scorr(k,ind)',1);
%     Sslope(k,i) = P(1); Sint(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         Sslope_err(k,i) = err(1); Sint_err(k) = err(2);
%     end
% end
% %%% Plot trends at specific depths
% subplot(3,2,i)
% for j = 1:length(depthplot)
%     iddepth = find(plotting.depth_grid == depthplot(j));
%     h(j) = plot(plotting.time_start, O2gaincorr(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, O2int(iddepth) + plotting.time_start*O2slope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([274 292]); %legend(h, depthstr)
% title(['O_2 concentration from WFP depths below winter ventilation' yearspan(i)])
% subplot(3,2,i+2)
% for j = 1:length(depthplot)
%     iddepth = find(plotting.depth_grid == depthplot(j));
%     h(j) = plot(plotting.time_start, plotting.T(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, Tint(iddepth) + plotting.time_start*Tslope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([2.3 4]); %legend(h, depthstr)
% title(['Temperature from WFP depths below winter ventilation' yearspan(i)])
% subplot(3,2,i+4)
% for j = 1:length(depthplot)
%     iddepth = find(plotting.depth_grid == depthplot(j));
%     h(j) = plot(plotting.time_start, Scorr(iddepth,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, Sint(iddepth) + plotting.time_start*Sslope(iddepth,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([34.86 34.93]); legend(h,depthstr,'location','southwest')
% title(['Salinity from WFP depths below winter ventilation' yearspan(i)])
% end
% 
% figure(4); clf
% for i = 1:2
% subplot(3,2,i)
%     ind = find(O2slope(:,i) < 0.1 & O2slope(:,i) > -0.1 & isnan(O2slope_err(:,i)) == 0);
% plot(plotting.depth_grid(ind),O2slope(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),zeros(size(ind)),'r--'); hold on;
% axis([300 2500 -20 20]); ylabel('dO_2/dt (\mumol kg^{-1}/yr)'); xlabel('Depth')
% title(['Oxygen drift by depth for ' yearspan(i)])
% end
% for i = 1:2
% subplot(3,2,i+2)
%     ind = find(Tslope(:,i) < 0.1 & Tslope(:,i) > -0.1 & isnan(Tslope_err(:,i)) == 0);
% plot(plotting.depth_grid(ind),Tslope(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),zeros(size(ind)),'r--'); hold on;
% axis([300 2500 -0.7 0.35]); ylabel('dT/dt (deg C/yr)'); xlabel('Depth')
% title(['Temperature drift by depth for ' yearspan(i)])
% end
% for i = 1:2
% subplot(3,2,i+4)
%     ind = find(Sslope(:,i) < 0.1 & Sslope(:,i) > -0.1 & isnan(Sslope_err(:,i)) == 0);
% plot(plotting.depth_grid(ind),Sslope(ind,i)*365,'k.'); hold on;
% plot(plotting.depth_grid(ind),zeros(size(ind)),'r--'); hold on;
% axis([300 2500 -0.07 0.035]); ylabel('dS/dt (psu/yr)'); xlabel('Depth')
% title(['Salinity drift by depth for ' yearspan(i)])
% end

%% Plotting on density surfaces - note that this doesn't work until drift in salinity is corrected
% figure(6); clf
% densplot = [1027.76:0.01:1027.86];
% %densstr = {'27.76','27.78','27.80','27.82','27.84','27.86'};
% clear C; top = nicecolor('rrywwwwww'); bot = nicecolor('rrykkkkkk');
% C = [linspace(top(1),bot(1),length(densplot))' linspace(top(2),bot(2),length(densplot))' linspace(top(3),bot(3),length(densplot))'];
% yearspan = {'2014-2015','2015-2016'};
% for i = 1:2
%     if i == 1
%         plotting = Yr1_wfpgrid_dens;
%         plotting.time_start = Yr1_wfpgrid_dens.time_start(Yr1_wfpgrid_dens.ind_pair);
%         O2gaincorr = plotting.O2conc*gain_wfp(1); %note need a better format for saving gain
%         %Initialize arrays to hold slopes
%         O2slope = NaN*ones(length(plotting.dens_grid),2);
%         O2slope_err = NaN*ones(length(plotting.dens_grid),2);
%         Tslope = NaN*ones(length(plotting.dens_grid),2);
%         Tslope_err = NaN*ones(length(plotting.dens_grid),2);
%     elseif i == 2
%         plotting = Yr2_wfpgrid_dens;
%         plotting.time_start = Yr2_wfpgrid_dens.time_start(Yr2_wfpgrid_dens.ind_pair);
%         O2gaincorr = plotting.O2conc*gain_wfp(2); %note need a better format for saving gain
%     end
% %%% Look for trends over time in individual density intervals
% %(Note that this could also be done on density surfaces)
% for k = 1:length(plotting.dens_grid)
%     ind = find(~isnan(O2gaincorr(k,:)));
%     [P,S] = polyfit(plotting.time_start(ind), O2gaincorr(k,ind)',1);
%     O2slope(k,i) = P(1); O2int(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         O2slope_err(k,i) = err(1); O2int_err(k) = err(2);
%     end
%     [P,S] = polyfit(plotting.time_start(ind), plotting.T(k,ind)',1);
%     Tslope(k,i) = P(1); Tint(k) = P(2);
%     if sum(size(S.R)) == 4
%         err = sqrt(diag((S.R)\inv(S.R'))./S.normr.^2./S.df);
%         Tslope_err(k,i) = err(1); Tint_err(k) = err(2);
%     end
% end
% %%% Plot trends at specific densities
% subplot(2,2,i)
% for j = 1:length(densplot)
%     iddens = find(plotting.dens_grid == densplot(j));
%     h(j) = plot(plotting.time_start, O2gaincorr(iddens,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, O2int(iddens) + plotting.time_start*O2slope(iddens,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([274 292]); %legend(h, densstr)
% title(['O_2 concentration from WFP density surfaces below winter ventilation' yearspan(i)])
% subplot(2,2,i+2)
% for j = 1:length(densplot)
%     iddens = find(plotting.dens_grid == densplot(j));
%     h(j) = plot(plotting.time_start, plotting.T(iddens,:),'.','color',C(j,:),'markersize',M); hold on;
%     plot(plotting.time_start, Tint(iddens) + plotting.time_start*Tslope(iddens,i),'color',C(j,:),'linewidth',1); hold on;
% end
% datetick('x'); ylim([2.3 4]); %legend(h, densstr)
% title(['Temperature from WFP density surfaces below winter ventilation' yearspan(i)])
% end
% 
% 
% figure(7); clf
% for i = 1:2
% subplot(2,1,i)
%     ind = find(O2slope(:,i) < 0.1 & O2slope(:,i) > -0.1 & isnan(O2slope_err(:,i)) == 0);
% errorbar(plotting.dens_grid(ind),O2slope(ind,1),O2slope_err(ind,1),'k.'); hold on;
% end



