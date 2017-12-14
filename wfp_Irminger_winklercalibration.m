%% Load discrete oxygen data from Irminger deployment cruises
loadWinklerIrminger

%% Plot Year 1 discrete casts at Irminger OOI site, and line up with corresponding gider profiles
C(1,:) = nicecolor('rrry'); C(2,:) = nicecolor('gggkb'); C(3,:) = nicecolor('bbc'); C(4,:) = nicecolor('rbm'); %color
time_window = 1; %all profiles within X days before or after the cast
depth_window = 20; %compare points in profile only if within X meters of Winkler sampling depth
tol_O2std = 5; %tolerance - only calculate gain if std of profile O2 measurements at a given depth are within X
tol_dist = 30; %limit of distance between cast and glider in order to use for calibration

for k = 1:2
    if k == 1
        cast_list = [5,6,7,9]; %casts in OOI site region (no test casts)
        wfp_plot = Yr1_wfpgrid;
        disc_plot = Yr1_disc;
    elseif k == 2
        cast_list = [6,9,10,11,12,13]; %casts in OOI site region (no test casts, and not 4 and 5 b/c before gliders in water) %casts 11,12
        wfp_plot = Yr2_wfpgrid;
        disc_plot = Yr2_disc;
    end
figure(10 + k); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=28;
    height=20;
    set(gcf,'units','centimeters','position',[x0,y0,width,height]) 
G = NaN*ones(length(disc_plot.cast),1);
    for i = 1:length(cast_list)
        subplot(2,length(cast_list)/2,i)
        ind = find(disc_plot.cast == cast_list(i));
            A = disc_plot.time(ind(1)) + disc_plot.day(ind(1));
            ind_time = find(wfp_plot.time_start < A + time_window & wfp_plot.time_start > A - time_window);
            dist = distlatlon(wfp_plot.lat(ind_time(1)),disc_plot.lat(ind(end)),wfp_plot.lon(ind_time(1))+360,disc_plot.lon360(ind(end)));
        plot(disc_plot.oxy(ind),disc_plot.depth(ind),'k*','color',C(1,:),'linewidth',4); hold on;
        plot(wfp_plot.O2conc(:,ind_time),wfp_plot.depth_grid,'k.'); hold on;
        set(gca,'YDir','reverse'); axis([240 315 0 2500]); ylabel('Depth (m)'); xlabel('O_2 concentration')
        title({['Cast ' num2str(cast_list(i)) ', Distance to mooring = ' num2str(dist,3) ' km']});
        for j = 1:length(ind) %loop over all discrete samples in cast
            [depth_dif,depth_ind] = min(abs(wfp_plot.depth_grid - disc_plot.depth(ind(j))));
            if depth_dif < depth_window %only use points within X meters of sampling depth
                O2raw = nanmean(wfp_plot.O2conc(depth_ind,ind_time));
                O2raw_std = nanstd(wfp_plot.O2conc(depth_ind,ind_time));
                if O2raw_std < tol_O2std & dist < tol_dist %don't calculate gain if the profiles at this depth/time are too variable
                    G(ind(j)) = disc_plot.oxy(ind(j))./O2raw;
                    plot(disc_plot.oxy(ind(j)),disc_plot.depth(ind(j)),'bo','markersize',10,'linewidth',2); hold on;
                end
            end
        end
    end   

% Calculate gain correction
gain(k) = nanmean(G);
gainstd(k) = nanstd(G);
gainnum(k) = sum(~isnan(G));

% Replot WFP data, now making gain correction
    for i = 1:length(cast_list)
        subplot(2,length(cast_list)/2,i)
        ind = find(disc_plot.cast == cast_list(i));
            A = disc_plot.time(ind(1)) + disc_plot.day(ind(1));
            ind_time = find(wfp_plot.time_start < A + time_window & wfp_plot.time_start > A - time_window);
            dist = distlatlon(wfp_plot.lat(ind_time(1)),disc_plot.lat(ind(end)),wfp_plot.lon(ind_time(1))+360,disc_plot.lon360(ind(end)));
        plot(wfp_plot.O2conc(:,ind_time)*nanmean(G),wfp_plot.depth_grid,'.','color',nicecolor('kw')); hold on;
        plot(wfp_plot.O2conc(:,ind_time)*(nanmean(G)-nanstd(G)),wfp_plot.depth_grid,'.','color',nicecolor('kww')); hold on;
        plot(wfp_plot.O2conc(:,ind_time)*(nanmean(G)+nanstd(G)),wfp_plot.depth_grid,'.','color',nicecolor('kww')); hold on;
        plot(disc_plot.oxy(ind),disc_plot.depth(ind),'k*','color',C(1,:),'linewidth',4); hold on;
    end  
    
end
    