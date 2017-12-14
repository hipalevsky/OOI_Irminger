%% Calculate O2 consumption rates at depths below stratified mixed layer

wfpmerge.time = [Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair); Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair)];
wfpmerge.depth_grid = Yr1_wfpgrid.depth_grid;
wfpmerge.T = [Yr1_wfpgrid.T Yr2_wfpgrid.T];
wfpmerge.S = [Yr1_wfpgrid.S - repmat(wfp_Sdrift(1)*(Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair) - min(Yr1_wfpgrid.time_start)),1,length(wfpmerge.depth_grid))'...
    Yr2_wfpgrid.S - repmat(wfp_Sdrift(2)*(Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair) - min(Yr2_wfpgrid.time_start)),1,length(wfpmerge.depth_grid))'];
wfpmerge.O2 = [Yr1_wfpgrid.O2conc*gain_wfp(1) - wfp_O2drift(1)*repmat(Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair) - min(Yr1_wfpgrid.time_start),1,length(wfpmerge.depth_grid))'...
    Yr2_wfpgrid.O2conc*gain_wfp(2) - wfp_O2drift(2)*repmat(Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair) - min(Yr2_wfpgrid.time_start),1,length(wfpmerge.depth_grid))'];

%Visualize sections
% figure(100); clf
% imagesc(wfpmerge.T); colormap(cmocean('thermal'))
% figure(101); clf
% imagesc(wfpmerge.S); colormap(cmocean('haline'))
% figure(102); clf
% imagesc(wfpmerge.O2); colormap(cmocean('thermal')); caxis([260 330]); colorbar
%%
%Pull out some depth levels to plot as lines through time series
depth_id = [11:20:151]; %intervals more like {3:4:251] once finalized
windowSize = 18; 
b = (1/windowSize)*ones(1,windowSize);
windowSize2 = 12;
b2 = (1/windowSize2)*ones(1,windowSize2);
windowSize3 = 20;
b3 = (1/windowSize3)*ones(1,windowSize3);
a = 1;
delay = nanmean(diff(wfpmerge.time))*windowSize/2;
delay2 = nanmean(diff(wfpmerge.time))*(windowSize*windowSize2)/2;
tol = 30; %window of variability to keep above/below mean (assume rest are outliers)
zerotol = [0.05 -0.01]; %tolerance for determining if derivative is different enough from zero (separate for up vs down)
yrbreak = datenum(2015,10,30);
ventstart = NaN*ones(length(depth_id),2); ventend = NaN*ones(length(depth_id),2);
ventstartO2 = NaN*ones(length(depth_id),2); ventendO2 = NaN*ones(length(depth_id),2);

% figure(103); clf
% for i = 1:length(depth_id)
%     subplot(ceil(length(depth_id)/2),2,i)
%     plot(wfpmerge.time, wfpmerge.T(depth_id(i),:),'k.'); hold on;
%         y = filter(b,a,wfpmerge.T(depth_id(i),:));
%         plot(wfpmerge.time, y,'b'); hold on;
%     title(['T at ' num2str(wfpmerge.depth_grid(depth_id(i))) ' m'])
%     datetick('x'); ylim([3.2 5.5])
% end
% 
% figure(104); clf
% for i = 1:length(depth_id)
%     subplot(ceil(length(depth_id)/2),2,i)
%     plot(wfpmerge.time, wfpmerge.S(depth_id(i),:),'k.'); hold on;
%     title(['S at ' num2str(wfpmerge.depth_grid(depth_id(i))) ' m'])
%     datetick('x',3); ylim([34.84 35]);
% end

figure(105); clf
for i = 1:length(depth_id)
clear indventstart indventend
    subplot(ceil(length(depth_id)/2),2,i)
    %Plot O2 over time
    plot(wfpmerge.time, wfpmerge.O2(depth_id(i),:),'k.'); hold on;
        meanval = nanmean(wfpmerge.O2(depth_id(i),:));
        ind = find(isnan(wfpmerge.O2(depth_id(i),:)) == 0 & wfpmerge.O2(depth_id(i),:) > meanval - tol & wfpmerge.O2(depth_id(i),:) < meanval + tol);
        y = filter(b,a,wfpmerge.O2(depth_id(i),ind));
        x = wfpmerge.time(ind);
        plot(x(windowSize:end) - delay, y(windowSize:end),'b','linewidth',2); hold on;
    title(['O2 at ' num2str(wfpmerge.depth_grid(depth_id(i))) ' m'])
    xlim([min(wfpmerge.time) max(wfpmerge.time)]); datetick('x',3,'keepticks'); ylim([280 320]);
    %Calculate and plot O2 saturation
%     O2sat = O2sol(wfpmerge.S(depth_id(i),:),wfpmerge.T(depth_id(i),:));
%     O2satfilt = filter(b,a,O2sat(ind));
%     plot(x(windowSize:end) - delay, O2satfilt(windowSize:end),'color',nicecolor('gkw'),'linewidth',2); hold on;
    %Plot derivatives
    %subplot(ceil(length(depth_id)),2,i + length(depth_id))
        dy = diff(y); dyfilt = filter(b2,a,dy);
    plot(x(windowSize+1:end) - delay, dyfilt(windowSize:end)*16+290,'m','linewidth',2); hold on;
    %plot(x(windowSize+1:end) - delay, dy(windowSize:end)*16+290,'r','linewidth',0.5); hold on;
    plot(x, 290*ones(size(x)),'k--','linewidth',1); hold on;
    %Look at sign changes in derivative
    dysign = sign(dyfilt);
    ind0 = find(dyfilt > 0 & dyfilt < zerotol(1)); ind0 = find(dyfilt < 0 & dyfilt > zerotol(2));
    dysign(ind0) = 0; %set to zero where within small tolerance of zero
    indpos = find(dysign == 1); indneg = find(dysign == -1);
    %plot(x(indpos) - delay, dyfilt(indpos)*16+290,'c.','linewidth',2); hold on;
    %plot(x(indneg) - delay, dyfilt(indneg)*16+290,'y.','linewidth',2); hold on;
    %Calculate ventilation beginning and ending dates
        int = 1; %number of points to use in identifying ventilation switch
        month = str2num(datestr(x,5));
        ventind = zeros(size(x));
        for k = windowSize:length(dyfilt)-int
            if sum(dysign(k-int:k) < 0 & sum(dysign(k:k+int)) > 0 & sum(ventind(k-int:k-1)) == 0) %Changes sign and hasn't already been identified
                if month(k)>=11 | month(k)<= 3
                    ventind(k) = 1; %where ventilation starts
                end
            elseif sum(dysign(k-int:k) > 0 & sum(dysign(k:k+int)) < 0 & sum(ventind(k-int:k-1)) == 0) %Changes sign and hasn't already been identified
                if month(k) >= 4 & month(k) <= 8
                    ventind(k) = -1; %where ventilation ends
                end
            end
        end
        indventstart = find(ventind == 1);
        indventend = find(ventind == -1);
        %Pull out single value for each year for start/end
            idyr1 = find(x < yrbreak); idyr2 = find(x > yrbreak);
        if length(intersect(indventstart,idyr1)) >= 1
            ventstart(i,1) = x(min(intersect(indventstart,idyr1))) - delay;
            ventstartO2(i,1) = y(min(intersect(indventstart,idyr1)));
        end
        if length(intersect(indventstart,idyr2)) >= 1
            ventstart(i,2) = x(min(intersect(indventstart,idyr2))) - delay;
            ventstartO2(i,2) = y(min(intersect(indventstart,idyr2)));
        end
        if length(intersect(indventend,idyr1)) >= 1
            ventend(i,1) = x(min(intersect(indventend,idyr1))) - delay;
            ventendO2(i,1) = y(min(intersect(indventend,idyr1)));
            ventend_ts(i,:) = x(min(intersect(indventend,idyr1)):min(intersect(indventend,idyr1)) + 150) - delay;
            ventend_tsO2(i,:) = y(min(intersect(indventend,idyr1)):min(intersect(indventend,idyr1)) + 150);
        end
        if length(intersect(indventend,idyr2)) >= 1
            ventend(i,2) = x(min(intersect(indventend,idyr2))) - delay; 
            ventendO2(i,2) = y(min(intersect(indventend,idyr2)));
        end
         %plot(x(indventstart) - delay, dyfilt(indventstart)*16+290,'r*','linewidth',2); hold on;
         %plot(x(indventend) - delay, dyfilt(indventend)*16+290,'g*','linewidth',2); hold on;       
        plot(ventstart(i,:), [290 290],'*','color',nicecolor('rk'),'linewidth',2); hold on;
        plot(ventend(i,:), [290 290],'*','color',nicecolor('gk'),'linewidth',2); hold on;
    title(['O2 change at ' num2str(wfpmerge.depth_grid(depth_id(i))) ' m'])
    xlim([min(wfpmerge.time) max(wfpmerge.time)]); datetick('x',3,'keepticks'); 
end
%%
%Respiration rates by depth intervals
figure(126); clf
O2change = (ventendO2(:,1) - ventstartO2(:,2));
tResp = (ventstart(:,2) - ventend(:,1));
plot(O2change,wfpmerge.depth_grid(depth_id),'k.','MarkerSize',20); %maybe should not keep ./tResp*365 (because don't want per year, but per season?)
set(gca,'ydir','reverse')
title('Respiration over 2015 stratified season')
ylabel('Depth (m)');
xlabel('O_2 decrease (\mumol kg^{-1})')
axis([0 18 150 950])

%Calculate respiration rates for monthly time intervals from ventend_ts and
%ventend_tsO2
    monlist = [str2num(datestr(ventend(:,1),5)) str2num(datestr(ventstart(:,2),5))];
    respRateMon = NaN*ones(length(depth_id),12);
for i = 1:length(depth_id)
    mon_ts = str2num(datestr(ventend_ts(i,:),5)); %months of full time series
    for j = min(monlist(:,1)):max(monlist(:,2))
        idmon = find(mon_ts == j);
        dt = ventend_ts(i,max(idmon)) - ventend_ts(i,min(idmon));
        if length(dt) == 1
            respRateMon(i,j) = -(ventend_tsO2(i,max(idmon)) - ventend_tsO2(i,min(idmon)))/dt*365; %umol/kg/yr
        end
    end
end

figure(107); clf
C = colormap(lbmap(length(depth_id),'RedBlue')); C = colormap(jet);
subplot(211)
for i = 1:length(depth_id)
    plot([4:10],respRateMon(i,4:10),'.-','color',C(i,:),'markersize',30,'linewidth',1); hold on;
end
subplot(212)
plot([4:10],nanmean(respRateMon(:,4:10)),'k.','markersize',30); hold on;

%% Rough estimate of seasonal respiration rate
thermocline_resp = nansum(O2change*20)/1000; %Needs to be updated b/c NaNs in O2change