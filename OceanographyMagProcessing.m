%% Data pipeline to create final figures for Oceanography mag paper
%Run code to read in all relevant data
glider_Irminger
mooring_extract_Irminger
IrmingerAssetStitching
IrmingerDeepRespiration

%% Close existing figures before making these new ones
close all

%% Plot surface mixed layer properties
figure(1); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=25;
    height=17;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])   
mindate = datenum(2014,9,05); maxdate = datenum(2016,11,25);
F = 12; %title font size

%Paired color scheme from colorbrewer
C_flmA = [31 120 180]/255;
C_flmB = [166 206 227]/255; 
C_sb = [202 178 214]/255; 
C_rid = [106 61 154]/255;
C_G02 = [51 160 44]/255;
C_G03 = [178 223 138]/255;
surfticks = [datenum(2014,10,1) datenum(2015,1,1) datenum(2015,4,1) datenum(2015,7,1) datenum(2015,10,1) datenum(2016,1,1) datenum(2016,4,1) datenum(2016,7,1)];
surfticks2 = [datenum(2014,[10:12],1) datenum(2015,[1:12],1) datenum(2016,[1:7],1)];
    L = 2.5;

%%% Create and plot cleaned SST and SSS data
Yr1_flmA.indML = find(Yr1_flmA.time_mat_ctd30 < datenum(2015,7,15)); 
Yr2_flmA.indML = find(Yr2_flmA.time_mat_ctd30 > datenum(2015,10,15) & Yr2_flmA.time_mat_ctd30 < datenum(2016,5,10));
Yr1_sb.indgood = find(Yr1_sb.time_mat_fl < datenum(2015,2,15));
Yr2_rid.indgood1 = find(Yr2_rid.time_mat < datenum(2016,1,25));
Yr2_rid.indgood2 = find(Yr2_rid.time_mat > (datenum(2016,6,8)));
    filteringpoints = 200;
Yr1_sb.SST_fl_filt = medfilt1(Yr1_sb.SST_fl(Yr1_sb.indgood),filteringpoints);
Yr2_rid.SST_dosta_filt1 = medfilt1(Yr2_rid.SST_dosta(Yr2_rid.indgood1),filteringpoints);
Yr2_rid.SST_dosta_filt2 = medfilt1(Yr2_rid.SST_dosta(Yr2_rid.indgood2),filteringpoints);
Yr1_sb.SSS_fl_filt = medfilt1(Yr1_sb.SSS_fl(Yr1_sb.indgood),filteringpoints);
Yr2_rid.SSS_dosta_filt1 = medfilt1(Yr2_rid.SSS_dosta(Yr2_rid.indgood1),filteringpoints);
Yr2_rid.SSS_dosta_filt2 = medfilt1(Yr2_rid.SSS_dosta(Yr2_rid.indgood2),filteringpoints);

subplot(311)
plot(Yr1_flmA.time_mat_ctd30(Yr1_flmA.indML(3:end)),Yr1_flmA.T_ctd30(Yr1_flmA.indML(3:end)),'-','color',C_flmA,'linewidth',L); hold on;
plot(Yr1_sb.time_mat_fl(Yr1_sb.indgood(2:end)),Yr1_sb.SST_fl_filt(2:end),'-','color',C_sb,'linewidth',L); hold on;
plot(Yr2_flmA.time_mat_ctd30(Yr2_flmA.indML),Yr2_flmA.T_ctd30(Yr2_flmA.indML),'-','color',C_flmA,'linewidth',L); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.indgood1(2:end)),Yr2_rid.SST_dosta_filt1(2:end),'-','color',C_sb,'linewidth',L); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.indgood2(2:end)),Yr2_rid.SST_dosta_filt2(2:end),'-','color',C_flmA,'linewidth',L); hold on;
xlim([mindate maxdate])
xticks(surfticks2);
xticklabels({'10/01/14' '' '' '01/01/15' '' '' '04/01/15' '' '' '07/01/15' '' '' '10/01/15' '' '' '01/01/16' '' '' '04/01/16' '' '' '07/01/16'})
%datetick('x',2,'keeplimits','keepticks');
%legend('Flanking Mooring A','Apex Mooring','location','southeast')
ylabel('^oC','Fontsize',10)
title('Sea surface temperature','fontsize',F)
%     subplot(212)
% plot(Yr1_flmA.time_mat_ctd30(Yr1_flmA.indML(3:end)),Yr1_flmA.S_ctd30(Yr1_flmA.indML(3:end)),'r.'); hold on;
% plot(Yr1_sb.time_mat_fl(Yr1_sb.indgood(2:end)),Yr1_sb.SSS_fl_filt(2:end),'k.'); hold on;
% plot(Yr2_flmA.time_mat_ctd30(Yr2_flmA.indML),Yr2_flmA.S_ctd30(Yr2_flmA.indML),'r.'); hold on;
% plot(Yr2_rid.time_mat(Yr2_rid.indgood1(2:end)),Yr2_rid.SSS_dosta_filt1(2:end),'k.'); hold on;
% plot(Yr2_rid.time_mat(Yr2_rid.indgood2(2:end)),Yr2_rid.SSS_dosta_filt2(2:end),'k.'); hold on;
% xlim([datenum(2014,9,1) datenum(2016,8,15)])
% ylim([34.5 35.5])
% datetick('x',12,'keeplimits');
% legend('Flanking Mooring A (~30 m)','Apex mooring MET (~1 m)','location','northwest')
% ylabel('Salinity','Fontsize',10)
% title('Sea surface salinity, Years 1-2 OOI Irminger Sea')

%%% Create time merged SST and SSS product
surfaceST = [Yr1_flmA.time_mat_ctd30(Yr1_flmA.indML(3:end)); Yr1_sb.time_mat_fl(Yr1_sb.indgood(2:end));...
    Yr2_flmA.time_mat_ctd30(Yr2_flmA.indML); Yr2_rid.time_mat(Yr2_rid.indgood1(2:end)); Yr2_rid.time_mat(Yr2_rid.indgood2(2:end))];
surfaceST(:,2) = [Yr1_flmA.T_ctd30(Yr1_flmA.indML(3:end)); Yr1_sb.SST_fl_filt(2:end);...
    Yr2_flmA.T_ctd30(Yr2_flmA.indML); Yr2_rid.SST_dosta_filt1(2:end); Yr2_rid.SST_dosta_filt2(2:end)];
surfaceST(:,3) = [Yr1_flmA.S_ctd30(Yr1_flmA.indML(3:end)); Yr1_sb.SSS_fl_filt(2:end);...
    Yr2_flmA.S_ctd30(Yr2_flmA.indML); Yr2_rid.SSS_dosta_filt1(2:end); Yr2_rid.SSS_dosta_filt2(2:end)];
%1st column is sorted times, second column is SST, third column is SSS
surfaceSTsort = sortrows(surfaceST,1);
O2sat = O2sol(surfaceSTsort(:,3),surfaceSTsort(:,2));

%Filter O2 saturation for smoother plotting
O2sat_time = surfaceSTsort(:,1);
windowSize = 50; 
b = (1/windowSize)*ones(1,windowSize); a = 1;
O2sat_smooth = medfilt1(O2sat,windowSize);

%%% Surface chlorophyll data
    Yr1_flmA.nightind_fl_ML = intersect(Yr1_flmA.nightind_fl,Yr1_flmA.indML);
    Yr2_flmA.nightind_fl_ML = intersect(Yr2_flmA.nightind_fl,Yr2_flmA.indML);
        Yr2_rid.indgoodchl = find(Yr2_rid.time_mat_fl < datenum(2016,2,5) | Yr2_rid.time_mat_fl > datenum(2016,4,1))
    Yr2_rid.nightind_fl_good = intersect(Yr2_rid.nightind_fl,Yr2_rid.indgoodchl);
    filteringpointschl = 12;
    Yr2_flmA.chlfilt = medfilt1(Yr2_flmA.chl(Yr2_flmA.nightind_fl_ML),filteringpointschl);
    Yr2_rid.chlfilt = medfilt1(Yr2_rid.chl(Yr2_rid.nightind_fl_good),filteringpointschl);
    Yr1_flmA.chlfilt = medfilt1(Yr1_flmA.chl(Yr1_flmA.nightind_fl_ML),filteringpointschl);
    Yr2_rid.timebreak = max(find(Yr2_rid.time_mat_fl(Yr2_rid.nightind_fl_good) < datenum(2016,3,1)));
    Yr2_flmA.timebreak = max(find(Yr2_flmA.time_mat_fl(Yr2_flmA.nightind_fl_ML) < datenum(2016,3,1)));

subplot(312)
plot(Yr1_flmA.time_mat_fl(Yr1_flmA.nightind_fl_ML(2:end)),Yr1_flmA.chlfilt(2:end),'-','color',C_flmA,'linewidth',L); hold on;
plot(Yr2_rid.time_mat_fl(Yr2_rid.nightind_fl_good(2:Yr2_rid.timebreak)),Yr2_rid.chlfilt(2:Yr2_rid.timebreak),'-','color',C_rid,'linewidth',L); hold on;
plot(Yr2_rid.time_mat_fl(Yr2_rid.nightind_fl_good(Yr2_rid.timebreak+1:end)),Yr2_rid.chlfilt(Yr2_rid.timebreak+1:end),'-','color',C_rid,'linewidth',L); hold on;
plot(Yr2_flmA.time_mat_fl(Yr2_flmA.nightind_fl_ML(2:Yr2_flmA.timebreak)),Yr2_flmA.chlfilt(2:Yr2_flmA.timebreak),'-','color',C_flmA,'linewidth',L); hold on;
plot(Yr2_flmA.time_mat_fl(Yr2_flmA.nightind_fl_ML(Yr2_flmA.timebreak+1:end)),Yr2_flmA.chlfilt(Yr2_flmA.timebreak+1:end),'-','color',C_flmA,'linewidth',L); hold on;
xlim([mindate maxdate])
ylim([-0.1 18])
% xticks(surfticks);
% datetick('x',2,'keeplimits','keepticks');
xticks(surfticks2);
xticklabels({'10/01/14' '' '' '01/01/15' '' '' '04/01/15' '' '' '07/01/15' '' '' '10/01/15' '' '' '01/01/16' '' '' '04/01/16' '' '' '07/01/16'})
%legend('Flanking Mooring A','Apex Mooring','location','northeast')
ylabel('\mug L^{-1}','Fontsize',10)
title('Mixed layer chlorophyll-a concentration','fontsize',F)

%%% Surface oxygen compilation
filteringpointsO2 = 20; %remove spikes in data

Yr1_flmA.indMLO2 = find(Yr1_flmA.time_mat < datenum(2015,7,15)); 
Yr1_flmA.O2ind = intersect(Yr1_flmA.nightind,Yr1_flmA.indMLO2);
Yr1_flmA.O2cleaned = medfilt1(Yr1_flmA.oxygen(Yr1_flmA.O2ind)*gain_out(1) -...
    Yr1_surfdrift(1)*(Yr1_flmA.time_mat(Yr1_flmA.O2ind) - min(Yr1_flmA.time_mat)) - Yr1_surfoffset(1),filteringpointsO2/2);

Yr2_flmA.indMLO2 = find(Yr2_flmA.time_mat > datenum(2015,10,15) & Yr2_flmA.time_mat < datenum(2016,5,10));
Yr2_flmA.O2ind = intersect(Yr2_flmA.nightind,Yr2_flmA.indMLO2);
Yr2_flmA.O2cleaned = medfilt1(Yr2_flmA.oxygen(Yr2_flmA.O2ind)*gain_out(4) - ...
    Yr2_surfdrift(1)*(Yr2_flmA.time_mat(Yr2_flmA.O2ind) - min(Yr2_flmA.time_mat)) - Yr2_surfoffset(1),filteringpointsO2);

Yr1_flmB.indMLO2 = find(Yr1_flmB.time_mat < datenum(2015,7,15)); 
Yr1_flmB.O2ind = intersect(Yr1_flmB.nightind,Yr1_flmB.indMLO2);
Yr1_flmB.O2cleaned = Yr1_flmB.oxygen(Yr1_flmB.O2ind)*gain_out(2) - ...
    Yr1_surfdrift(2)*(Yr1_flmB.time_mat(Yr1_flmB.O2ind) - min(Yr1_flmB.time_mat)) - Yr1_surfoffset(2);
%Note: not including flmB in Year2 because only a tiny amount of data (part of Nov15)

Yr1_sb.O2cleaned = medfilt1((Yr1_sb.oxygen(Yr1_sb.nightind)*gain_out(3) - ...
    Yr1_surfdrift(3)*(Yr1_sb.time_mat(Yr1_sb.nightind) - min(Yr1_sb.time_mat)) - Yr1_surfoffset(3)),filteringpointsO2);

Yr2_sb.idgood = find(Yr2_sb.time_mat < datenum(2016,1,25));
Yr2_sb.O2ind = intersect(Yr2_sb.idgood,Yr2_sb.nightind);
Yr2_sb.O2cleaned = medfilt1((Yr2_sb.oxygen(Yr2_sb.O2ind)*gain_out(6) - ...
    Yr2_surfdrift(2)*(Yr2_sb.time_mat(Yr2_sb.O2ind) - min(Yr2_sb.time_mat)) - Yr2_surfoffset(2)),filteringpointsO2);

Yr2_rid.isgood = find(Yr2_rid.time_mat < datenum(2016,1,25) | Yr2_rid.time_mat > datenum(2016,4,1));
Yr2_rid.O2ind = intersect(Yr2_rid.isgood,Yr2_rid.nightind);
Yr2_rid.O2cleaned = Yr2_rid.oxygen(Yr2_rid.O2ind)*gain_out(7) -...
    Yr2_surfdrift(3)*(Yr2_rid.time_mat(Yr2_rid.O2ind) - min(Yr2_rid.time_mat)) - Yr2_surfoffset(3);

Yr1_GL002_grid.O2cleaned = Yr1_GL002_grid.O2conc(3,:)*nanmean(G1) -...
    gliderO2drift(1)*(Yr1_GL002_grid.time_start - min(Yr1_GL002_grid.time_start))';
Yr2_GL002_grid.O2cleaned = Yr2_GL002_grid.O2conc(3,:)*Yr2_gain(1) -...
    gliderO2drift(2)*(Yr2_GL002_grid.time_start - min(Yr2_GL002_grid.time_start))';
Yr2_GL003_grid.O2cleaned = Yr2_GL003_grid.O2conc(3,:)*Yr2_gain(2) -...
    gliderO2drift(3)*(Yr2_GL003_grid.time_start - min(Yr2_GL003_grid.time_start))';

%Note: not including profiling glider b/c can't correct for drift by
%comparing w/ WFP b/c not enough depth overlap

%WFP in winter to validate absolute O2 conc
Yr1_wfpgrid.time_winind = find(Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair) > datenum(2015,1,1) & Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair) < datenum(2015,3,1));
Yr1_wfpgrid.time_win = Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair(Yr1_wfpgrid.time_winind));
Yr1_wfpgrid.O2corr_MLwin = Yr1_wfpgrid.O2conc(11,Yr1_wfpgrid.time_winind).*gain_wfp(1) -...
    wfp_O2drift(1)*(Yr1_wfpgrid.time_start(Yr1_wfpgrid.ind_pair(Yr1_wfpgrid.time_winind)) - min(Yr1_wfpgrid.time_start))';
Yr2_wfpgrid.time_winind = find(Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair) > datenum(2015,12,1) & Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair) < datenum(2016,2,1));
Yr2_wfpgrid.time_win = Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair(Yr2_wfpgrid.time_winind));
Yr2_wfpgrid.O2corr_MLwin = Yr2_wfpgrid.O2conc(11,Yr2_wfpgrid.time_winind).*gain_wfp(2) -...
    wfp_O2drift(2)*(Yr2_wfpgrid.time_start(Yr2_wfpgrid.ind_pair(Yr2_wfpgrid.time_winind)) - min(Yr2_wfpgrid.time_start))';

O2sat_timebreak = [max(find(O2sat_time < datenum(2015,8,1))) max(find(O2sat_time < datenum(2016,6,1)))];
Yr2_rid.timebreak = [max(find(Yr2_rid.time_mat(Yr2_rid.O2ind) < datenum(2016,3,1)))];
Yr2_flmA.timebreak = [max(find(Yr2_flmA.time_mat(Yr2_flmA.O2ind) < datenum(2016,3,1)))];

subplot(313)
plot(O2sat_time(2:O2sat_timebreak(1)),O2sat_smooth(2:O2sat_timebreak(1)),'k-','linewidth',L/2); hold on;

plot(Yr1_flmA.time_mat(Yr1_flmA.O2ind(2:end)), (Yr1_flmA.O2cleaned(2:end)), '-','color',C_flmA,'linewidth',L); hold on;
plot(Yr1_flmB.time_mat(Yr1_flmB.O2ind(2:end)), (Yr1_flmB.O2cleaned(2:end)), '-','color',C_flmB,'linewidth',L); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.O2ind(2:Yr2_rid.timebreak)),(Yr2_rid.O2cleaned(2:Yr2_rid.timebreak)),'-','color',C_rid,'linewidth',L); hold on; %For legend
plot(Yr1_sb.time_mat(Yr1_sb.nightind(2:end)), (Yr1_sb.O2cleaned(2:end)), '-','color',C_sb,'linewidth',L); hold on;
plot(Yr1_GL002_grid.time_start, (Yr1_GL002_grid.O2cleaned),'-','color',C_G02,'linewidth',L); hold on; %25m
plot(Yr2_GL003_grid.time_start, (Yr2_GL003_grid.O2cleaned),'-','color',C_G03,'linewidth',L); hold on;  %For legend
%plot(Yr1_wfpgrid.time_win,Yr1_wfpgrid.O2corr_MLwin,'b.'); hold on;

plot(O2sat_time(O2sat_timebreak(1)+2:O2sat_timebreak(2)),O2sat_smooth(O2sat_timebreak(1)+2:O2sat_timebreak(2)),'k-','linewidth',L/2); hold on;
plot(O2sat_time(O2sat_timebreak(2)+2:end),O2sat_smooth(O2sat_timebreak(2)+2:end),'k-','linewidth',L/2); hold on;

plot(Yr2_flmA.time_mat(Yr2_flmA.O2ind(2:Yr2_flmA.timebreak)), (Yr2_flmA.O2cleaned(2:Yr2_flmA.timebreak)), '-','color',C_flmA,'linewidth',L); hold on;
plot(Yr2_flmA.time_mat(Yr2_flmA.O2ind(Yr2_flmA.timebreak+2:end)), (Yr2_flmA.O2cleaned(Yr2_flmA.timebreak+2:end)), '-','color',C_flmA,'linewidth',L); hold on;
plot(Yr2_sb.time_mat(Yr2_sb.O2ind(2:end)), (Yr2_sb.O2cleaned(2:end)), '-','color',C_sb,'linewidth',L); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.O2ind(2:Yr2_rid.timebreak)),(Yr2_rid.O2cleaned(2:Yr2_rid.timebreak)),'-','color',C_rid,'linewidth',L); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.O2ind(Yr2_rid.timebreak+1:end)),(Yr2_rid.O2cleaned(Yr2_rid.timebreak+1:end)),'-','color',C_rid,'linewidth',L); hold on;
plot(Yr2_GL002_grid.time_start, (Yr2_GL002_grid.O2cleaned),'-','color',C_G02,'linewidth',L); hold on; %25m
plot(Yr2_GL003_grid.time_start, (Yr2_GL003_grid.O2cleaned),'-','color',C_G03,'linewidth',L); hold on;  %25m
%plot(Yr2_wfpgrid.time_win,Yr2_wfpgrid.O2corr_MLwin,'b.'); hold on;


legend('Equilibrium O_2','Flanking Mooring A','Flanking Mooring B','Apex Mooring (~12m)','Apex Mooring (~1m)','Glider 002','Glider 003',...
    'location','east');
xlim([mindate maxdate]);
% xticks(surfticks);
% datetick('x',2,'keeplimits','keepticks');
xticks(surfticks2);
xticklabels({'10/01/14' '' '' '01/01/15' '' '' '04/01/15' '' '' '07/01/15' '' '' '10/01/15' '' '' '01/01/16' '' '' '04/01/16' '' '' '07/01/16'})
ylim([260 410]); ylabel('\mumol kg^{-1}')
title(['Mixed layer dissolved oxygen concentration'],'fontsize',F)

%% Calculate O2 supersaturation at maximum
[O2max_flmA,I] = max(Yr1_flmA.O2cleaned)
maxdate_flmA = datestr(Yr1_flmA.time_mat(Yr1_flmA.O2ind(I)))
densmax_flmA = sw_dens0(Yr1_flmA.S_ctd30(Yr1_flmA.O2ind(I)),Yr1_flmA.S_ctd30(Yr1_flmA.O2ind(I)));
[O2max_flmB,I] = max(Yr1_flmB.O2cleaned)
maxdate_flmB = datestr(Yr1_flmB.time_mat(Yr1_flmB.O2ind(I)))

O2satmaxO2_date = find(O2sat_time > datenum(2015,5,10) & O2sat_time < datenum(2015,5,12));
O2satmaxO2 = nanmean(O2sat_smooth(O2satmaxO2_date))
flmA_Aprilbaseline_date = find(Yr1_flmA.time_mat(Yr1_flmA.O2ind) > datenum(2015,4,1) & Yr1_flmA.time_mat(Yr1_flmA.O2ind) < datenum(2015,4,10));
flmA_AprilbaselineO2 = nanmean(Yr1_flmA.O2cleaned(flmA_Aprilbaseline_date))
flmB_Aprilbaseline_date = find(Yr1_flmB.time_mat(Yr1_flmB.O2ind) > datenum(2015,4,1) & Yr1_flmB.time_mat(Yr1_flmB.O2ind) < datenum(2015,4,10));
flmB_AprilbaselineO2 = nanmean(Yr1_flmB.O2cleaned(flmB_Aprilbaseline_date))
%For back of the envelope calculation of O2 inventory above saturation,
%need MLD

%% Calculate O2 supersaturation over bloom period (2 weeks before to 2 weeks after O2 maximum)
bloom_beg = datenum(maxdate_flmA) - 10;
bloom_end = datenum(maxdate_flmA) + 20;

bloom_dates_flmA = find(Yr1_flmA.time_mat(Yr1_flmA.O2ind) > bloom_beg & Yr1_flmA.time_mat(Yr1_flmA.O2ind) < bloom_end);
bloom_dates_flmB = find(Yr1_flmA.time_mat(Yr1_flmB.O2ind) > bloom_beg & Yr1_flmB.time_mat(Yr1_flmB.O2ind) < bloom_end);
bloom_dates_ST = find(surfaceSTsort(:,1) > bloom_beg & surfaceSTsort(:,1) < bloom_end);

O2_bloommean_flmA = nanmean(Yr1_flmA.O2cleaned(bloom_dates_flmA));
O2_bloommean_flmB = nanmean(Yr1_flmB.O2cleaned(bloom_dates_flmB));
O2_bloommean_equil = O2sol(nanmean(surfaceSTsort(bloom_dates_ST,3)),nanmean(surfaceSTsort(bloom_dates_ST,2)));


%% Calculate filtered SST (to determine minimum value and date)

windowSize = 100;
b = (1/windowSize)*ones(1,windowSize);
a = 1;
delay = nanmean(diff(surfaceSTsort(:,1)))*windowSize/2;

SST_filt = filter(b,a,surfaceSTsort(:,2));
SST_filt_date = surfaceSTsort(:,1) - delay;
    ind1 = find(SST_filt_date > (min(surfaceSTsort(:,1) + 5)) & SST_filt_date < datenum(2015,7,1));
[SSTmin_Yr1,I] = min(SST_filt(ind1));
SSTmindate_Yr1 = datestr(SST_filt_date(ind1(I)));
    ind2 = find(SST_filt_date > datenum(2015,7,1));
[SSTmin_Yr2,I] = min(SST_filt(ind2));
SSTmindate_Yr2 = datestr(SST_filt_date(ind2(I)));
[SSTmax,I] = max(SST_filt(ind2));
    SSTmaxdate = datestr(SST_filt_date(ind2(I)));

%% Determine mixed layer depth during time of May 2015 spring bloom/O2 peak
datemin = datenum(2015,5,1);
datemax = datenum(2015,6,1);
windowSize = 100;
b = (1/windowSize)*ones(1,windowSize);
a = 1;


figure(5); clf
    ind = find(Yr1_flmA.d28_time_mat > datemin & Yr1_flmA.d28_time_mat < datemax);
plot(Yr1_flmA.d28_time_mat(ind), filter(b,a,Yr1_flmA.d28_T(ind)), '.','color',nicecolor('ccbw')); hold on;
    ind = find(Yr1_flmA.d48_time_mat > datemin & Yr1_flmA.d48_time_mat < datemax);
plot(Yr1_flmA.d48_time_mat(ind), filter(b,a,Yr1_flmA.d48_T(ind)), '.','color',nicecolor('cbwk')); hold on;
    ind = find(Yr1_flmA.d78_time_mat > datemin & Yr1_flmA.d78_time_mat < datemax);
plot(Yr1_flmA.d78_time_mat(ind), filter(b,a,Yr1_flmA.d78_T(ind)), '.','color',nicecolor('cbbbwkk')); hold on;
    ind = find(Yr1_flmA.d119_time_mat > datemin & Yr1_flmA.d119_time_mat < datemax);
plot(Yr1_flmA.d119_time_mat(ind), filter(b,a,Yr1_flmA.d119_T(ind)), '.','color',nicecolor('bbk')); hold on;
    ind = find(Yr1_flmA.d170_time_mat > datemin & Yr1_flmA.d170_time_mat < datemax);
plot(Yr1_flmA.d170_time_mat(ind), filter(b,a,Yr1_flmA.d170_T(ind)), '.','color',nicecolor('bkk')); hold on;
    ind = find(SST_filt_date > datemin & SST_filt_date < datemax);
plot(SST_filt_date(ind), SST_filt(ind), 'k.'); hold on;
plot(SST_filt_date(ind), SST_filt(ind) - 0.2, 'k--','linewidth',2); hold on;
axis([datemin + 5 datemax 3.6 5])
datetick('x','keepticks')
%%% Based on this, determine that MLD is ~30 m at O2 maximum (ephemeral
%%% stratification)

%%% Density at max
O2peak = (O2max_flmA - flmA_AprilbaselineO2)*densmax_flmA.*30/(10^6); %O2 change in umol/kg* density in kg/m3* MLD in m *(1 mol/10^6 umol)

%% Plot WFP corrected O2 data and overlain MLD

%Use dataset combined over both years, with corrected O2 (from IrmingerDeepRespiration)

wfpmerge.pdens = sw_pden(wfpmerge.S, wfpmerge.T, repmat(sw_pres(wfpmerge.depth_grid,60)',1,length(wfpmerge.time)), zeros(size(wfpmerge.S)));
wfpmerge.ptemp = sw_ptmp(wfpmerge.S, wfpmerge.T, repmat(sw_pres(wfpmerge.depth_grid,60)',1,length(wfpmerge.time)), zeros(size(wfpmerge.S)));
wfpmerge.O2sat = O2sol(wfpmerge.S,wfpmerge.T);
wfpmerge.DO2 = (wfpmerge.O2./wfpmerge.O2sat - 1)*100;

%%% Calculate MLD within WFP data
wfpmerge.mld = NaN*ones(length(wfpmerge.time),1); %initialize array to hold mld data
wfpmerge.mld_T = NaN*ones(length(wfpmerge.time),1); %initialize array to hold mld data
wfpmerge.SST = NaN*ones(length(wfpmerge.time),1); %initialize array to hold SST data
wfpmerge.SSS = NaN*ones(length(wfpmerge.time),1); %initialize array to hold SSS data
mld_id = NaN*ones(length(wfpmerge.time),1);
    tol = 0.5; %use surface properties from x time after start of profile
    criterion_dens = 0.03;
    criterion_T = 0.2;

wfpmerge.ML_sdT = NaN*ones(length(wfpmerge.time),1); %initialize array to hold stdev of temp in ML
wfpmerge.ML_sdS = NaN*ones(length(wfpmerge.time),1); %initialize array to hold stdev of salinity in ML
wfpmerge.ML_sdDens = NaN*ones(length(wfpmerge.time),1); %initialize array to hold stdev of density in ML

%Calculate surface properties to compare with WFP data
for i = 1:length(wfpmerge.time)
    id_surf = find(surfaceSTsort(:,1) > wfpmerge.time(i) & surfaceSTsort(:,1) < (wfpmerge.time(i) + tol));
    wfpmerge.SST(i) = nanmean(surfaceSTsort(id_surf,2));
    wfpmerge.SSS(i) = nanmean(surfaceSTsort(id_surf,3));
end
wfpmerge.surfdens = sw_dens0(wfpmerge.SSS, wfpmerge.SST);

%Calculate MLD from WFP data using density criterion
for k = 1:length(wfpmerge.time)
    if ~isnan(wfpmerge.surfdens(k)) == 1
       idtherm = min(find(wfpmerge.pdens(:,k) > (wfpmerge.surfdens(k) + criterion_dens))); %find depth index for first depth interval below MLD
       idML = find(wfpmerge.pdens(:,k) < (wfpmerge.surfdens(k) + criterion_dens)); %find depths within ML
       if length(idtherm) == 1 & length(idML) > 0
            wfpmerge.mld(k) = mean(wfpmerge.depth_grid(idtherm-1:idtherm)); %set MLD midway between the two depth intervals bracketing the MLD bottom
       else
           wfpmerge.mld(k) = NaN;
       end
    end
end

%Calculate MLD from WFP data using temperature criterion
for k = 1:length(wfpmerge.time)
    if ~isnan(wfpmerge.SST(k)) == 1
       idtherm = min(find(wfpmerge.T(:,k) < (wfpmerge.SST(k) - criterion_T))); %find depth index for first depth interval below MLD
       idML = find(wfpmerge.T(:,k) > (wfpmerge.SST(k) - criterion_T)); %find depths within ML
       if length(idtherm) == 1 & length(idML) > 0
            wfpmerge.mld_T(k) = mean(wfpmerge.depth_grid(idtherm-1:idtherm)); %set MLD midway between the two depth intervals bracketing the MLD bottom
            wfpmerge.ML_sdT(k) = nanstd(wfpmerge.T(k,idML));
            wfpmerge.ML_sdS(k) = nanstd(wfpmerge.S(k,idML));
            wfpmerge.ML_sdDens(k) = nanstd(wfpmerge.pdens(k,idML));
%             if nanstd(wfpmerge.O2(k,idML)) > 2
%                 wfpmerge.mld_T(k) = NaN;
%             end
       else
           wfpmerge.mld_T(k) = NaN;
       end
    end
end

%% Calculate date of first ventilation (i.e. 1st ML to reach given depth in a year)
maxdepth = 1100; %don't include ventilation calculated below this depth

for k = 1:2
    if k == 1
        ids = find(wfpmerge.time < datenum(2015,7,1));
    elseif k == 2
        ids = find(wfpmerge.time > datenum(2015,7,1));
    end
A = [wfpmerge.time(ids) wfpmerge.mld_T(ids)]; A = sortrows(A,1);
    idnonan = ~isnan(A(:,2));
A = A(idnonan,:);
clear firstvent
firstvent(1,:) = A(1,:); ventind = 1;
for i = 2:length(A)
    if A(i,2) > firstvent(ventind,2)
        firstvent(ventind+1,:) = A(i,:);
        ventind = ventind + 1;
    end
end
    idkeep = find(firstvent(:,2) < maxdepth);
    if k == 1
        wfpmerge.firstventYr1 = firstvent(idkeep,:);
    elseif k == 2
        wfpmerge.firstventYr2 = firstvent(idkeep,:);
    end
end

%%% Interpolate dates of first ventilation onto even depth grid
wfpmerge.ventdateYr1_grid = interp1(wfpmerge.firstventYr1(:,2),wfpmerge.firstventYr1(:,1),wfpmerge.depth_grid);
wfpmerge.ventdateYr2_grid = interp1(wfpmerge.firstventYr2(:,2),wfpmerge.firstventYr2(:,1),wfpmerge.depth_grid);

    
%% Calculate respiration rates for seasonal thermocline
%Pull out some depth levels to plot as lines through time series
depth_id = [1:171];
depth_id_plot = [11:60:161];
windowSize = 10; 
b = (1/windowSize)*ones(1,windowSize);
a = 1;
delay = nanmean(diff(wfpmerge.time))*windowSize/2;
tol = 30; %window of variability to keep above/below mean (assume rest are outliers)
zerotol = [0.05 -0.01]; %tolerance for determining if derivative is different enough from zero (separate for up vs down)
yrbreak = datenum(2015,10,30);
ventstart = NaN*ones(length(depth_id),2); ventend = NaN*ones(length(depth_id),2);
ventstartO2 = NaN*ones(length(depth_id),2); ventendO2 = NaN*ones(length(depth_id),2);
minval = 280; %assume any values below this are erroneous
counter = 0; %initialize counter for plotting
    C_vent = nicecolor('rrr');
    C_endvent = nicecolor('yyyw');

figure(3); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=22;
    height=5;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])   
for i = 1:length(depth_id)
clear indventstart indventend
    %Filter data through time series
        meanval = nanmean(wfpmerge.O2(depth_id(i),:));
        indval = find(isnan(wfpmerge.O2(depth_id(i),:)) == 0 & wfpmerge.O2(depth_id(i),:) > meanval - tol &...
            wfpmerge.O2(depth_id(i),:) < meanval + tol & wfpmerge.O2(depth_id(i),:) > minval);
        indtime = find(wfpmerge.time < datenum(2015,8,15) | wfpmerge.time > datenum(2015,8,25));
        ind = intersect(indval,indtime);
        y = filter(b,a,wfpmerge.O2(depth_id(i),ind));
        x = wfpmerge.time(ind);
    %Determine date and O2 value of maximum and minimum
        indfullyr = find(x > wfpmerge.ventdateYr1_grid(depth_id(i)) & x < datenum(2015,7,1));
        if length(indfullyr) > 1
            [O2max(i),I] = max(y(indfullyr));
            O2maxdate(i) = x(indfullyr(I)) - delay;
        else
            O2max(i) = NaN; O2maxdate(i) = NaN;
        end
        indpostmax = find(x > datenum(2015,7,1) & x < wfpmerge.ventdateYr2_grid(depth_id(i)));
        if length(indpostmax) > 1
            [O2min(i),I] = min(y(indpostmax));
            O2mindate(i) = x(indpostmax(I)) - delay;
        else
            O2min(i) = NaN; O2mindate(i) = NaN;
        end
        indsecondyr = find(x > datenum(2016,2,1));
        if length(indsecondyr) > 2
            [O2max2(i),I] = max(y(indsecondyr));
            O2max2date(i) = x(indsecondyr(I)) - delay;
        else
             O2max2(i) = NaN; O2max2date(i) = NaN;
        end
     %Plot for subset of depths
     if ismember(depth_id(i),depth_id_plot) == 1
         counter = counter + 1;
            subplot(1,ceil(length(depth_id_plot)),counter)
        plot(wfpmerge.time(ind), wfpmerge.O2(depth_id(i),ind),'k.'); hold on;
        plot(x(windowSize:end) - delay, y(windowSize:end),'b','linewidth',2); hold on;
        plot([wfpmerge.ventdateYr1_grid(depth_id(i)) wfpmerge.ventdateYr1_grid(depth_id(i))],[280 320],'--','color',C_vent,'linewidth',1); hold on;
        plot([wfpmerge.ventdateYr2_grid(depth_id(i)) wfpmerge.ventdateYr2_grid(depth_id(i))],[280 320],'--','color',C_vent,'linewidth',1); hold on;
            plot(O2maxdate(i), O2max(i),'bo','markerfacecolor',C_endvent,'markersize',7); hold on;
            plot(O2mindate(i), O2min(i),'bo','markerfacecolor',nicecolor('ccw'),'markersize',7); hold on;
            %plot(O2max2date(i), O2max2(i),'m.','markersize',15); hold on;
        title(['Oxygen at ' num2str(wfpmerge.depth_grid(depth_id(i))) ' m'])
        datetick('x',12,'keepticks'); xlim([min(wfpmerge.time) max(wfpmerge.time)]);
        if counter == 1
            ylabel('\mumol kg^{-1}')
        end
     end
end

%% Plotting oxygen sections from wire-following profiler
figure(2); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=22;
    height=14;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])   
    F = 12; %title font size
    C1 = cmocean('Thermal'); C2 = cmocean('Dense'); C3 = cmocean('Balance'); %C = colormap(flipud(lbmap(101,'RedBlue')));
    mindepth = 0; maxdepth = 1650; %maxdepth = 2450;
    cints = 60; %number of contour intervals
    L = 1.75;
    
    subplot(311)
        cmin = 3.3; cmax = 4.5; %cmin = 2; cmax = 5;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
    [X,Y] = meshgrid(wfpmerge.time,Yr1_wfpgrid.depth_grid);
    contourf(X,Y,squeeze(wfpmerge.ptemp),cvec,'linecolor','none'); hold on;
        axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
    plot(wfpmerge.firstventYr1(:,1), wfpmerge.firstventYr1(:,2),'-','color',C_vent,'linewidth',L); hold on;
    plot(wfpmerge.firstventYr2(:,1), wfpmerge.firstventYr2(:,2),'-','color',C_vent,'linewidth',L); hold on;
    plot(O2maxdate, wfpmerge.depth_grid(depth_id), '-','color',C_endvent,'linewidth',L); hold on;
    colormap(C2); hcb = colorbar; set(gca,'YDir','reverse'); ylabel('Depth (m)');
    title(['Potential temperature'],'Fontsize',F); title(hcb,'^oC')
    %datetick('x',2,'keeplimits'); 
    xticks(surfticks2);
    xticklabels({'10/01/14' '' '' '01/01/15' '' '' '04/01/15' '' '' '07/01/15' '' '' '10/01/15' '' '' '01/01/16' '' '' '04/01/16' '' '' '07/01/16'})

%     subplot(311)
%         cmin = 27.5; cmax = 28;
%     cvec = [cmin:(cmax-cmin)/cints:cmax];
%     [X,Y] = meshgrid(wfpmerge.time,Yr1_wfpgrid.depth_grid);
%     contourf(X,Y,squeeze(wfpmerge.pdens-1000),cvec,'linecolor','none'); hold on;
%         axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
%     plot(wfpmerge.firstventYr1(:,1), wfpmerge.firstventYr1(:,2), 'r-','linewidth',L); hold on;
%     plot(wfpmerge.firstventYr2(:,1), wfpmerge.firstventYr2(:,2), 'r-','linewidth',L); hold on;
%     plot(O2maxdate, wfpmerge.depth_grid(depth_id), 'm-','linewidth',L); hold on;
%     colormap(C2); hcb = colorbar; set(gca,'YDir','reverse'); datetick('x',2,'keeplimits'); ylabel('Depth (m)');
%     title(['Potential density'],'Fontsize',F); title(hcb,'\sigma_\Theta')    
    
    subplot(312)
        cmin = 275; cmax = 310; %cmin = 260; cmax = 315;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
    contourf(X,Y,wfpmerge.O2,cvec,'linecolor','none'); hold on;
        axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
    %plot(wfpmerge.time, wfpmerge.mld_T,'g.'); hold on;
    plot(wfpmerge.firstventYr1(:,1), wfpmerge.firstventYr1(:,2),'-','color',C_vent,'linewidth',L); hold on;
    plot(wfpmerge.firstventYr2(:,1), wfpmerge.firstventYr2(:,2),'-','color',C_vent,'linewidth',L); hold on;
    plot(O2maxdate, wfpmerge.depth_grid(depth_id),'-','color',C_endvent,'linewidth',L); hold on;
    colormap(C2); hcb = colorbar; set(gca,'YDir','reverse'); ylabel('Depth (m)');
    %datetick('x',2,'keeplimits'); 
    xticks(surfticks2);
    xticklabels({'10/01/14' '' '' '01/01/15' '' '' '04/01/15' '' '' '07/01/15' '' '' '10/01/15' '' '' '01/01/16' '' '' '04/01/16' '' '' '07/01/16'})
    title(['Dissolved oxygen concentration'],'Fontsize',F); title(hcb,'\mumol kg^{-1}')
    
    subplot(313)
        cmin = -13; cmax = 0; %cmin = -14; cmax = 0;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
    contourf(X,Y,wfpmerge.DO2,cvec,'linecolor','none'); hold on;
        axis([min(wfpmerge.time) max(wfpmerge.time) mindepth maxdepth]); caxis([cmin cmax]);
    %plot(wfpmerge.time, wfpmerge.mld_T,'g.'); hold on;
    plot(wfpmerge.firstventYr1(:,1), wfpmerge.firstventYr1(:,2),'-','color',C_vent,'linewidth',L); hold on;
    plot(wfpmerge.firstventYr2(:,1), wfpmerge.firstventYr2(:,2),'-','color',C_vent,'linewidth',L); hold on;
    plot(O2maxdate, wfpmerge.depth_grid(depth_id),'-','color',C_endvent,'linewidth',L); hold on;
    %plot(O2max2date, wfpmerge.depth_grid(depth_id), 'm.','linewidth',L); hold on;
    colormap(C2); hcb = colorbar; set(gca,'YDir','reverse'); ylabel('Depth (m)');
    title(['Dissolved oxygen saturation state (\DeltaO_2)'],'Fontsize',F); title(hcb,'%')
    %datetick('x',2,'keeplimits');
    xticks(surfticks2);
    xticklabels({'10/01/14' '' '' '01/01/15' '' '' '04/01/15' '' '' '07/01/15' '' '' '10/01/15' '' '' '01/01/16' '' '' '04/01/16' '' '' '07/01/16'})

    
%% Plot respiration rates

figure(4); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=10;
    height=8;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])  
    M = 8;
    %subplot(131)
plot(O2max - O2min, wfpmerge.depth_grid(depth_id), 'k.','markersize',M); hold on;
set(gca,'YDir','reverse');
title('Total stratified season respiration'); xlabel('O_2 decrease (\mumol kg^{-1})');
ylabel('Depth (m)'); ylim([180 1000])

%     subplot(132)
% plot(O2mindate - O2maxdate, wfpmerge.depth_grid(depth_id), 'k.'); hold on;
% set(gca,'YDir','reverse');
% title('Length of stratified season'); xlabel('Days')
%     subplot(133)
% plot((O2max - O2min)./(O2mindate - O2maxdate)*365, wfpmerge.depth_grid(depth_id), 'k.'); hold on;
% set(gca,'YDir','reverse');
% title('Respiration rate'); xlabel('\mumol O_2 kg^{-1} yr^{-1}')

    meandens = nanmean(wfpmerge.pdens');
ThermResp = (O2max - O2min).*meandens(depth_id)*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000000; %respiration per depth interval (mol/m2)
    intdepth = 1000; %choose integration depth for base of seasonal thermocline (top is 200 m because that's where data starts)
    id_basetherm = find(wfpmerge.depth_grid == intdepth);
ThermResp_Int = nansum(ThermResp(1:id_basetherm)); %mol O2 m-2 respired and ventilated during winter

%% Calculate uncertainty in dissolved oxygen

%Compare Year 2 Apex Mooring measurements from different depths
[C,Isb,Irid] = intersect(Yr2_sb.time_mat(Yr2_sb.O2ind),Yr2_rid.time_mat(Yr2_rid.O2ind));
O2diff = Yr2_sb.O2cleaned(Isb) - Yr2_rid.O2cleaned(Irid);

%Uncertainty in offset corrections between gliders and other assets (stdev)
Yr1_surfoffset_delta;
Yr2_surfoffset_delta;

%Propagated uncertainty in fixed depth gain corrections
    O2surf_mean_Yr1 = [nanmean(Yr1_flmA.O2cleaned) nanmean(Yr1_flmB.O2cleaned) nanmean(Yr1_sb.O2cleaned)]; %use to covert to percentages
    O2surf_mean_Yr2 = [nanmean(Yr2_flmA.O2cleaned) nanmean(Yr2_rid.O2cleaned) nanmean(Yr2_sb.O2cleaned)]; %use to covert to percentages
Yr1_percenterr_total = sqrt((2*Yr1_surfoffset_delta./O2surf_mean_Yr1.^2) + (gliders_percentaccuracy(1)/100).^2);
Yr2_percenterr_total = sqrt((2*Yr2_surfoffset_delta./O2surf_mean_Yr2.^2) + (gliders_percentaccuracy(3)/100).^2);



