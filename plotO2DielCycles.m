%% Script to look at O2 diel cycles in OOI Irminger Sea data
% Written for possible collaboration with Nathan Briggs based on
% conversation at Irminger Sea workshop
% 5 December 2017, H. Palevsky

% Run this after OceanographyMagProcessing.m

figure(11); clf
    subplot(313)
plot(Yr2_sb.time_mat(Yr2_sb.O2ind(2:end)), (Yr2_sb.O2cleaned(2:end)), '.','color',C_sb,'linewidth',L); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.O2ind(2:Yr2_rid.timebreak)),(Yr2_rid.O2cleaned(2:Yr2_rid.timebreak)),'.','color',C_rid,'linewidth',L); hold on;
xlim([datenum(2015,8,15) datenum(2015,9,15)]);
datetick('x','keeplimits')
ylabel('\mumol kg^{-1}')
title('OOI Irminger Sea surface nighttime O_2 data from late summer 2015')

    subplot(311)
plot(Yr1_flmA.time_mat(Yr1_flmA.O2ind(2:end)), (Yr1_flmA.O2cleaned(2:end)), '.','color',C_flmA,'linewidth',L); hold on;
%plot(Yr1_flmB.time_mat(Yr1_flmB.O2ind(2:end)), (Yr1_flmB.O2cleaned(2:end)), '.','color',C_flmB,'linewidth',L); hold on;
plot(Yr1_sb.time_mat(Yr1_sb.nightind(2:end)), (Yr1_sb.O2cleaned(2:end)), '.','color',C_sb,'linewidth',L); hold on;
xlim([datenum(2014,10,15) datenum(2014,11,1)]);
datetick('x','keeplimits')
ylabel('\mumol kg^{-1}')
title('OOI Irminger Sea surface nighttime O_2 data from early fall 2014')

    subplot(312)
plot(Yr1_flmA.time_mat(Yr1_flmA.O2ind(2:end)), (Yr1_flmA.O2cleaned(2:end)), '.','color',C_flmA,'linewidth',L); hold on;
plot(Yr1_flmB.time_mat(Yr1_flmB.O2ind(2:end)), (Yr1_flmB.O2cleaned(2:end)), '.','color',C_flmB,'linewidth',L); hold on;
xlim([datenum(2015,6,23) datenum(2015,7,16)]);
datetick('x','keeplimits')
ylabel('\mumol kg^{-1}')
title('OOI Irminger Sea surface nighttime O_2 data from end of spring bloom 2015')

%% Plot respiration rates (Figure 6 from Oceanography paper) in units of per day, and along with rough C flux curve from Nathan

figure(10); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=10;
    height=8;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])  
    M = 8;

%Calculate respiration rate as mg C m-3 d-1 (matching nathan's units)
    meandens = nanmean(wfpmerge.pdens');
plot(((O2max - O2min)./(O2mindate - O2maxdate))*nanmean(meandens/1000)/1.4*12, wfpmerge.depth_grid(depth_id), 'k.'); hold on;
set(gca,'YDir','reverse');
title('Respiration rate'); xlabel('mg C m^{-3} d^{-1}')

%Calculate respiration rate as mg C m-2 d-1 (matching what I think nathan's
%units actually are)
ThermResp = (O2max - O2min)./(O2mindate - O2maxdate).*meandens(depth_id)*(wfpmerge.depth_grid(2) - wfpmerge.depth_grid(1))/1000/1.4; %respiration per depth interval (mmol C/m2/d)
    intdepths = [200:5:1000]; %choose integration depth for base of seasonal thermocline (top is 200 m because that's where data starts)
    id_basetherm = find(wfpmerge.depth_grid == 1000);
for i = 1:length(intdepths)
    id_toptherm = find(wfpmerge.depth_grid == intdepths(i));
    ThermResp_Int(i) = nansum(ThermResp(id_toptherm:id_basetherm)); %mmol C m-2 d-1 respired and ventilated during winter
end

figure(11); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=10;
    height=8;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])  
    M = 8;

%Calculate respiration rate as mmol C m-2 d-1
plot(ThermResp_Int, intdepths, 'k.'); hold on;
set(gca,'YDir','reverse');
title('Carbon Flux (normalized to zero at 1000 m)'); xlabel('mmol C m^{-2} d^{-1}')

%% Calculate daily production and respiration from O2 diel cycles
[sunrise,sunset] = SunriseSunset(60,-39,0);
A = [sunrise; sunset];
noonDaily = mean(A,1);
noonAnnualMean = mean(noonDaily); %for now, just use annual mean, not accounting for seasonal variation in exact local apparent noon
noonOffset = noonAnnualMean/24 - 0.5; %correct so that apparent noon is at clock noon

fitCutoff = 0.25; %minimim r^2 to keep linear fit through data in dielCycles

tic
dailyOutput_Yr1_flmA = dielCycles(Yr1_flmA.time_mat(Yr1_flmA.O2ind(2:end)) - noonOffset, (Yr1_flmA.O2cleaned(2:end)), 0.5, fitCutoff);
dailyOutput_Yr1_flmB = dielCycles(Yr1_flmB.time_mat(Yr1_flmB.O2ind(2:end)) - noonOffset, (Yr1_flmB.O2cleaned(2:end)), 0.5, fitCutoff);
dailyOutput_Yr1_sb = dielCycles(Yr1_sb.time_mat(Yr1_sb.nightind(2:10:end)) - noonOffset, (Yr1_sb.O2cleaned(2:10:end)), 0.5, fitCutoff);
dailyOutput_Yr2_flmA = dielCycles(Yr2_flmA.time_mat(Yr2_flmA.O2ind(2:end)) - noonOffset, Yr2_flmA.O2cleaned(2:end), 0.5, fitCutoff);
dailyOutput_Yr2_sb = dielCycles(Yr2_sb.time_mat(Yr2_sb.O2ind(2:10:end)) - noonOffset, Yr2_sb.O2cleaned(2:10:end), 0.5, fitCutoff);
dailyOutput_Yr2_rid = dielCycles(Yr2_rid.time_mat(Yr2_rid.O2ind(2:10:end)) - noonOffset, Yr2_rid.O2cleaned(2:10:end), 0.5, fitCutoff);
toc

%% Merge output from all sensors
combinedDailyOutput = [dailyOutput_Yr1_flmA; dailyOutput_Yr1_flmB; dailyOutput_Yr1_sb; dailyOutput_Yr2_flmA; dailyOutput_Yr2_sb; dailyOutput_Yr2_rid];
combinedDailyOutput = sortrows(combinedDailyOutput);

    filt_time = 3; %days
    pts_per_filt_time = 1;
    
    ind_P = find(isnan(combinedDailyOutput(:,5)) == 0); %only keep rows that have GPP rate
output_filt = arrayloessfilter(combinedDailyOutput(ind_P,:), filt_time, pts_per_filt_time);
output_filt_broad = arrayloessfilter(combinedDailyOutput(ind_P,:), filt_time*3, pts_per_filt_time);

%% Plot GPP and R rates determined from diel cycles
figure(12); clf
    M = 10;
    L = 1.5;
    subplot(311)
plot(dailyOutput_Yr1_flmA(:,1), dailyOutput_Yr1_flmA(:,5), '.', 'color', C_flmA, 'markersize', M); hold on;
plot(dailyOutput_Yr1_flmB(:,1), dailyOutput_Yr1_flmB(:,5), '.', 'color', C_flmB, 'markersize', M); hold on;
plot(dailyOutput_Yr1_sb(:,1), dailyOutput_Yr1_sb(:,5), '.', 'color', C_sb, 'markersize', M); hold on;
plot(dailyOutput_Yr2_flmA(:,1), dailyOutput_Yr2_flmA(:,5), '.', 'color', C_flmA, 'markersize', M); hold on;
plot(dailyOutput_Yr2_sb(:,1), dailyOutput_Yr2_sb(:,5), '.', 'color', C_sb, 'markersize', M); hold on;
plot(dailyOutput_Yr2_rid(:,1), dailyOutput_Yr2_rid(:,5), '.', 'color', C_rid, 'markersize', M); hold on;
plot(output_filt(:,1), output_filt(:,5), '.', 'color', 'k', 'markersize', M*1.5); hold on;
plot(output_filt(:,1), output_filt(:,5), '--', 'color', 'k'); hold on;
ylabel('\mumol O_2 kg^{-1} d^{-1}');
title('Gross primary production, determined from O_2 diel cycles')
ylim([-1 60])
xlim([datenum(2014,8,1) datenum(2016,8,1)])
datetick('x',3,'keeplimits')
    subplot(312)
plot(dailyOutput_Yr1_flmA(:,1), -1*dailyOutput_Yr1_flmA(:,4), '.', 'color', C_flmA, 'markersize', M); hold on;
plot(dailyOutput_Yr1_flmB(:,1), -1*dailyOutput_Yr1_flmB(:,4), '.', 'color', C_flmB, 'markersize', M); hold on;
plot(dailyOutput_Yr1_sb(:,1), -1*dailyOutput_Yr1_sb(:,4), '.', 'color', C_sb, 'markersize', M); hold on;
plot(dailyOutput_Yr2_flmA(:,1), -1*dailyOutput_Yr2_flmA(:,4), '.', 'color', C_flmA, 'markersize', M); hold on;
plot(dailyOutput_Yr2_sb(:,1), -1*dailyOutput_Yr2_sb(:,4), '.', 'color', C_sb, 'markersize', M); hold on;
plot(dailyOutput_Yr2_rid(:,1), -1*dailyOutput_Yr2_rid(:,4), '.', 'color', C_rid, 'markersize', M); hold on;
plot(output_filt(:,1), -1*output_filt(:,4), '.', 'color', 'k', 'markersize', M*1.5); hold on;
plot(output_filt(:,1), -1*output_filt(:,4), '--', 'color', 'k'); hold on;
ylabel('\mumol O_2 kg^{-1} d^{-1}');
title('Respiration rate, determined from O_2 diel cycles')
xlim([datenum(2014,8,1) datenum(2016,8,1)])
ylim([-1 100])
datetick('x',3,'keeplimits')
    subplot(313)
plot([datenum(2014,8,1) datenum(2016,8,1)],[0 0],'b--'); hold on;
plot(output_filt_broad(:,1), output_filt_broad(:,5) + output_filt_broad(:,4), '.', 'color', 'k', 'markersize', M*1.5); hold on;
plot(output_filt_broad(:,1), output_filt_broad(:,5) + output_filt_broad(:,4), '--', 'color', 'k'); hold on;
ylabel('\mumol O_2 kg^{-1} d^{-1}');
title('NCP (GPP - R)')
xlim([datenum(2014,8,1) datenum(2016,8,1)])
ylim([-20 5])
datetick('x',3,'keeplimits')
