function [dailyOutput] = dielCycles(time, O2, localNoon, fitCutoff)

%%% INPUTS
% time - time in Matlab datenumformat
% O2 - measured O2, same length as time (umol/kg)
% localNoon - time of local apparent noon (value between 0 and 1 - i.e. clock noon = 0.5)
% fitCutoff - minimim r^2 to keep linear fit through data

%%% OUTPUT
% dailyOutput - array of:
    % day (Matlab datetime)
    % O2 at noon yesterday, extrapolated from night O2 fit (umol/kg)
    % O2 at noon today, extrapolated from night O2 fit (umol/kg)
    % Respiration rate (umol/kg/day)
    % O2 primary production = O2 expected at noon today, extrapolated from
        % last night's respiration - O2 expected at noon today, extrapolated
        % from tommorow night's respiration, where difference is due to
        % today's photosynthetic production (umol/kg/day)

year = str2num(datestr(time,10)); %year of sampling
julianday = time - datenum(year,0,0); %Julian day w/o year to create annual composite
timeofday = julianday - floor(julianday); %time of day (UTC, from 0 = midnight to 1 = midnight of subsequent day)
day = datenum(year,0,0) + floor(julianday); %year-day for each point
daylist = unique(day); %all unique year-days

idMorning = find(timeofday < localNoon);
idAfternoon = find(timeofday > localNoon);
minPointsNeeded = 2; %minimum number of points from each of evening and morning

%Initialize array to hold daily output:
    % Day, O2 at noon yesterday, O2 at noon today, respiration rate,
    % production rate
dailyOutput = NaN*ones(length(daylist),5);
dailyOutput(:,1) = daylist;

for i = 2:length(daylist) %loop starting on 2nd day and going to final day of record
    today = find(day == daylist(i));
    yesterday = find(day == (daylist(i)-1));
    thisMorning = intersect(today,idMorning);
    yesterdayAfternoon = intersect(yesterday,idAfternoon);
    if length(thisMorning) >= minPointsNeeded & length(yesterdayAfternoon) >= minPointsNeeded %if there are enough data points to work with
%%%%%% Next step is to fit a line through the O2 data for each night
        xx = [time(yesterdayAfternoon); time(thisMorning)];
        yy = [O2(yesterdayAfternoon); O2(thisMorning)];
        p = polyfit(xx, yy, 1);
        [rho,df,rho_sig95] = correlate(xx, yy);
%%%%%% For development, plot a figure showing the O2 data each night
%         figure(100 + i); clf
%         plot(xx,yy,'k.'); hold on;
%         plot(xx, p(1)*xx + p(2),'r--'); hold on;
%         datetick('x')
%         title([datestr(daylist(i)) ', r^2 = ' num2str(rho.^2)])
        if rho.^2 > fitCutoff & p(1) < 0 %only keep values if there is a reasonable fit and if O2 is decreasing over nighttime
            %Extrapolate the linear fit to noon yesterday and noon today
            dailyOutput(i,2) = p(1)*(localNoon + daylist(i)-1) + p(2);
            dailyOutput(i,3) = p(1)*(localNoon + daylist(i)) + p(2);
            %Keep the slope (which is the respiration rate, in umol/kg/day)
            dailyOutput(i,4) = p(1);
%         else
%             disp(['Data do not meet linear fit (r^2 < ' num2str(fitCutoff) ') for ' datestr(daylist(i))])
        end
        if dailyOutput(i,1) - dailyOutput(i-1,1) == 1 %if consecutive days in sequence
            dailyOutput(i-1,5) = dailyOutput(i,2) - dailyOutput(i-1,3); %take difference between O2 predicted at noon yesterday and yesterday's fit to noon - difference is production
        end
%     else
%         disp(['Insufficient data for ' datestr(daylist(i))])
    end
end

end