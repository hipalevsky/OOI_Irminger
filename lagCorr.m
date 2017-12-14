function [O2_lagcorr] = lagCorr(O2_in,time_in,tau,timetol,O2min)

% Correct for time lag in profiles of oxygen data (though could
% theoretically be used for other properties)
% O2corr = O2raw + tau*(dO2/dt) -equation 1, Nicholson et al. 2008
    % tau is timescale of lag (30 s from that paper)
% Implemented by H. Palevsky, 4/7/2017

%%%% INPUTS
% O2_in - timeseries of O2 data from profiles
% time_in - timestamp corresponding to O2_in (matlab timestamp, where unit
% = 1 day)
% tau - lag timescale in seconds
% timetol - maximum time between measurements in seconds
% O2min - minimum value of O2_in to use as good data
    
%%%% OUTPUT
% O2_lagcorr - same length as O2_in, with NaNs in spots without viable
% lag-corrected data

secinday = 60*60*24;
    
O2_lagcorr = NaN*(ones(size(O2_in)));
    dt = diff(time_in)*secinday;
    dO2 = diff(O2_in);
    indkeep = find(dt < timetol & O2_in(2:end) > O2min & O2_in(1:end-1) > O2min); %remove time gaps at profile breaks and anomalous low oxygens
        twindow = nearest(tau/nanmean(dt(indkeep))); %calculate number of points to use in moving mean to calculate dO2/dt
O2_lagcorr(indkeep) = O2_in(indkeep) + tau*movmean((dO2(indkeep)./dt(indkeep)),[twindow 0]); %use trailing mean dO2/dt over twindow (based on tau)

end