%% Extract Year 1 data

%Wire-following profiler, Year 1, CTD
filename = ['deployment0001_GI02HYPM-WFP02-04-CTDPFL000-recovered_wfp-ctdpf_ckl_wfp_instrument_recovered_20140912T050002-20150812T104141.998194.nc']; ncdisp(filename)
%internal_timestamp = ncread(filename,'internal_timestamp'); %long_name   = 'Internal Timestamp, UTC' ' units = 'seconds since 1900-01-01' %note that there are other timestamp options
    Yr1_wfp.time_ctd = ncread(filename,'time'); %same as internal_timestamp
    Yr1_wfp.lon_ctd = ncread(filename,'lon');
    Yr1_wfp.lat_ctd = ncread(filename,'lat');
    Yr1_wfp.pracsal_ctd = ncread(filename,'ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_ctd = ncread(filename,'ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    Yr1_wfp.temperature_ctd = ncread(filename,'ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'

%Wire-following profiler, Year 1, DOSTA    
filename = ['deployment0001_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20140912T050204-20150812T103930.nc']; ncdisp(filename)
    Yr1_wfp.time_dosta = ncread(filename,'time');
    Yr1_wfp.lon_dosta = ncread(filename,'lon');
    Yr1_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr1_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr1_wfp.oxygen = ncread(filename,'dosta_ln_wfp_abs_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr1_wfp.optode_temperature = ncread(filename,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'
        %Note that these points fall close to, but not exactly on, a 1:1
        %line with the CTD temperature. Points with zero value of optode
        %temperature may be indicator of bad oxygen data points.
    %Convert to matlab time
    Yr1_wfp.time_dosta_mat = convertTime(Yr1_wfp.time_dosta);

%Wire-following profiler, Year 1, Fluorometer    
filename = ['deployment0001_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20140912T050204-20150812T103930.nc']; ncdisp(filename)
    Yr1_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
    Yr1_wfp.lon_flord = ncread(filename,'lon');
    Yr1_wfp.lat_flord = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr1_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr1_wfp.backscatter = ncread(filename,'flort_kn_bback_total'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr1_wfp.scat_total = ncread(filename,'scat_seawater'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr1_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
        %Note that there are other scattering based results that I don't
        %really understand - look back at these if I actually want to use
        %the backscatter or scat_total results
    %Convert to matlab time
    Yr1_wfp.time_flord_mat = convertTime(Yr1_wfp.time_flord);
        
%% Extract Year 2 data
%Wire-following profiler, Year 2, DOSTA    
filename = ['deployment0002_GI02HYPM-WFP02-03-DOSTAL000-recovered_wfp-dosta_ln_wfp_instrument_recovered_20150817T030206-20160628T060527.nc']; ncdisp(filename)
    Yr2_wfp.time_dosta = ncread(filename,'time');
    Yr2_wfp.lon_dosta = ncread(filename,'lon');
    Yr2_wfp.lat_dosta = ncread(filename,'lat');
    %CTD data - Note that these appear to be directly taken from corresponding points in
    %CTD file, and so could be used without having to pull from the CTD data
    Yr2_wfp.temperature_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr2_wfp.pracsal_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr2_wfp.pressure_dosta = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Optode data
    Yr2_wfp.oxygen = ncread(filename,'dosta_ln_wfp_abs_oxygen'); %standard_name = 'moles_of_oxygen_per_unit_mass_in_sea_water' units = 'umol kg-1'
    Yr2_wfp.optode_temperature = ncread(filename,'optode_temperature'); %long_name = 'Optode Temperature' units = 'deg_C'   
    %Convert to matlab time
    Yr2_wfp.time_dosta_mat = convertTime(Yr2_wfp.time_dosta);

%Wire-following profiler, Year 2, Fluorometer    
filename = ['deployment0002_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150817T030206-20160628T060527.nc']; ncdisp(filename)
    %All of these are the same as the DOSTA points
%     Yr2_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
%     Yr2_wfp.lon_flord = ncread(filename,'lon');
%     Yr2_wfp.lat_flord = ncread(filename,'lat');
%     Yr2_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
%     Yr2_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
%     Yr2_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr2_wfp.backscatter = ncread(filename,'flort_kn_bback_total'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr2_wfp.scat_total = ncread(filename,'scat_seawater'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr2_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'  

 %% Assign profile indices prior to gridding
Yr1_wfp.depth_dosta = sw_dpth(Yr1_wfp.pressure_dosta,Yr1_wfp.lat_dosta);
    [Yr1_wfp.profile_index,Yr1_wfp.updown_index] = profileIndex(Yr1_wfp.depth_dosta);

Yr2_wfp.depth_dosta = sw_dpth(Yr2_wfp.pressure_dosta,Yr2_wfp.lat_dosta);
    [Yr2_wfp.profile_index,Yr2_wfp.updown_index] = profileIndex(Yr2_wfp.depth_dosta);

    
%% Apply a lag correction prior to gridding (use if not pairing up and down profiles)
% O2corr = O2raw + tau*(dO2/dt) -equation 1, Nicholson et al. 2008
    % tau is timescale of lag (30 s from that paper)
timetol = 10*60; %in seconds
O2min = 50; %throw out any O2 data below this cutoff
tau = 30; %seconds of lag

[Yr1_wfp.oxygen_lagcorr] = lagCorr(Yr1_wfp.oxygen,Yr1_wfp.time_dosta_mat,tau,timetol,O2min);
[Yr2_wfp.oxygen_lagcorr] = lagCorr(Yr2_wfp.oxygen,Yr2_wfp.time_dosta_mat,tau,timetol,O2min);

%% Calculate density in raw profiles to enable gridding on density surfaces
Yr1_wfp.pdens = sw_pden(Yr1_wfp.pracsal_dosta, Yr1_wfp.temperature_dosta, Yr1_wfp.pressure_dosta, 0);
Yr2_wfp.pdens = sw_pden(Yr2_wfp.pracsal_dosta, Yr2_wfp.temperature_dosta, Yr2_wfp.pressure_dosta, 0);

%% Grid data to consistent depth intervals for each profile
depth_grid = [150:5:2600];
dens_grid = [1027.5:0.01:1028];
therm_grid = [1.1:0.05:5.6];
secinday = 60*60*24;

%All profiles for year 1
scivars = [Yr1_wfp.temperature_dosta, Yr1_wfp.pracsal_dosta, Yr1_wfp.oxygen, Yr1_wfp.optode_temperature...
        Yr1_wfp.backscatter, Yr1_wfp.scat_total, Yr1_wfp.chla];
[Yr1_wfpgrid] = glider_grid(Yr1_wfp.time_dosta,Yr1_wfp.lat_dosta,Yr1_wfp.lon_dosta,Yr1_wfp.depth_dosta,Yr1_wfp.profile_index,Yr1_wfp.updown_index',scivars,depth_grid);
    Yr1_wfpgrid.depth_grid = depth_grid;
Yr1_wfpgrid.time_start = convertTime(Yr1_wfpgrid.time_start);
Yr1_wfpgrid.duration = Yr1_wfpgrid.duration/secinday;
Yr1_wfpgrid.updown = Yr1_wfpgrid.profile_direction;
%Grid on density surfaces (note that currently not correcting for drifting
%conductivity sensor)
[Yr1_wfpgrid_dens] = glider_grid_dens(Yr1_wfp.time_dosta,Yr1_wfp.lat_dosta,Yr1_wfp.lon_dosta,...
    Yr1_wfp.pdens,Yr1_wfp.profile_index,Yr1_wfp.updown_index',[scivars, Yr1_wfp.depth_dosta],dens_grid);
    Yr1_wfpgrid_dens.dens_grid = dens_grid;
Yr1_wfpgrid_dens.time_start = convertTime(Yr1_wfpgrid_dens.time_start);
Yr1_wfpgrid_dens.duration = Yr1_wfpgrid_dens.duration/secinday;
Yr1_wfpgrid_dens.updown = Yr1_wfpgrid_dens.profile_direction;
% Grid on isotherms
    [~,deepind] = unique(Yr1_wfp.temperature_dosta); %remove duplicates so can use interp1 function
[Yr1_wfpgrid_therm] = glider_grid_dens(Yr1_wfp.time_dosta(deepind),Yr1_wfp.lat_dosta(deepind),Yr1_wfp.lon_dosta(deepind),...
    Yr1_wfp.temperature_dosta(deepind),Yr1_wfp.profile_index(deepind),Yr1_wfp.updown_index(deepind)',[scivars(deepind,:), Yr1_wfp.depth_dosta(deepind)],therm_grid);
    Yr1_wfpgrid_therm.therm_grid = therm_grid;
Yr1_wfpgrid_therm.time_start = convertTime(Yr1_wfpgrid_therm.time_start);
Yr1_wfpgrid_therm.duration = Yr1_wfpgrid_therm.duration/secinday;
Yr1_wfpgrid_therm.updown = Yr1_wfpgrid_therm.profile_direction;

%All profiles for year 2
scivars = [Yr2_wfp.temperature_dosta, Yr2_wfp.pracsal_dosta, Yr2_wfp.oxygen, Yr2_wfp.optode_temperature...
        Yr2_wfp.backscatter, Yr2_wfp.scat_total, Yr2_wfp.chla];
[Yr2_wfpgrid] = glider_grid(Yr2_wfp.time_dosta,Yr2_wfp.lat_dosta,Yr2_wfp.lon_dosta,Yr2_wfp.depth_dosta,Yr2_wfp.profile_index,Yr2_wfp.updown_index',scivars,depth_grid);
    Yr2_wfpgrid.depth_grid = depth_grid;
Yr2_wfpgrid.time_start = convertTime(Yr2_wfpgrid.time_start);
Yr2_wfpgrid.duration = Yr2_wfpgrid.duration/secinday;
Yr2_wfpgrid.updown = Yr2_wfpgrid.profile_direction;  
%Grid on density surfaces (note that currently not correcting for drifting
%conductivity sensor)
[Yr2_wfpgrid_dens] = glider_grid_dens(Yr2_wfp.time_dosta,Yr2_wfp.lat_dosta,Yr2_wfp.lon_dosta,...
    Yr2_wfp.pdens,Yr2_wfp.profile_index,Yr2_wfp.updown_index',[scivars, Yr2_wfp.depth_dosta],dens_grid);
    Yr2_wfpgrid_dens.dens_grid = dens_grid;
Yr2_wfpgrid_dens.time_start = convertTime(Yr2_wfpgrid_dens.time_start);
Yr2_wfpgrid_dens.duration = Yr2_wfpgrid_dens.duration/secinday;
Yr2_wfpgrid_dens.updown = Yr2_wfpgrid_dens.profile_direction;
% Grid on isotherms
    [~,deepind] = unique(Yr2_wfp.temperature_dosta); %remove duplicates so can use interp1 function
[Yr2_wfpgrid_therm] = glider_grid_dens(Yr2_wfp.time_dosta(deepind),Yr2_wfp.lat_dosta(deepind),Yr2_wfp.lon_dosta(deepind),...
    Yr2_wfp.temperature_dosta(deepind),Yr2_wfp.profile_index(deepind),Yr2_wfp.updown_index(deepind)',[scivars(deepind,:), Yr2_wfp.depth_dosta(deepind)],therm_grid);
    Yr2_wfpgrid_therm.therm_grid = therm_grid;
Yr2_wfpgrid_therm.time_start = convertTime(Yr2_wfpgrid_therm.time_start);
Yr2_wfpgrid_therm.duration = Yr2_wfpgrid_therm.duration/secinday;
Yr2_wfpgrid_therm.updown = Yr2_wfpgrid_therm.profile_direction;

%% Take mean of paired up and down profiles
    tol = 1; %only combine profiles where time_start is < 1 day apart
[Yr1_wfpgrid.scivars_pair,Yr1_wfpgrid.ind_pair] = profilePairMean(Yr1_wfpgrid,tol);
[Yr2_wfpgrid.scivars_pair,Yr2_wfpgrid.ind_pair] = profilePairMean(Yr2_wfpgrid,tol);

[Yr1_wfpgrid_dens.scivars_pair,Yr1_wfpgrid_dens.ind_pair] = profilePairMean(Yr1_wfpgrid_dens,tol);
[Yr2_wfpgrid_dens.scivars_pair,Yr2_wfpgrid_dens.ind_pair] = profilePairMean(Yr2_wfpgrid_dens,tol);
[Yr1_wfpgrid_therm.scivars_pair,Yr1_wfpgrid_therm.ind_pair] = profilePairMean(Yr1_wfpgrid_therm,tol);
[Yr2_wfpgrid_therm.scivars_pair,Yr2_wfpgrid_therm.ind_pair] = profilePairMean(Yr2_wfpgrid_therm,tol);
    
%% Unpack scivars in gridded form
%When using scivars, gets all profiles (both up and down)
%When using scivars_pair, takes mean of paired up and down profiles

% Year 1
Yr1_wfpgrid.T = squeeze(Yr1_wfpgrid.scivars_pair(:,1,:));
Yr1_wfpgrid.S = squeeze(Yr1_wfpgrid.scivars_pair(:,2,:));
Yr1_wfpgrid.O2conc = squeeze(Yr1_wfpgrid.scivars_pair(:,3,:));
Yr1_wfpgrid.optode_temperature = squeeze(Yr1_wfpgrid.scivars_pair(:,4,:));
Yr1_wfpgrid.backscatter = squeeze(Yr1_wfpgrid.scivars_pair(:,5,:));
Yr1_wfpgrid.scat_total = squeeze(Yr1_wfpgrid.scivars_pair(:,6,:));
Yr1_wfpgrid.chla = squeeze(Yr1_wfpgrid.scivars_pair(:,7,:));
Yr1_wfpgrid.pdens = sw_dens0(Yr1_wfpgrid.S,Yr1_wfpgrid.T); 
Yr1_wfpgrid.press = sw_pres(repmat(Yr1_wfpgrid.depth_grid,length(Yr1_wfpgrid.profile_ind),1)',...
        repmat(Yr1_wfpgrid.lat,1,length(Yr1_wfpgrid.depth_grid))');

Yr1_wfpgrid_dens.T = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,1,:));
Yr1_wfpgrid_dens.S = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,2,:));
Yr1_wfpgrid_dens.O2conc = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,3,:));
Yr1_wfpgrid_dens.optode_temperature = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,4,:));
Yr1_wfpgrid_dens.backscatter = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,5,:));
Yr1_wfpgrid_dens.scat_total = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,6,:));
Yr1_wfpgrid_dens.chla = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,7,:));  
Yr1_wfpgrid_dens.depth = squeeze(Yr1_wfpgrid_dens.scivars_pair(:,8,:));  

Yr1_wfpgrid_therm.T = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,1,:));
Yr1_wfpgrid_therm.S = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,2,:));
Yr1_wfpgrid_therm.O2conc = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,3,:));
Yr1_wfpgrid_therm.optode_temperature = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,4,:));
Yr1_wfpgrid_therm.backscatter = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,5,:));
Yr1_wfpgrid_therm.scat_total = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,6,:));
Yr1_wfpgrid_therm.chla = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,7,:));  
Yr1_wfpgrid_therm.depth = squeeze(Yr1_wfpgrid_therm.scivars_pair(:,8,:)); 
    
% Year 2
Yr2_wfpgrid.T = squeeze(Yr2_wfpgrid.scivars_pair(:,1,:));
Yr2_wfpgrid.S = squeeze(Yr2_wfpgrid.scivars_pair(:,2,:));
Yr2_wfpgrid.O2conc = squeeze(Yr2_wfpgrid.scivars_pair(:,3,:));
Yr2_wfpgrid.optode_temperature = squeeze(Yr2_wfpgrid.scivars_pair(:,4,:));
Yr2_wfpgrid.backscatter = squeeze(Yr2_wfpgrid.scivars_pair(:,5,:));
Yr2_wfpgrid.scat_total = squeeze(Yr2_wfpgrid.scivars_pair(:,6,:));
Yr2_wfpgrid.chla = squeeze(Yr2_wfpgrid.scivars_pair(:,7,:));
Yr2_wfpgrid.pdens = sw_dens0(Yr2_wfpgrid.S,Yr2_wfpgrid.T); 
Yr2_wfpgrid.press = sw_pres(repmat(Yr2_wfpgrid.depth_grid,length(Yr2_wfpgrid.profile_ind),1)',...
        repmat(Yr2_wfpgrid.lat,1,length(Yr2_wfpgrid.depth_grid))');      

Yr2_wfpgrid_dens.T = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,1,:));
Yr2_wfpgrid_dens.S = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,2,:));
Yr2_wfpgrid_dens.O2conc = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,3,:));
Yr2_wfpgrid_dens.optode_temperature = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,4,:));
Yr2_wfpgrid_dens.backscatter = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,5,:));
Yr2_wfpgrid_dens.scat_total = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,6,:));
Yr2_wfpgrid_dens.chla = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,7,:));  
Yr2_wfpgrid_dens.depth = squeeze(Yr2_wfpgrid_dens.scivars_pair(:,8,:));   

Yr2_wfpgrid_therm.T = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,1,:));
Yr2_wfpgrid_therm.S = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,2,:));
Yr2_wfpgrid_therm.O2conc = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,3,:));
Yr2_wfpgrid_therm.optode_temperature = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,4,:));
Yr2_wfpgrid_therm.backscatter = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,5,:));
Yr2_wfpgrid_therm.scat_total = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,6,:));
Yr2_wfpgrid_therm.chla = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,7,:));  
Yr2_wfpgrid_therm.depth = squeeze(Yr2_wfpgrid_therm.scivars_pair(:,8,:)); 
    
%Calculate O2 saturation    
    O2equil = O2sol(Yr1_wfpgrid.S,Yr1_wfpgrid.T);
Yr1_wfpgrid.O2sat = (Yr1_wfpgrid.O2conc./O2equil - 1)*100;   
    O2equil = O2sol(Yr2_wfpgrid.S,Yr2_wfpgrid.T);
Yr2_wfpgrid.O2sat = (Yr2_wfpgrid.O2conc./O2equil - 1)*100;
    O2equil = O2sol(Yr1_wfpgrid_dens.S,Yr1_wfpgrid_dens.T);
Yr1_wfpgrid_dens.O2sat = (Yr1_wfpgrid_dens.O2conc./O2equil - 1)*100;
    O2equil = O2sol(Yr2_wfpgrid_dens.S,Yr2_wfpgrid_dens.T);
Yr2_wfpgrid_dens.O2sat = (Yr2_wfpgrid_dens.O2conc./O2equil - 1)*100;
    O2equil = O2sol(Yr1_wfpgrid_therm.S,Yr1_wfpgrid_therm.T);
Yr1_wfpgrid_therm.O2sat = (Yr1_wfpgrid_therm.O2conc./O2equil - 1)*100;
    O2equil = O2sol(Yr2_wfpgrid_therm.S,Yr2_wfpgrid_therm.T);
Yr2_wfpgrid_therm.O2sat = (Yr2_wfpgrid_therm.O2conc./O2equil - 1)*100;

%% If using separate up and down profiles, show comparison (don't do this if using profilePairMean above)
%plotUpDownProfileComparisonWFP

%% Visualize gridded data
for i = 1:2
    if i == 1
        plotting = Yr1_wfpgrid;
    elseif i == 2
        plotting = Yr2_wfpgrid;
    end

figure(i); clf
set(gcf,'color','w')
x0=1;
y0=1;
width=28;
height=20;
set(gcf,'units','centimeters','position',[x0,y0,width,height]) 
    subplot(4,2,1)
imagesc(plotting.T); colorbar; title('Temperature');
    subplot(4,2,2)
imagesc(plotting.S); colorbar; title('Salinity');
    subplot(4,2,3)
imagesc(plotting.pdens); colorbar; title('Density');
    subplot(4,2,4)
imagesc(plotting.O2conc); colorbar; caxis([240 300]); title('O_2 concentration');
    subplot(4,2,5)
imagesc(plotting.O2sat); colorbar; caxis([-25 0]); title('O_2 sat');
    subplot(4,2,6)
imagesc(plotting.backscatter); colorbar; caxis([0 0.002]); title('Backscatter');
    subplot(4,2,7)
imagesc(plotting.scat_total); colorbar; title('Scat Total');
    subplot(4,2,8)
imagesc(plotting.chla); colorbar; caxis([0 0.3]); title('Chlorophyll a');

end

%% Calibrate profiler data using discrete Winkler measurements

wfp_Irminger_winklercalibration
gain_wfp = gain;
gain_wfpstd = gainstd(1:2);
gain_wfpnum = gainnum;

%% Plot corrected data

for i = 1:2
    if i == 1
        plotdata = Yr1_wfpgrid;
    elseif i == 2
        plotdata = Yr2_wfpgrid;
    end

figure(20+i); clf
set(gcf,'color','w')
    x0=1;
    y0=1;
    width=28;
    height=20;
    set(gcf,'units','centimeters','position',[x0,y0,width,height])   
    F = 12; %title font size
    C = cmocean('Thermal'); %C = colormap(flipud(lbmap(101,'RedBlue')));
    mindepth = 0; maxdepth = 2600;
    cints = 60; %number of contour intervals
    
    subplot(411)
        cmin = 1; cmax = 6;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
    [X,Y] = meshgrid(plotdata.time_start(plotdata.ind_pair),plotdata.depth_grid);
    contourf(X,Y,squeeze(plotdata.T),cvec,'linecolor','none'); hold on;
        axis([min(plotdata.time_start) max(plotdata.time_start) mindepth maxdepth]);
    colormap(C); hcb = colorbar; set(gca,'YDir','reverse'); datetick('x',2,'keeplimits'); ylabel('Depth (m)');
    title(['Wire-following profiler Year ' num2str(i)],'Fontsize',F); title(hcb,'Temperature (^oC)')
    
    subplot(412)
        cmin = 250; cmax = 315;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
    contourf(X,Y,plotdata.O2conc*gain(i),cvec,'linecolor','none'); hold on;
        axis([min(plotdata.time_start) max(plotdata.time_start) mindepth maxdepth]); caxis([cmin cmax]);
    colormap(C); hcb = colorbar; set(gca,'YDir','reverse'); datetick('x',2,'keeplimits'); ylabel('Depth (m)');
    title(['Wire-following profiler Year ' num2str(i)],'Fontsize',F); title(hcb,'Corrected O_2 conc (from Winklers)')
    
    subplot(413)
        cmin = 20; cmax = 40;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
    contourf(X,Y,plotdata.O2conc*gain(i) - plotdata.O2conc,cvec,'linecolor','none'); hold on;
        axis([min(plotdata.time_start) max(plotdata.time_start) mindepth maxdepth]); caxis([cmin cmax]);
    colormap(C); hcb = colorbar; set(gca,'YDir','reverse'); datetick('x',2,'keeplimits'); ylabel('Depth (m)');
    title(['Wire-following profiler Year ' num2str(i)],'Fontsize',F); title(hcb,'Corrected - Uncorrected O_2 conc')
    
    subplot(414)
            cmin = -18; cmax = 0;
    cvec = [cmin:(cmax-cmin)/cints:cmax];
        O2equil = O2sol(plotdata.S,plotdata.T);
        O2sat = ((plotdata.O2conc*gain(i))./O2equil - 1)*100;
    contourf(X,Y,O2sat,cvec,'linecolor','none'); hold on;
        axis([min(plotdata.time_start) max(plotdata.time_start) mindepth maxdepth]);
    colormap(C); hcb = colorbar; set(gca,'YDir','reverse'); datetick('x',2,'keeplimits'); ylabel('Depth (m)'); caxis([cmin cmax]);
    title(['Wire-following profiler Year ' num2str(i)],'Fontsize',F); title(hcb,'O_2 % saturation')

end

%% Non-wire following mooring data extraction

%% Pulling out CTD data from flanking moorings
filename = ['deployment0001_GI03FLMA-RIM01-02-CTDMOG040-recovered_inst-ctdmo_ghqr_instrument_recovered_20140912T201501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.time_ctd30 = ncread(filename,'time');
    Yr1_flmA.lat_ctd30 = ncread(filename,'lat');
    Yr1_flmA.lon_ctd30 = ncread(filename,'lon'); 
    Yr1_flmA.T_ctd30 = ncread(filename,'ctdmo_seawater_temperature'); 
    Yr1_flmA.S_ctd30 = ncread(filename,'practical_salinity'); 
    Yr1_flmA.P_ctd30 = ncread(filename,'ctdmo_seawater_pressure'); %median pressure is 15.3
    %Convert to matlab time
    Yr1_flmA.time_mat_ctd30 = convertTime(Yr1_flmA.time_ctd30);

filename = ['deployment0002_GI03FLMA-RIM01-02-CTDMOG040-recovered_inst-ctdmo_ghqr_instrument_recovered_20150819T000001-20160716T234501.nc']; ncdisp(filename)    
    Yr2_flmA.time_ctd30 = ncread(filename,'time');
    Yr2_flmA.lat_ctd30 = ncread(filename,'lat');
    Yr2_flmA.lon_ctd30 = ncread(filename,'lon'); 
    Yr2_flmA.T_ctd30 = ncread(filename,'ctdmo_seawater_temperature'); 
    Yr2_flmA.S_ctd30 = ncread(filename,'practical_salinity'); 
    Yr2_flmA.P_ctd30 = ncread(filename,'ctdmo_seawater_pressure'); %median pressure is 37
    %Convert to matlab time
    Yr2_flmA.time_mat_ctd30 = convertTime(Yr2_flmA.time_ctd30);

filename = ['deployment0002_GI03FLMB-RIM01-02-CTDMOG060-recovered_inst-ctdmo_ghqr_instrument_recovered_20150821T171501-20160717T234501.nc']; ncdisp(filename)    
    Yr2_flmB.time_ctd30 = ncread(filename,'time');
    Yr2_flmB.lat_ctd30 = ncread(filename,'lat');
    Yr2_flmB.lon_ctd30 = ncread(filename,'lon'); 
    Yr2_flmB.T_ctd30 = ncread(filename,'ctdmo_seawater_temperature'); 
    Yr2_flmB.S_ctd30 = ncread(filename,'practical_salinity'); 
    Yr2_flmB.P_ctd30 = ncread(filename,'ctdmo_seawater_pressure'); %median pressure is 33
    %Convert to matlab time
    Yr2_flmB.time_mat_ctd30 = convertTime(Yr2_flmB.time_ctd30);

%% O2 data from flanking moorings - will have to line up manually with CTD    
%Flanking Mooring A - Optode at ~30 m - No CTD data in file so have to line
%up myself (or wait for update from OOI)
filename = ['deployment0001_GI03FLMA-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20140912T201501-20150818T103001.nc']; ncdisp(filename)
    Yr1_flmA.oxygenraw = ncread(filename,'estimated_oxygen_concentration'); %note that 'dosta_abcdjm_sio_abs_oxygen' is T and S corrected, but all NaN b/c issue w/ CTD lineup
    Yr1_flmA.time = ncread(filename,'time');
    Yr1_flmA.lat = ncread(filename,'lat');
    Yr1_flmA.lon = ncread(filename,'lon');
    %Convert to matlab time
    Yr1_flmA.time_mat = convertTime(Yr1_flmA.time);
    %Correct oxygen
    Yr1_flmA.oxygen = aaoptode_salpresscorr(Yr1_flmA.oxygenraw,Yr1_flmA.T_ctd30,Yr1_flmA.S_ctd30,Yr1_flmA.P_ctd30,S0);
    
filename = ['deployment0002_GI03FLMA-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20150819T000001-20160418T201501.nc']; ncdisp(filename)
    Yr2_flmA.oxygenraw = ncread(filename,'estimated_oxygen_concentration'); %note that 'dosta_abcdjm_sio_abs_oxygen' is T and S corrected, but all NaN b/c issue w/ CTD lineup
    Yr2_flmA.time = ncread(filename,'time'); 
    Yr2_flmA.lat = ncread(filename,'lat');
    Yr2_flmA.lon = ncread(filename,'lon');
    %Convert to matlab time
    Yr2_flmA.time_mat = convertTime(Yr2_flmA.time);
    %Correct oxygen
    Yr2_flmA.oxygen = aaoptode_salpresscorr(Yr2_flmA.oxygenraw,interp1(Yr2_flmA.time_ctd30,Yr2_flmA.T_ctd30,Yr2_flmA.time),...
        interp1(Yr2_flmA.time_ctd30,Yr2_flmA.S_ctd30,Yr2_flmA.time),interp1(Yr2_flmA.time_ctd30,Yr2_flmA.P_ctd30,Yr2_flmA.time),S0);
    
%Flanking Mooring B - Optode at ~30 m - No CTD data in file so have to line
%up myself (or wait for update from OOI)
filename = ['deployment0001_GI03FLMB-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20140916T133001-20150820T124501.nc']; ncdisp(filename)
    Yr1_flmB.oxygen = ncread(filename,'estimated_oxygen_concentration'); %note that 'dosta_abcdjm_sio_abs_oxygen' is T and S corrected, but all NaN b/c issue w/ CTD lineup
    Yr1_flmB.time = ncread(filename,'time');
    Yr1_flmB.lat = ncread(filename,'lat');
    Yr1_flmB.lon = ncread(filename,'lon');
    %Convert to matlab time
    Yr1_flmB.time_mat = convertTime(Yr1_flmB.time);
    %Note - don't have CTD data to to T,S,P correction for O2, but doesn't
    %seem to make a difference in corrected data
    
filename = ['deployment0002_GI03FLMB-RIS01-03-DOSTAD000-recovered_host-dosta_abcdjm_sio_instrument_recovered_20150821T173001-20151122T220001.nc']; ncdisp(filename)
    Yr2_flmB.oxygen = ncread(filename,'estimated_oxygen_concentration'); %note that 'dosta_abcdjm_sio_abs_oxygen' is T and S corrected, but all NaN b/c issue w/ CTD lineup
    Yr2_flmB.time = ncread(filename,'time');  
    Yr2_flmB.lat = ncread(filename,'lat');
    Yr2_flmB.lon = ncread(filename,'lon');
    %Convert to matlab time
    Yr2_flmB.time_mat = convertTime(Yr2_flmB.time);

%% Apex Surface Mooring

% Surface buoy - Year 1
filename = ['deployment0001_GI01SUMO-SBD11-04-DOSTAD000-recovered_host-dosta_abcdjm_dcl_instrument_recovered_20140910T190020.500000-20150310T050302.177000.nc']; ncdisp(filename)
    Yr1_sb.oxygen = ncread(filename,'dosta_abcdjm_dcl_metbk_abs_oxygen'); %umol kg-1, long_name = 'DO - Temp Sal Corrected (METBK)'
    Yr1_sb.oxygen_uncorr = ncread(filename,'estimated_oxygen_concentration');
    Yr1_sb.time = ncread(filename,'time');
    Yr1_sb.lat = ncread(filename,'lat');
    Yr1_sb.lon = ncread(filename,'lon');
    %CTD data - Note that these appear to be directly taken from corresponding points in CTD file, but don't last for much of year - need to check CTD file
    Yr1_sb.temperature_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-temp');
    Yr1_sb.pracsal_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-salinity');
    Yr1_sb.pressure_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-pressure');
    %METBK data
    Yr1_sb.SSS_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf');
    Yr1_sb.SST_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature');
    %Convert to matlab time
    Yr1_sb.time_mat = convertTime(Yr1_sb.time);
    
% Surface buoy - Year 2
filename = ['deployment0002_GI01SUMO-SBD11-04-DOSTAD000-recovered_host-dosta_abcdjm_dcl_instrument_recovered_20150815T193013.458000-20160130T083304.921000.nc']; ncdisp(filename)
    Yr2_sb.oxygen = ncread(filename,'dosta_abcdjm_dcl_metbk_abs_oxygen'); %umol kg-1, long_name = 'DO - Temp Sal Corrected (METBK)'
    Yr2_sb.oxygen_uncorr = ncread(filename,'estimated_oxygen_concentration');
    Yr2_sb.time = ncread(filename,'time');
    Yr2_sb.lat = ncread(filename,'lat');
    Yr2_sb.lon = ncread(filename,'lon');
    %CTD data - Note that these appear to be directly taken from corresponding points in CTD file, but don't last for much of year - need to check CTD file
    Yr2_sb.temperature_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-temp');
    Yr2_sb.pracsal_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-salinity');
    Yr2_sb.pressure_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-pressure');
    %METBK data
    Yr2_sb.SSS_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf');
    Yr2_sb.SST_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature');
    %Convert to matlab time
    Yr2_sb.time_mat = convertTime(Yr2_sb.time);
    
% Near surface instrument frame - Year 2
filename = ['deployment0002_GI01SUMO-RID16-06-DOSTAD000-recovered_host-dosta_abcdjm_dcl_instrument_recovered_20150815T193009.886000-20160718T234802.768000.nc']; ncdisp(filename)
    Yr2_rid.oxygen = ncread(filename,'dosta_abcdjm_dcl_metbk_abs_oxygen'); %umol kg-1, long_name = 'DO - Temp Sal Corrected (METBK)'
    Yr2_rid.oxygen_uncorr = ncread(filename,'estimated_oxygen_concentration');
    Yr2_rid.time = ncread(filename,'time');
    Yr2_rid.lat = ncread(filename,'lat');
    Yr2_rid.lon = ncread(filename,'lon');
    %CTD data - Note that these appear to be directly taken from corresponding points in CTD file, but don't last for much of year - need to check CTD file
    Yr2_rid.temperature_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-temp');
    Yr2_rid.pracsal_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-salinity');
    Yr2_rid.pressure_dosta = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-pressure');
    %METBK data
    Yr2_rid.SSS_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf');
    Yr2_rid.SST_dosta = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature');
    %Convert to matlab time
    Yr2_rid.time_mat = convertTime(Yr2_rid.time);
    
% Mooring riser - year 2
% None of these have usable oxygen data, but might be a data ingestion issue
% % 40 m
% filename = ['deployment0002_GI01SUMO-RII11-02-DOSTAD031-recovered_host-dosta_abcdjm_ctdbp_p_dcl_instrument_recovered_20150815T200003-20151113T170003.nc']; ncdisp(filename)
%     time = ncread(filename,'time'); %runs from 14 Aug 2015 - 12 Nov 2015
%     oxygen_uncorr = ncread(filename,'dosta_analog_tc_oxygen'); %all NaNs
% % 80 m
% filename = ['deployment0002_GI01SUMO-RII11-02-DOSTAD032-recovered_host-dosta_abcdjm_ctdbp_p_dcl_instrument_recovered_20150815T200003-20160413T230003.nc']; ncdisp(filename)
%     time = ncread(filename,'time'); %runs from 14 Aug 2015 - 12 Apr 2016
%     oxygen_uncorr = ncread(filename,'dosta_analog_tc_oxygen'); %all NaNs
% % 130 m
% filename = ['deployment0002_GI01SUMO-RII11-02-DOSTAD033-recovered_host-dosta_abcdjm_ctdbp_p_dcl_instrument_recovered_20150815T200003-20160617T000003.nc']; ncdisp(filename)
%     time = ncread(filename,'time'); %runs from 14 Aug 2015 - 16 June 2016
%     oxygen_uncorr = ncread(filename,'dosta_analog_tc_oxygen'); %all NaNs

%% Adding chlorophyll data
%Near surface instrument frame - only has year 2 data, lasts most of year
filename = ['deployment0002_GI01SUMO-RID16-02-FLORTD000-recovered_host-flort_dj_dcl_instrument_recovered_20150815T193016.099000-20160718T234805.764000.nc']; ncdisp(filename)
    Yr2_rid.time_fl = ncread(filename,'time');
    Yr2_rid.lat_fl = ncread(filename,'lat');
    Yr2_rid.lon_fl = ncread(filename,'lon'); 
    Yr2_rid.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr2_rid.bcksct = ncread(filename,'optical_backscatter'); 
    %CTD data
    Yr2_rid.temperature_fl = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-temp'); 
    Yr2_rid.pracsal_fl = ncread(filename,'ctdbp_cdef_dcl_instrument_recovered-practical_salinity'); 
    %Convert to matlab time
    Yr2_rid.time_mat_fl = convertTime(Yr2_rid.time_fl);
    
% Year 1 surface buoy
filename = ['deployment0001_GI01SUMO-SBD12-02-FLORTD000-recovered_host-flort_dj_dcl_instrument_recovered_20140910T190014.252000-20150301T093044.036000.nc']; ncdisp(filename)    
    Yr1_sb.time_fl = ncread(filename,'time');
    Yr1_sb.lat_fl = ncread(filename,'lat');
    Yr1_sb.lon_fl = ncread(filename,'lon'); 
    Yr1_sb.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr1_sb.bcksct = ncread(filename,'optical_backscatter'); 
    %METBK data
    Yr1_sb.SST_fl = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature'); 
    Yr1_sb.SSS_fl = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf'); 
    %Convert to matlab time
    Yr1_sb.time_mat_fl = convertTime(Yr1_sb.time_fl);

% Year 2 surface buoy   
filename = ['deployment0002_GI01SUMO-SBD12-02-FLORTD000-recovered_host-flort_dj_dcl_instrument_recovered_20150815T193013.074000-20160520T020257.404000.nc']; ncdisp(filename)    
    Yr2_sb.time_fl = ncread(filename,'time');
    Yr2_sb.lat_fl = ncread(filename,'lat');
    Yr2_sb.lon_fl = ncread(filename,'lon'); 
    Yr2_sb.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr2_sb.bcksct = ncread(filename,'optical_backscatter'); 
    %METBK data
    Yr2_sb.SST_fl = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature'); 
    Yr2_sb.SSS_fl = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf'); 
    %Convert to matlab time
    Yr2_sb.time_mat_fl = convertTime(Yr2_sb.time_fl);
    
%% Flanking mooring chlorophyll (note still not integrated with CTD data)
filename = ['deployment0001_GI03FLMA-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.time_fl = ncread(filename,'time');
    Yr1_flmA.lat_fl = ncread(filename,'lat');
    Yr1_flmA.lon_fl = ncread(filename,'lon'); 
    Yr1_flmA.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr1_flmA.bcksct = ncread(filename,'optical_backscatter'); 
    Yr1_flmA.bcksct_total = ncread(filename,'total_volume_scattering_coefficient'); 
    %Convert to matlab time
    Yr1_flmA.time_mat_fl = convertTime(Yr1_flmA.time_fl);
    
filename = ['deployment0001_GI03FLMB-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20140916T173001-20150820T124501.nc']; ncdisp(filename)    
    Yr1_flmB.time_fl = ncread(filename,'time');
    Yr1_flmB.lat_fl = ncread(filename,'lat');
    Yr1_flmB.lon_fl = ncread(filename,'lon'); 
    Yr1_flmB.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr1_flmB.bcksct = ncread(filename,'optical_backscatter'); 
    %Convert to matlab time
    Yr1_flmB.time_mat_fl = convertTime(Yr1_flmB.time_fl);
    
filename = ['deployment0002_GI03FLMA-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20150819T000001-20160418T201501.nc']; ncdisp(filename)    
    Yr2_flmA.time_fl = ncread(filename,'time');
    Yr2_flmA.lat_fl = ncread(filename,'lat');
    Yr2_flmA.lon_fl = ncread(filename,'lon'); 
    Yr2_flmA.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr2_flmA.bcksct = ncread(filename,'optical_backscatter'); 
    %Convert to matlab time
    Yr2_flmA.time_mat_fl = convertTime(Yr2_flmA.time_fl);
    
filename = ['deployment0002_GI03FLMB-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20150821T173001-20151122T220001.nc']; ncdisp(filename)    
    Yr2_flmB.time_fl = ncread(filename,'time');
    Yr2_flmB.lat_fl = ncread(filename,'lat');
    Yr2_flmB.lon_fl = ncread(filename,'lon'); 
    Yr2_flmB.chl = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr2_flmB.bcksct = ncread(filename,'optical_backscatter'); %Optical Backscatter (Red Wavelengths) is a measure of the amount of red light (630-740 nm wavelengths) scattered in the backward direction due to suspended matter within seawater, providing a proxy for turbidity and suspended solids. Units: m-1 --> Note that this is calculated from total volume scattering coefficient - seawater scattering coefficient
    %Convert to matlab time
    Yr2_flmB.time_mat_fl = convertTime(Yr2_flmB.time_fl);
    
%% Year 1 Flanking Mooring A - CTDs from varying depths
    
filename = ['deployment0001_GI03FLMA-RIM01-02-CTDMOG045-recovered_inst-ctdmo_ghqr_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.d170_P = ncread(filename,'ctdmo_seawater_pressure');
    Yr1_flmA.d170_T = ncread(filename,'ctdmo_seawater_temperature');
    Yr1_flmA.d170_S = ncread(filename,'practical_salinity');
    Yr1_flmA.d170_time = ncread(filename,'time');
        %Convert to matlab time
    Yr1_flmA.d170_time_mat = convertTime(Yr1_flmA.d170_time);

filename = ['deployment0001_GI03FLMA-RIM01-02-CTDMOG041-recovered_inst-ctdmo_ghqr_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.d28_P = ncread(filename,'ctdmo_seawater_pressure');
    Yr1_flmA.d28_T = ncread(filename,'ctdmo_seawater_temperature');
    Yr1_flmA.d28_S = ncread(filename,'practical_salinity');
    Yr1_flmA.d28_time = ncread(filename,'time');
            %Convert to matlab time
    Yr1_flmA.d28_time_mat = convertTime(Yr1_flmA.d28_time);
    
filename = ['deployment0001_GI03FLMA-RIM01-02-CTDMOG042-recovered_inst-ctdmo_ghqr_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.d48_P = ncread(filename,'ctdmo_seawater_pressure');
    Yr1_flmA.d48_T = ncread(filename,'ctdmo_seawater_temperature');
    Yr1_flmA.d48_S = ncread(filename,'practical_salinity');
    Yr1_flmA.d48_time = ncread(filename,'time');
            %Convert to matlab time
    Yr1_flmA.d48_time_mat = convertTime(Yr1_flmA.d48_time);
    
filename = ['deployment0001_GI03FLMA-RIM01-02-CTDMOG043-recovered_inst-ctdmo_ghqr_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.d78_P = ncread(filename,'ctdmo_seawater_pressure');
    Yr1_flmA.d78_T = ncread(filename,'ctdmo_seawater_temperature');
    Yr1_flmA.d78_S = ncread(filename,'practical_salinity');
    Yr1_flmA.d78_time = ncread(filename,'time');
            %Convert to matlab time
    Yr1_flmA.d78_time_mat = convertTime(Yr1_flmA.d78_time);
    
filename = ['deployment0001_GI03FLMA-RIM01-02-CTDMOG044-recovered_inst-ctdmo_ghqr_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.d119_P = ncread(filename,'ctdmo_seawater_pressure');
    Yr1_flmA.d119_T = ncread(filename,'ctdmo_seawater_temperature');
    Yr1_flmA.d119_S = ncread(filename,'practical_salinity');
    Yr1_flmA.d119_time = ncread(filename,'time');
            %Convert to matlab time
    Yr1_flmA.d119_time_mat = convertTime(Yr1_flmA.d119_time);

%% Separate out day and night indices
lat_irming = 60;
lon_irming = -39;
UTCoffset = 0;
tol = 1; %hrs before/after sunrise/sunset to count as daylight

[Yr1_flmA.dayind,Yr1_flmA.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr1_flmA.time_mat,tol);
[Yr1_flmB.dayind,Yr1_flmB.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr1_flmB.time_mat,tol);
[Yr2_flmA.dayind,Yr2_flmA.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_flmA.time_mat,tol);
[Yr2_flmB.dayind,Yr2_flmB.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_flmB.time_mat,tol);
[Yr1_sb.dayind,Yr1_sb.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr1_sb.time_mat,tol);
[Yr2_sb.dayind,Yr2_sb.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_sb.time_mat,tol);
[Yr2_rid.dayind,Yr2_rid.nightind] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_rid.time_mat,tol);

[Yr2_rid.dayind_fl,Yr2_rid.nightind_fl] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_rid.time_mat_fl,tol);
[Yr1_flmA.dayind_fl,Yr1_flmA.nightind_fl] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr1_flmA.time_mat_fl,tol);
[Yr1_flmB.dayind_fl,Yr1_flmB.nightind_fl] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr1_flmB.time_mat_fl,tol);
[Yr2_flmA.dayind_fl,Yr2_flmA.nightind_fl] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_flmA.time_mat_fl,tol);
[Yr2_flmB.dayind_fl,Yr2_flmB.nightind_fl] = indexDayNight(lat_irming,lon_irming,UTCoffset,Yr2_flmB.time_mat_fl,tol);

%% Plot chlorophyll fluorescence data
figure(1); clf
% plot(Yr1_sb.time_mat_fl,Yr1_sb.chl,'b.'); hold on;
% plot(Yr2_sb.time_mat_fl,Yr2_sb.chl,'b.'); hold on;
plot(Yr2_rid.time_mat_fl(Yr2_rid.nightind_fl),Yr2_rid.chl(Yr2_rid.nightind_fl),'k.'); hold on;
plot(Yr1_flmA.time_mat_fl(Yr1_flmA.nightind_fl),Yr1_flmA.chl(Yr1_flmA.nightind_fl),'r.'); hold on;
plot(Yr2_flmA.time_mat_fl(Yr2_flmA.nightind_fl),Yr2_flmA.chl(Yr2_flmA.nightind_fl),'r.'); hold on;
%plot(Yr1_flmB.time_mat_fl(Yr1_flmB.nightind_fl),Yr1_flmB.chl(Yr1_flmB.nightind_fl),'c.'); hold on;
%plot(Yr2_flmB.time_mat_fl(Yr2_flmB.nightind_fl),Yr2_flmB.chl(Yr2_flmB.nightind_fl),'c.'); hold on;
xlim([datenum(2014,9,1) datenum(2016,8,15)])
ylim([-0.1 15])
datetick('x',12,'keeplimits');
legend('Apex mooring (~12 m)','Flanking Mooring A (~30 m)','location','northwest')
ylabel('Chlorophyll (\mug/L)','Fontsize',10)
title('Surface chlorophyll (factory calibration), nighttime-only filter, Years 1-2 OOI Irminger Sea')

figure(2); clf
plot(Yr1_sb.time_mat_fl,Yr1_sb.bcksct,'b.'); hold on;
plot(Yr2_sb.time_mat_fl,Yr2_sb.bcksct,'b.'); hold on;
plot(Yr2_rid.time_mat_fl,Yr2_rid.bcksct,'k.'); hold on;
plot(Yr1_flmA.time_mat_fl,Yr1_flmA.bcksct,'r.'); hold on;
plot(Yr2_flmA.time_mat_fl,Yr2_flmA.bcksct,'r.'); hold on;
plot(Yr1_flmB.time_mat_fl,Yr1_flmB.bcksct,'c.'); hold on;
plot(Yr2_flmB.time_mat_fl,Yr2_flmB.bcksct,'c.'); hold on;
xlim([datenum(2014,9,1) datenum(2016,8,15)])
%ylim([-0.1 15])
datetick('x',12,'keeplimits');
legend('Apex mooring (~12 m)','Flanking Mooring A (~30 m)','location','northwest')
ylabel('Backscatter (m^{-1})','Fontsize',10)
title('Particle backscatter (factory calibration), Years 1-2 OOI Irminger Sea')

%Next steps:
% 1) Consider if it is worth calibrating by comparing with satellite
% product

%% Visualize what exists for non-wire following
figure(1); clf %Flanking moorings
for i = 1:2
    subplot(2,1,i)
%Year 1
plot(Yr1_flmA.time_mat(Yr1_flmA.nightind),Yr1_flmA.oxygen(Yr1_flmA.nightind),'k.'); hold on;
plot(Yr1_flmB.time_mat(Yr1_flmB.nightind),Yr1_flmB.oxygen(Yr1_flmB.nightind),'b.'); hold on;
%Year 2
plot(Yr2_flmA.time_mat(Yr2_flmA.nightind),Yr2_flmA.oxygen(Yr2_flmA.nightind),'k.'); hold on;
plot(Yr2_flmB.time_mat(Yr2_flmB.nightind),Yr2_flmB.oxygen(Yr2_flmB.nightind),'b.'); hold on;
if i == 1
    xlim([min(Yr1_flmA.time_mat) max(Yr2_flmA.time_mat)])
elseif i == 2
    plot(Yr1_flmA.time_mat(Yr1_flmA.dayind),Yr1_flmA.oxygen(Yr1_flmA.dayind),'r.'); hold on;
    plot(Yr1_flmB.time_mat(Yr1_flmB.dayind),Yr1_flmB.oxygen(Yr1_flmB.dayind),'m.'); hold on;
    xlim([datenum(2015,4,1) datenum(2015,4,10)])
end
datetick('x','keeplimits'); title('Flanking Mooring Oxygen (uncorrected)'); legend('A','B');

end
%%
figure(2); clf %Surface buoy
    subplot(211)
plot(Yr1_sb.time_mat(Yr1_sb.dayind),Yr1_sb.oxygen(Yr1_sb.dayind),'b.'); hold on; %Year 1
plot(Yr2_sb.time_mat(Yr2_sb.dayind),Yr2_sb.oxygen(Yr2_sb.dayind),'b.'); hold on; %Year 2    
plot(Yr1_sb.time_mat(Yr1_sb.nightind),Yr1_sb.oxygen(Yr1_sb.nightind),'k.'); hold on; %Year 1
plot(Yr2_sb.time_mat(Yr2_sb.nightind),Yr2_sb.oxygen(Yr2_sb.nightind),'k.'); hold on; %Year 2
datetick('x'); title('Surface Buoy Oxygen (T,S,P corrected)')
    subplot(212)
plot(Yr1_sb.time_mat,Yr1_sb.SST_dosta,'k.'); hold on; %Year 1
plot(Yr2_sb.time_mat,Yr2_sb.SST_dosta,'k.'); hold on; %Year 2
datetick('x'); title('Surface Buoy SST')

%Note that zooming in on September shows what may be real diurnal signals

%%
figure(3); clf %Near surface instrument frame
for i = 1:2
    subplot(3,1,1 + 2*(i-1))
plot(Yr2_rid.time_mat(Yr2_rid.nightind),Yr2_rid.oxygen(Yr2_rid.nightind),'k.'); hold on; %Year 2
if i == 1
    xlim([min(Yr2_rid.time_mat) max(Yr2_rid.time_mat)])
elseif i == 2
    plot(Yr2_rid.time_mat(Yr2_rid.dayind),Yr2_rid.oxygen(Yr2_rid.dayind),'b.'); hold on; %Year 2
    plot(Yr2_rid.time_mat(Yr2_rid.nightind),Yr2_rid.oxygen(Yr2_rid.nightind),'k.'); hold on; %Year 2
    xlim([datenum(2016,5,1) datenum(2016,5,10)])
end
datetick('x'); title('Near surface instrument frame oxygen (T,S,P corrected)')
end
    subplot(312)
plot(Yr2_rid.time_mat,Yr2_rid.SST_dosta,'k.'); hold on; %Year 2
datetick('x'); title('METBK SST')


%% Visualize all of what exists
Yr1_flmA.indO2 = find(isnan(Yr1_flmA.oxygen) == 0);
Yr2_flmA.indO2 = find(isnan(Yr2_flmA.oxygen) == 0);
Yr1_flmB.indO2 = find(isnan(Yr1_flmB.oxygen) == 0);
Yr2_flmB.indO2 = find(isnan(Yr2_flmB.oxygen) == 0);
Yr1_sb.indO2 = find(isnan(Yr1_sb.oxygen) == 0);
Yr2_sb.indO2 = find(isnan(Yr2_sb.oxygen) == 0);
Yr2_rid.indO2 = find(isnan(Yr2_rid.oxygen) == 0);
Yr1_wfp.indO2 = find(isnan(Yr1_wfp.oxygen) == 0);
Yr2_wfp.indO2 = find(isnan(Yr2_wfp.oxygen) == 0);

figure(10); clf
plot(Yr1_flmA.time_mat(Yr1_flmA.indO2),1*ones(size(Yr1_flmA.indO2)),'k.'); hold on;
plot(Yr2_flmA.time_mat(Yr2_flmA.indO2),1*ones(size(Yr2_flmA.indO2)),'k.'); hold on;
plot(Yr1_flmB.time_mat(Yr1_flmB.indO2),2*ones(size(Yr1_flmB.indO2)),'m.'); hold on;
plot(Yr2_flmB.time_mat(Yr2_flmB.indO2),2*ones(size(Yr2_flmB.indO2)),'m.'); hold on;
plot(Yr1_sb.time_mat(Yr1_sb.indO2),3*ones(size(Yr1_sb.indO2)),'b.'); hold on;
plot(Yr2_sb.time_mat(Yr2_sb.indO2),3*ones(size(Yr2_sb.indO2)),'b.'); hold on;
plot(Yr2_rid.time_mat(Yr2_rid.indO2),4*ones(size(Yr2_rid.indO2)),'c.'); hold on;
plot(Yr1_wfp.time_dosta_mat(Yr1_wfp.indO2),5*ones(size(Yr1_wfp.indO2)),'g.'); hold on;
plot(Yr2_wfp.time_dosta_mat(Yr2_wfp.indO2),5*ones(size(Yr2_wfp.indO2)),'g.'); hold on;
%Add in gliders
    L = 2;
plot([datenum(2014,10,1) datenum(2015,4,1)],[6 6],'r-','linewidth',L); hold on;
plot([datenum(2015,8,20) datenum(2015,11,20)],[6 6],'r-','linewidth',L); hold on;
plot([datenum(2015,8,20) datenum(2016,3,15)],[6.2 6.2],'r--','linewidth',L); hold on;
% Labels
xlim([datenum(2014,9,1) datenum(2016,8,1)]); datetick('x',20,'keeplimits');
set(gca,'ydir','reverse'); ylim([0.7 6.5]); yticks([1:6]); yticklabels({'Flanking mooring A (30 m)','Flanking mooring B (30 m)',...
    'Apex mooring surface buoy (1 m)','Apex Mooring instrument frame (12 m)','Wire-following profiler (150-2000 m)','Gliders (0-1000 m, dashed = 0-200 m)'});
