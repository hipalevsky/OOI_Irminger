%% OOI Irminger Sea - extracting fluorometer data from moorings in Yrs 1-2
% H. Palevsky, Nov. 2017

%% Extract profiler mooring fluorometer data
%Wire-following profiler, Year 1, Fluorometer    
filename = ['deployment0001_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20140912T050204-20150812T103930.nc']; ncdisp(filename)
    %Time-aligned data from CTD
    Yr1_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
    Yr1_wfp.lon_flord = ncread(filename,'lon');
    Yr1_wfp.lat_flord = ncread(filename,'lat');
    Yr1_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr1_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr1_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr1_wfp.backscatter = ncread(filename,'flort_kn_bback_total'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr1_wfp.scat_total = ncread(filename,'scat_seawater'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr1_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'
    %Convert to matlab time
    Yr1_wfp.time_flord_mat = convertTime(Yr1_wfp.time_flord);

%Wire-following profiler, Year 2, Fluorometer    
filename = ['deployment0002_GI02HYPM-WFP02-01-FLORDL000-recovered_wfp-flord_l_wfp_instrument_recovered_20150817T030206-20160628T060527.nc']; ncdisp(filename)
    %Time-aligned data from CTD
    Yr2_wfp.time_flord = ncread(filename,'time'); %Note that this is the same as time_dosta
    Yr2_wfp.lon_flord = ncread(filename,'lon');
    Yr2_wfp.lat_flord = ncread(filename,'lat');
    Yr2_wfp.temperature_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_temperature'); %standard_name = 'sea_water_temperature' units = 'deg_C'
    Yr2_wfp.pracsal_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_sci_water_pracsal'); %standard_name = 'sea_water_practical_salinity'
    Yr2_wfp.pressure_flord = ncread(filename,'ctdpf_ckl_wfp_instrument_recovered-ctdpf_ckl_seawater_pressure'); %standard_name = 'sea_water_pressure' units = 'dbar'
    %Fluorometer data
    Yr2_wfp.backscatter = ncread(filename,'flort_kn_bback_total'); %long_name = 'Optical Backscatter' units = 'm-1'
    Yr2_wfp.scat_total = ncread(filename,'scat_seawater'); %long_name = 'Total Scattering Coefficient of Pure Seawater' units = 'm-1'
    Yr2_wfp.chla = ncread(filename,'fluorometric_chlorophyll_a'); %long_name = 'Chlorophyll-a Concentration' units = 'ug L-1'  
    %Convert to matlab time
    Yr2_wfp.time_flord_mat = convertTime(Yr2_wfp.time_flord);
    
%% Assign profile indices to profiler mooring data
Yr1_wfp.depth_flord = sw_dpth(Yr1_wfp.pressure_flord,Yr1_wfp.lat_flord);
    [Yr1_wfp.profile_index,Yr1_wfp.updown_index] = profileIndex(Yr1_wfp.depth_flord);

Yr2_wfp.depth_flord = sw_dpth(Yr2_wfp.pressure_flord,Yr2_wfp.lat_flord);
    [Yr2_wfp.profile_index,Yr2_wfp.updown_index] = profileIndex(Yr2_wfp.depth_flord);
    
%% Fluorometer data from fixed depth sensors on the Apex Surface Mooring
%Near surface instrument frame - only has year 2 data, lasts most of year
filename = ['deployment0002_GI01SUMO-RID16-02-FLORTD000-recovered_host-flort_dj_dcl_instrument_recovered_20150815T193016.099000-20160718T234805.764000.nc']; ncdisp(filename)
    Yr2_rid.time_fl = ncread(filename,'time');
    Yr2_rid.lat_fl = ncread(filename,'lat');
    Yr2_rid.lon_fl = ncread(filename,'lon'); 
    Yr2_rid.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr2_rid.backscatter = ncread(filename,'optical_backscatter'); 
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
    Yr1_sb.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr1_sb.backscatter = ncread(filename,'optical_backscatter'); 
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
    Yr2_sb.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    Yr2_sb.backscatter = ncread(filename,'optical_backscatter'); 
    %METBK data
    Yr2_sb.SST_fl = ncread(filename,'metbk_a_dcl_instrument_recovered-sea_surface_temperature'); 
    Yr2_sb.SSS_fl = ncread(filename,'metbk_a_dcl_instrument_recovered-met_salsurf'); 
    %Convert to matlab time
    Yr2_sb.time_mat_fl = convertTime(Yr2_sb.time_fl);
    
%% Flanking mooring fluorometer data
%%% note that there is an issue with OOI's processing algorithm so CTD data
%%% is not integrated into the output for other sensors
filename = ['deployment0001_GI03FLMA-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20140913T001501-20150818T103001.nc']; ncdisp(filename)    
    Yr1_flmA.time_fl = ncread(filename,'time');
    Yr1_flmA.lat_fl = ncread(filename,'lat');
    Yr1_flmA.lon_fl = ncread(filename,'lon'); 
    Yr1_flmA.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    %Yr1_flmA.backscatter = ncread(filename,'optical_backscatter'); all nan
    Yr1_flmA.scat_total = ncread(filename,'total_volume_scattering_coefficient'); 
    %Convert to matlab time
    Yr1_flmA.time_mat_fl = convertTime(Yr1_flmA.time_fl);
    
filename = ['deployment0001_GI03FLMB-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20140916T173001-20150820T124501.nc']; ncdisp(filename)    
    Yr1_flmB.time_fl = ncread(filename,'time');
    Yr1_flmB.lat_fl = ncread(filename,'lat');
    Yr1_flmB.lon_fl = ncread(filename,'lon'); 
    Yr1_flmB.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    %Yr1_flmB.backscatter = ncread(filename,'optical_backscatter'); all nan
    Yr1_flmB.scat_total = ncread(filename,'total_volume_scattering_coefficient'); 
    %Convert to matlab time
    Yr1_flmB.time_mat_fl = convertTime(Yr1_flmB.time_fl);
    
filename = ['deployment0002_GI03FLMA-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20150819T000001-20160418T201501.nc']; ncdisp(filename)    
    Yr2_flmA.time_fl = ncread(filename,'time');
    Yr2_flmA.lat_fl = ncread(filename,'lat');
    Yr2_flmA.lon_fl = ncread(filename,'lon'); 
    Yr2_flmA.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    %Yr2_flmA.backscatter = ncread(filename,'optical_backscatter'); all nan
    Yr2_flmA.scat_total = ncread(filename,'total_volume_scattering_coefficient'); 
    %Convert to matlab time
    Yr2_flmA.time_mat_fl = convertTime(Yr2_flmA.time_fl);
    
filename = ['deployment0002_GI03FLMB-RIS01-05-FLORTD000-recovered_host-flort_dj_sio_instrument_recovered_20150821T173001-20151122T220001.nc']; ncdisp(filename)    
    Yr2_flmB.time_fl = ncread(filename,'time');
    Yr2_flmB.lat_fl = ncread(filename,'lat');
    Yr2_flmB.lon_fl = ncread(filename,'lon'); 
    Yr2_flmB.chla = ncread(filename,'fluorometric_chlorophyll_a'); 
    %Yr2_flmB.backscatter = ncread(filename,'optical_backscatter'); %Optical Backscatter (Red Wavelengths) is a measure of the amount of red light (630-740 nm wavelengths) scattered in the backward direction due to suspended matter within seawater, providing a proxy for turbidity and suspended solids. Units: m-1 --> Note that this is calculated from total volume scattering coefficient - seawater scattering coefficient
    Yr2_flmB.scat_total = ncread(filename,'total_volume_scattering_coefficient'); 
    %Convert to matlab time
    Yr2_flmB.time_mat_fl = convertTime(Yr2_flmB.time_fl);