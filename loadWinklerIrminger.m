%% Load discrete oxygen data from Irminger deployment cruises
%Year 1
    %For calculating precision, note that first 12 samples are 6 pairs of
    %duplicates (test samples)
[num,txt,~]=xlsread('IrmingerYr1_KN221_CTDWaterSamplingData.xlsx');
date = txt(2:end,3); %pull out text dates
Yr1_disc.day = NaN*ones(length(txt)-1,1);
for i = 1:length(num)
    if cellfun('isempty',date(i)) == 1
        Yr1_disc.day(i) = NaN;
    else Yr1_disc.day(i) = datenum(date(i));
    end
end

Yr1_disc.cast = num(:,1);
Yr1_disc.time = num(:,3);
Yr1_disc.lat = num(:,4);
Yr1_disc.lon360 = -1*num(:,5) + 360; %degrees W - make negative to match convention
Yr1_disc.depth = num(:,8); 
Yr1_disc.oxy = num(:,10)/(0.0223916); %convert from mL/L to mmol/L; O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)

%Calculate uncertainty in Year 1 samples (take stdev of all 6 duplicates,
%and use mean stdev value)
for i = 1:6
    dups = Yr1_disc.oxy(2*i-1:2*i);
    duperr(i) = std(dups);
end
Yr1_disc.oxy_err = mean(duperr);


%Year 2
[num,txt,~]=xlsread('IrmingerYr2_AT30_CTDWaterSamplingData.xlsx');
date = txt(2:end,3); %pull out text dates
Yr2_disc.day = NaN*ones(length(txt)-1,1);
for i = 1:length(num)
    if cellfun('isempty',date(i)) == 1
        Yr2_disc.day(i) = NaN;
    else Yr2_disc.day(i) = datenum(date(i));
    end
end

Yr2_disc.cast = num(:,1);
Yr2_disc.time = num(:,3);
Yr2_disc.lat = num(:,4);
Yr2_disc.lon360 = -1*num(:,5) + 360; %degrees W - make negative to match convention
Yr2_disc.depth = num(:,8); 
Yr2_disc.oxy = num(:,10)/(0.0223916); %convert from mL/L to mmol/L; O2 = 22391.6 ml/mol (Garcia and Gordon, 1992)
Yr2_disc.nitrate = num(:,12); %uM

%Calculate uncertainty in Year 1 samples (take stdev of all 6 duplicates,
%and use mean stdev value)
for i = 1:6
    dups = Yr2_disc.oxy(2*i-1:2*i);
    duperr(i) = std(dups);
end
for i = 18:22
    dups = Yr2_disc.oxy(2*i:2*i+1)
    duperr(i-11) = std(dups)
end
Yr2_disc.oxy_err = mean(duperr);