DorR = 'd'; %all are telemetered data (no gliders were recovered)

%GL001_D1
rootDir = 'D://OOI_downloads/rawdata.oceanobservatories.org/files/GI05MOAS-GL001/D00001';
    [GL001_D1_L2,GL001_D1_meta_L2,GL001_D1_L1,GL001_D1_meta_L1,GL001_D1_L0,GL001_D1_meta_L0] = ws_process(rootDir,DorR);

%GL002_D1 
rootDir = 'D://OOI_downloads/rawdata.oceanobservatories.org/files/GI05MOAS-GL002/D00001';
    [GL002_D1_L2,GL002_D1_meta_L2,GL002_D1_L1,GL002_D1_meta_L1,GL002_D1_L0,GL002_D1_meta_L0] = ws_process(rootDir,DorR);

%GL002_D2
rootDir = 'D://OOI_downloads/rawdata.oceanobservatories.org/files/GI05MOAS-GL002/D00002';
    [GL002_D2_L2,GL002_D2_meta_L2,GL002_D2_L1,GL002_D2_meta_L1,GL002_D2_L0,GL002_D2_meta_L0] = ws_process(rootDir,DorR);

%GL003_D2
rootDir = 'D://OOI_downloads/rawdata.oceanobservatories.org/files/GI05MOAS-GL003/D00002';
    [GL003_D2_L2,GL003_D2_meta_L2,GL003_D2_L1,GL003_D2_meta_L1,GL003_D2_L0,GL003_D2_meta_L0] = ws_process(rootDir,DorR);

%PG001_D1
rootDir = 'D://OOI_downloads/rawdata.oceanobservatories.org/files/GI05MOAS-PG001/D00001';
    [PG001_D1_L2,PG001_D1_meta_L2,PG001_D1_L1,PG001_D1_meta_L1,PG001_D1_L0,PG001_D1_meta_L0] = ws_process(rootDir,DorR);