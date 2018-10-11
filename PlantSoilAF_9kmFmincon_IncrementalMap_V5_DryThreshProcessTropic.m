clc
clear
load coast   
clat = lat;                                                        
clon = long ; 
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\9km Half Year')
% load('RFIflags_NaNremoved_201504_201509_9km')
% load('RFIflags_NaNremoved_201510_201603_9km')
% load('surfaceflags_NaNremoved_201504_201509_9km')
% RFIFlags = RFI_detect_TBh_201504_201509+RFI_detect_TBh_201510_201603+...
%     RFI_detect_TBv_201504_201509+RFI_detect_TBv_201510_201603;
% RFIFlags(RFIFlags>100)=NaN;
% RFIFlags(isfinite(RFIFlags))=1;
% MountainFlag = surfaceflag_mountainous_201504_201509;
% MountainFlag(MountainFlag>1)=NaN;
% MountainFlag(isfinite(MountainFlag))=1;
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy')
load('BootConfidenceAuAfSA')
% RFIMitigate = RFI_mitigate_TBh_201504_201509+RFI_mitigate_TBh_201510_201603+...
%     RFI_mitigate_TBv_201504_201509+RFI_mitigate_TBv_201510_201603;
% [ks f] = ksdensity(RFIFlags(:),'bandwidth',1);
% plot(f,ks)
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\General Data')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP Data/General Data')
% load('GPMMonSum_201504_201703_9km')
% GPMMonSum_201504_201703 = nansum(GPMMonSum_201504_201703,3);
% GPMMonSum_201504_201703(GPMMonSum_201504_201703==0)=NaN;
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
load('ForestFlag')
load('SMAPCenterCoordinates9KM')
load('Ancillary_9km_3dB','IGBP_9km_3dB','IGBP_Names','zlidar_9km_3dB','water_9km_3dB')
load('MODISTreeCoverEASE2_2016_9km_Africa')
TreeCoverEASE2_2016_9km=TreeCoverEASE2_2016_9km(12:1012,58:858);
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy')
% load('IncrementMedian_HalfDeg_4Inc_DryEnd_NullRemove_SMThreshEst_Globe')

load('IncrementMedian_HalfDeg_4Inc_DryEnd_NullRemove_SMThreshEstNoNull')
IncMedianDry = IncMedian;
% % PixCount(DryDownCount<10)=NaN;
IncMedian(IncCount<50)=NaN; 
% DryDownCount(DryDownCount<10)=NaN;
% IncAfDry = IncAf./0.11;
load('IncrementMedian_HalfDeg_4Inc_DryEnd_NullRemove_SMThreshEstNoNull')
IncMedianDry = IncMedian;
IncMedianDry(IncCount<50)=NaN;
load('IncrementMedian_HalfDeg_4Inc_WetEnd_NullRemove_SMThreshEstNoNull')
IncMedianDry(IncCount<50)=NaN;

% DryDownCount(DryDownCount<10)=NaN;
atrans = 0.35;
dll = 0.5; 
% AfricaRow = 300:1300;
% AfricaCol = 1700:2500;
WaterAfrica = water_9km_3dB;
ForestFlag = ForestFlag;
IGBPAfrica = IGBP_9km_3dB;
LidarAfrica = zlidar_9km_3dB;
load('SMThresholdEstimate_NullRemove_Globe')

% scatter(LidarAfrica(:),TreeCoverEASE2_2016_9km(:))
% xlabel('Lidar Height (m)')
% ylabel('Tree Cover (%)')

% IGBPWoodyFraction = IGBPAfrica(:,:,2)+IGBPAfrica(:,:,8);
%Remove peninsula
% IncAfMedianDry(1:6,120:end) = NaN;
% IncMedianWet(1:6,120:end) = NaN;
% IncAfMedianDry(7:11,123:end) = NaN;
% IncMedianWet(7:11,123:end) = NaN;

lon = SMAPCenterLongitudes;
lat = SMAPCenterLatitudes;  
latAve = nanmean(lat,2);
lonAve = nanmean(lon,1);

alatGlobe = [ -60 : dll : 60 ];                              
alonGlobe = [ -179 : dll : 179 ];  
nlatGlobe = length(alatGlobe);                                    
nlonGlobe = length(alonGlobe); 
alatGlobe = flipud(transpose(repmat(alatGlobe,nlonGlobe,1))); %lat matrix         
alonGlobe = repmat(alonGlobe,nlatGlobe,1);                  %lon matrix 

alatAf = [ -36 : dll : 38 ];                              
alonAf = [ -18 : dll : 52 ]; 
nlatAf = length(alatAf);                                    
nlonAf = length(alonAf); 
alatAf = flipud(transpose(repmat(alatAf,nlonAf,1))); %lat matrix         
alonAf = repmat(alonAf,nlatAf,1);                  %lon matrix
rowAf = find(alatGlobe(:,1)==alatAf(1)):find(alatGlobe(:,1)==alatAf(end));
colAf = find(alonGlobe(1,:)==alonAf(1)):find(alonGlobe(1,:)==alonAf(end));

alatSA = [ -57 : dll : 13 ];                              
alonSA = [ -82 : dll : -32 ]; 
nlatSA = length(alatSA);                                    
nlonSA = length(alonSA); 
alatSA = flipud(transpose(repmat(alatSA,nlonSA,1))); %lat matrix         
alonSA = repmat(alonSA,nlatSA,1);                  %lon matrix 
rowSA = find(alatGlobe(:,1)==alatSA(1)):find(alatGlobe(:,1)==alatSA(end));
colSA = find(alonGlobe(1,:)==alonSA(1)):find(alonGlobe(1,:)==alonSA(end));

alatAu = [ -45 : dll : -11 ];                              
alonAu = [ 111 : dll : 155 ]; 
nlatAu = length(alatAu);                                    
nlonAu = length(alonAu); 
alatAu = flipud(transpose(repmat(alatAu,nlonAu,1))); %lat matrix         
alonAu = repmat(alonAu,nlatAu,1);                  %lon matrix
rowAu = find(alatGlobe(:,1)==alatAu(1)):find(alatGlobe(:,1)==alatAu(end));
colAu = find(alonGlobe(1,:)==alonAu(1)):find(alonGlobe(1,:)==alonAu(end));        
%lon matrix     

Ncr = 100;
Ncb = 500;
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBu = [red;flipud(blu)];

% ForestFlagMat = NaN(nlat,nlon);
% IGBPMat = NaN(nlat,nlon,16);
% TreeCoverMat = NaN(nlat,nlon);
% WaterMat = NaN(nlat,nlon);
% RFIMat = NaN(nlat,nlon);
% MountMat = NaN(nlat,nlon);
load('GlobeAncillaryHalfDegree')
load('GlobeAncillaryHalfDegreeTwoYearGPM')
PrecipMat=PrecipMat./2;

load('GlobeAncillaryHalfDegreeFlags')
load('GlobeAncillaryHalfDegreeWater')
% WaterMat(WaterMat>0.05)=nan;
% WaterMat(WaterMat<=0.05)=1;
% BareFlag = IGBPMat(:,:,16);
% BareFlag(BareFlag==1)=nan;
% BareFlag(BareFlag==0)=1;

% TreeCoverMat1 = TreeCoverMat;
% TreeCoverMat1 = TreeCoverMat1.*RFIMat.*MountMat.*WaterMat.*BareFlag;
% TreeCoverMat1(TreeCoverMat1>50)=nan;

% %%% THIS IS ONLY FOR TROPICAL ANALYSIS!
% TreeCoverMat1(1:95,446:end)=nan;
% TreeCoverMat1(1:90,440:end)=nan;
% TreeCoverMat1(1:80,435:end)=nan;
% TreeCoverMat1(1:64,432:end)=nan;
% TreeCoverMat1(1:71,429:end)=nan;
% TreeCoverMat1(1:65,426:end)=nan;
% TreeCoverMat1(1:62,425:end)=nan;
% TreeCoverMat1(1:48,340:360)=nan;
% TreeCoverMat1(1:48,383:end)=nan;
% 
% TreeCoverMat1Af = TreeCoverMat1(rowAf,colAf);
% TreeCoverMat1SA = TreeCoverMat1(rowSA,colSA);
% TreeCoverMat1Au = TreeCoverMat1(rowAu,colAu);
% TreeCoverMat1Af(isfinite(TreeCoverMat1Af))=1;
% TreeCoverMat1SA(isfinite(TreeCoverMat1SA))=1;
% TreeCoverMat1Au(isfinite(TreeCoverMat1Au))=1;

GPMAf = PrecipMat(rowAf,colAf);
GPMSA = PrecipMat(rowSA,colSA);
GPMAu = PrecipMat(rowAu,colAu);

TreeCoverAf = TreeCoverMat(rowAf,colAf);
TreeCoverSA = TreeCoverMat(rowSA,colSA);
TreeCoverAu = TreeCoverMat(rowAu,colAu);

% AfTotalPix = nansum(TreeCoverMat1Af(:));
% SATotalPix = nansum(TreeCoverMat1SA(:));
% AuTotalPix = nansum(TreeCoverMat1Au(:));

IGBPMatAverageGlobe = NaN(nlatGlobe,nlonGlobe,16);
for i = 1:size(IGBPMatAverageGlobe,3)
    A = IGBPMat(:,:,i);
    A(A<0.5)=NaN;
    A(A>=0.5)=1;
    IGBPMatAverageGlobe(:,:,i) = A*i;
end
IGBPMatAverageGlobe = nanmean(IGBPMatAverageGlobe,3);
% PrecipMat = NaN(nlat,nlon);
% for ilat = 1 : nlat     
%     sprintf(['irow = ' sprintf('%0.2f',(ilat/nlat)) ])
%          for ilon = 1 : nlon   
% % Find SMAP EASE2 Rows/Cols
%     [DegRow] = find(latAve>= alat(ilat,ilon)-(dll/2)   & ...
%              latAve<= alat(ilat,ilon)+(dll/2));
%     [DegCol] = find(lonAve>= alon(ilat,ilon)-(dll/2)   & ...
%              lonAve<= alon(ilat,ilon)+(dll/2)   ); 
%          % Forest Mat  
%          ForestBlock = ForestFlag(DegRow,DegCol);
%            if nansum(ForestBlock(:))==0
%            ForestFlagMat(ilat,ilon) = 1;
%            else
%            ForestFlagMat(ilat,ilon) = nan;           
%            end
%          
           % Mountain Mat
%          A1 = MountainFlag(DegRow,DegCol);
%          MountMat(ilat,ilon) = nanmean(A1(:));
           
           % RFI Mat
%          A1 = RFIFlags(DegRow,DegCol);
%          RFIMat(ilat,ilon) = nanmean(A1(:));           
           
%          % Tree Cover Mat
%          A1 = TreeCoverEASE2_2016_9km(DegRow,DegCol);
%          TreeCoverMat(ilat,ilon) = nanmean(A1(:));
% 
%          A1 = WaterAfrica(DegRow,DegCol);
%          WaterMat(ilat,ilon) = nanmean(A1(:));
%          Precipitation
%          A2 = GPMMonSum_201504_201703(DegRow,DegCol);
%          PrecipMat(ilat,ilon) = nanmean(A2(:));
%          
%          % IGBP Mat
%           for iG = 1:16
%              A = IGBPAfrica(DegRow,DegCol,iG);
%              IGBPMat(ilat,ilon,iG) = nanmean(A(:));
%           end
%           
%          end
% end
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
% save -v7.3 GlobeAncillaryHalfDegreeFlags.mat MountMat RFIMat
%
%
% scatter(PrecipMat(:),SMThreshAf(:))
% IGBPDominantMat=NaN(nlat,nlon);
% for iG = 1:16
%  A = IGBPMat(:,:,iG);
%  A(A<0.5)=NaN;
%  A(A>=0.5)=iG;
% IGBPDominantMat = IGBPDominantMat.*A;
% end
% WaterMat(WaterMat>5)=nan;
% WaterMat(isfinite(WaterMat))=1;
% WaterMat(7:11,123:end) = NaN;
% WaterMat(1:6,120:end) = NaN;

ForestOnlyMat = ForestFlagMat;
ForestOnlyMat(ForestOnlyMat==1)=0;
ForestOnlyMat(isnan(ForestOnlyMat))=1;
ForestOnlyMat(ForestOnlyMat==0)=NaN;

IncMedianWetAf = IncMedianDry(rowAf,colAf);
IncMedianWetSA = IncMedianDry(rowSA,colSA);
IncMedianWetAu = IncMedianDry(rowAu,colAu);

IncMedianWetAf(1:51,124:end)=nan;
IncMedianWetAf(1:46,119:end)=nan;
IncMedianWetAf(1:36,114:end)=nan;
IncMedianWetAf(1:29,109:end)=nan;
IncMedianWetAf(1:4,1:35)=nan;
IncMedianWetAf(1:5,78:end)=nan;
IncMedianWetAf(1:3,63:67)=nan;
IncMedianWetAf(1:23,102:end)=nan;
IncMedianWetAf(1:20,94:end)=nan;

IncMedianWetNoForestAf = IncMedianWetAf.*ForestFlagMat(rowAf,colAf);
IncMedianWetForestAf = IncMedianWetAf.*ForestOnlyMat(rowAf,colAf);
IncMedianWetNoForestAu = IncMedianWetAu.*ForestFlagMat(rowAu,colAu);
IncMedianWetForestAu = IncMedianWetAu.*ForestOnlyMat(rowAu,colAu);
IncMedianWetNoForestSA = IncMedianWetSA.*ForestFlagMat(rowSA,colSA);
IncMedianWetForestSA = IncMedianWetSA.*ForestOnlyMat(rowSA,colSA);

% IncMedianWetNoForest = IncMedianWet.*ForestFlagMat;
% IncMedianWetForest = IncMedianWet.*ForestOnlyMat;

% IncMedianWetNoForestAf = IncMedianWetNoForest(rowAf,colAf);
% IncMedianWetForestAf = IncMedianWetForest(rowAf,colAf);
% IncMedianWetNoForestSA = IncMedianWetNoForest(rowSA,colSA);
% IncMedianWetForestSA = IncMedianWetForest(rowSA,colSA);
% IncMedianWetNoForestAu = IncMedianWetNoForest(rowAu,colAu);
% IncMedianWetForestAu = IncMedianWetForest(rowAu,colAu);
% 
% IncMedianWetNoForestAf(1:51,124:end)=nan;
% IncMedianWetNoForestAf(1:46,119:end)=nan;
% IncMedianWetNoForestAf(1:36,114:end)=nan;
% IncMedianWetNoForestAf(1:29,109:end)=nan;
% IncMedianWetNoForestAf(1:4,1:35)=nan;
% IncMedianWetNoForestAf(1:5,78:end)=nan;
% IncMedianWetNoForestAf(1:3,63:67)=nan;
% IncMedianWetNoForestAf(1:23,102:end)=nan;
% IncMedianWetNoForestAf(1:20,94:end)=nan;
% 
% IncMedianWetForestAf(1:51,124:end)=nan;
% IncMedianWetForestAf(1:46,119:end)=nan;
% IncMedianWetForestAf(1:36,114:end)=nan;
% IncMedianWetForestAf(1:29,109:end)=nan;
% IncMedianWetForestAf(1:4,1:35)=nan;
% IncMedianWetForestAf(1:5,78:end)=nan;
% IncMedianWetForestAf(1:3,63:67)=nan;
% IncMedianWetForestAf(1:23,102:end)=nan;
% IncMedianWetForestAf(1:20,94:end)=nan;

IncMedianWetAllAf = IncMedianDry(rowAf,colAf);
IncMedianWetAllSA = IncMedianDry(rowSA,colSA);
IncMedianWetAllAu = IncMedianDry(rowAu,colAu);

IncMedianWetAllAf(1:51,124:end)=nan;
IncMedianWetAllAf(1:46,119:end)=nan;
IncMedianWetAllAf(1:36,114:end)=nan;
IncMedianWetAllAf(1:29,109:end)=nan;
IncMedianWetAllAf(1:4,1:35)=nan;
IncMedianWetAllAf(1:5,78:end)=nan;
IncMedianWetAllAf(1:3,63:67)=nan;
IncMedianWetAllAf(1:23,102:end)=nan;
IncMedianWetAllAf(1:20,94:end)=nan;

% % Estimating number of pixels with pulse response
% IncMedianWetNoForestAf1 = IncMedianWetNoForestAf.*TreeCoverMat1Af;
% IncMedianWetNoForestAf1(IncMedianWetNoForestAf1>0)=nan;
% IncMedianWetNoForestAf1(isfinite(IncMedianWetNoForestAf1))=1;
% AfPulsePix = nansum(IncMedianWetNoForestAf1(:));
% AfPulsePer = AfPulsePix/AfTotalPix;
% 
% IncMedianWetNoForestSA1 = IncMedianWetNoForestSA.*TreeCoverMat1SA;
% IncMedianWetNoForestSA1(IncMedianWetNoForestSA1>0)=nan;
% IncMedianWetNoForestSA1(isfinite(IncMedianWetNoForestSA1))=1;
% SAPulsePix = nansum(IncMedianWetNoForestSA1(:));
% SAPulsePer = SAPulsePix/SATotalPix;
% 
% IncMedianWetNoForestAu1 = IncMedianWetNoForestAu.*TreeCoverMat1Au;
% IncMedianWetNoForestAu1(IncMedianWetNoForestAu1>0)=nan;
% IncMedianWetNoForestAu1(isfinite(IncMedianWetNoForestAu1))=1;
% AuPulsePix = nansum(IncMedianWetNoForestAu1(:));
% AuPulsePer = AuPulsePix/AuTotalPix;
% TotalPulsePer = (AfPulsePix+SAPulsePix+AuPulsePix)/(AuTotalPix+SATotalPix+AfTotalPix);



[~,~,~,~,...
   TopLeftCornerLatAf,~,~,~] = InterpEdgesCorners(alatAf);
TopLeftCornerLatAf(:,1) = TopLeftCornerLatAf(:,2);
TopLeftCornerLatAf(1,:) = TopLeftCornerLatAf(2,:)+0.5;

[~,~,~,~,...
   TopLeftCornerLonAf,~,~,~] = InterpEdgesCorners(alonAf);
TopLeftCornerLonAf(1,:) = TopLeftCornerLonAf(2,:);
TopLeftCornerLonAf(:,1) = TopLeftCornerLonAf(:,2)-0.5;

[~,~,~,~,...
   TopLeftCornerLatSA,~,~,~] = InterpEdgesCorners(alatSA);
TopLeftCornerLatSA(:,1) = TopLeftCornerLatSA(:,2);
TopLeftCornerLatSA(1,:) = TopLeftCornerLatSA(2,:)+0.5;

[~,~,~,~,...
   TopLeftCornerLonSA,~,~,~] = InterpEdgesCorners(alonSA);
TopLeftCornerLonSA(1,:) = TopLeftCornerLonSA(2,:);
TopLeftCornerLonSA(:,1) = TopLeftCornerLonSA(:,2)-0.5;

[~,~,~,~,...
   TopLeftCornerLatAu,~,~,~] = InterpEdgesCorners(alatAu);
TopLeftCornerLatAu(:,1) = TopLeftCornerLatAu(:,2);
TopLeftCornerLatAu(1,:) = TopLeftCornerLatAu(2,:)+0.5;

[~,~,~,~,...
   TopLeftCornerLonAu,~,~,~] = InterpEdgesCorners(alonAu);
TopLeftCornerLonAu(1,:) = TopLeftCornerLonAu(2,:);
TopLeftCornerLonAu(:,1) = TopLeftCornerLonAu(:,2)-0.5;

SizeSA = size(IncMedianWetForestSA);
SizeAf = size(IncMedianWetForestAf);
SizeAu = size(IncMedianWetForestAu);

% Total Area from water fraction 
WaterMatAf = WaterMat(rowAf,colAf);
WaterMatAu = WaterMat(rowAu,colAu);
WaterMatSA = WaterMat(rowSA,colSA);
WaterMatAf(WaterMatAf>0.05)=nan;
WaterMatAf(WaterMatAf<=0.05)=1;
WaterMatAu(WaterMatAu>0.05)=nan;
WaterMatAu(WaterMatAu<=0.05)=1;
WaterMatSA(WaterMatSA>0.05)=nan;
WaterMatSA(WaterMatSA<=0.05)=1;
WaterMatAf(1:51,124:end)=nan;
WaterMatAf(1:46,119:end)=nan;
WaterMatAf(1:36,114:end)=nan;
WaterMatAf(1:29,109:end)=nan;
WaterMatAf(1:4,1:35)=nan;
WaterMatAf(1:5,78:end)=nan;
WaterMatAf(1:3,63:67)=nan;
WaterMatAf(1:23,102:end)=nan;
WaterMatAf(1:20,94:end)=nan;

IGBPAf = IGBPMatAverageGlobe(rowAf,colAf);
IGBPAu = IGBPMatAverageGlobe(rowAu,colAu);
IGBPSA = IGBPMatAverageGlobe(rowSA,colSA);
IGBPAfNonVegMask = IGBPAf;
IGBPAfNonVegMask(IGBPAfNonVegMask==11|IGBPAfNonVegMask==13|IGBPAfNonVegMask==15|IGBPAfNonVegMask==16)=nan;
IGBPAfNonVegMask(isfinite(IGBPAfNonVegMask))=1;
IGBPAuNonVegMask = IGBPAu;
IGBPAuNonVegMask(IGBPAuNonVegMask==11|IGBPAuNonVegMask==13|IGBPAuNonVegMask==15|IGBPAuNonVegMask==16)=nan;
IGBPAuNonVegMask(isfinite(IGBPAuNonVegMask))=1;
IGBPSANonVegMask = IGBPSA;
IGBPSANonVegMask(IGBPSANonVegMask==11|IGBPSANonVegMask==13|IGBPSANonVegMask==15|IGBPSANonVegMask==16)=nan;
IGBPSANonVegMask(isfinite(IGBPSANonVegMask))=1;

WaterMatAf = WaterMatAf.*IGBPAfNonVegMask;
WaterMatAu = WaterMatAu.*IGBPAuNonVegMask;
WaterMatSA = WaterMatSA.*IGBPSANonVegMask;

Total = nansum([WaterMatAf(:); WaterMatAu(:); WaterMatSA(:)]);

% Only consider the non-veg regions
TreeCoverAfCount = TreeCoverAf.*WaterMatAf;
TreeCoverAuCount = TreeCoverAu.*WaterMatAu;
TreeCoverSACount = TreeCoverSA.*WaterMatSA;
GPMAfCount = GPMAf.*WaterMatAf;
GPMAuCount = GPMAu.*WaterMatAu;
GPMSACount = GPMSA.*WaterMatSA;
IncMedianWetAllAfCount = IncMedianWetAllAf.*WaterMatAf;
IncMedianWetAllAuCount = IncMedianWetAllAu.*WaterMatAu;
IncMedianWetAllSACount = IncMedianWetAllSA.*WaterMatSA;

% IGBP
IGBPMatAf = IGBPMat(rowAf,colAf,:);
IGBPMatSA = IGBPMat(rowSA,colSA,:);
IGBPMatAu = IGBPMat(rowAu,colAu,:);
XticksMod = [7,10,14,9,8];
WetIncIGBP = nan(5,4000);
IncIGBPConf = nan(size(XticksMod,2),3);
for j = 1:size(XticksMod,2)
    iB = XticksMod(j);
    IncDiff1Af = IncMedianWetAf;
    IncDiff1SA = IncMedianWetSA;
    IncDiff1Au = IncMedianWetAu;
    
    IGBP1Af = IGBPMatAf(:,:,iB);
    IGBP1Af(IGBP1Af<0.5)=NaN;
    IGBP1Af(IGBP1Af>=0.5)=1;   
    IGBP1SA = IGBPMatSA(:,:,iB);
    IGBP1SA(IGBP1SA<0.5)=NaN;
    IGBP1SA(IGBP1SA>=0.5)=1;    
    IGBP1Au = IGBPMatAu(:,:,iB);
    IGBP1Au(IGBP1Au<0.5)=NaN;
    IGBP1Au(IGBP1Au>=0.5)=1;     
    
    IncDiff1 = [IncDiff1Af(:).*IGBP1Af(:); IncDiff1SA(:).*IGBP1SA(:); ...
        IncDiff1Au(:).*IGBP1Au(:)];
    IncTotalCount = [WaterMatAf(:).*IGBP1Af(:); WaterMatSA(:).*IGBP1SA(:); ...
        WaterMatAu(:).*IGBP1Au(:)];
    Boot90IGBPCount = [Boot90Af(:).*IGBP1Af(:); Boot90SA(:).*IGBP1SA(:); ...
        Boot90Au(:).*IGBP1Au(:)];
    
    % Remove non pulse response due to bootstrap confidence
    Boot90IGBPCount(Boot90IGBPCount>0)=nan;
    
    % Remove non veg and outside region pixels
    Boot90IGBPCount(isnan(IncTotalCount))=nan;
    
    % Compute Total and Pulse Response
    Boot90IGBPCount(isfinite(Boot90IGBPCount))=1;
    IncTotalCount(isfinite(IncTotalCount))=1;    
    IncIGBPConf(j,2) = nansum(Boot90IGBPCount);
    IncIGBPConf(j,3) = nansum(IncTotalCount);
    IncIGBPConf(j,1) = nansum(Boot90IGBPCount)/nansum(IncTotalCount);
    
    IncDiff1(isnan(IncDiff1))=[];
    WetIncIGBP(j,1:size(IncDiff1,1)) = IncDiff1';
end
IncIGBPConf(:,1) = IncIGBPConf(:,1)*100;
% PulseBehaviorTotal = [IncMedianWetAllAf(:); IncMedianWetAllAu(:); IncMedianWetAllSA(:)];
% PulseBehaviorTotal(isnan(PulseBehaviorTotal))=[];
% PulseBehaviorTotal(isfinite(PulseBehaviorTotal))=1;
% PulseResponsePercentage = nansum(PulseBehaviorTotal(:))/Total;

%%
figure
set(gcf,'position',[100 100 800 800])
set(gcf,'color','white')
%%%%%%%%%%%%%%%%%% South America %%%%%%%%%%%%%%%%%%%%%%
MSA = subplot(4,8,[3 5 11 13]);
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,IncMedianWetForestSA)
alpha(atrans)
hold on
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,IncMedianWetNoForestSA)
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% colorbar
caxis([-5 25])
colormap(cmapRedBu)
set(gca,'color',[0.7 0.7 0.7])
caxis([-5 25])
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',15)


% 1) Set size of first tile
AMSA = get(MSA,'Position'); % top right fig
set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(3)-0.05 AMSA(4)]);
AMSA = get(MSA,'Position'); % top right fig

% 2) Set aspect ratio
% set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(4)*(SizeSA(2)/SizeSA(1)) AMSA(4)]);
set(MSA,'Position',[AMSA(1)-0.02 AMSA(2)+0.085 AMSA(3) AMSA(3)*(SizeSA(1)/SizeSA(2))]);
AMSA = get(MSA,'Position'); % top right fig

%%%%%%%%%%%%%%%%%%%% Africa %%%%%%%%%%%%%%%%%%%%%%%
% MAf = subplot(3,2,[2 4]);
% figure
MAf = subplot(4,8,[6 8 14 16]);
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,IncMedianWetForestAf)
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,IncMedianWetNoForestAf)
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
set(gca,'ytick',{})
set(gca,'xtick',{})
c = colorbar;
c.Label.String = '\DeltaVWC/\DeltaSM  [(kg m^{-2})/(m^3 m^{-3})]';
c.Label.Rotation = 270;
pos = get(c,'Position')
set(gca,'color',[0.7 0.7 0.7])
c.Label.Position = [pos(1)+3 pos(2)+9];
shading flat
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
caxis([-5 25])
colormap(cmapRedBu)

% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',17)

AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1)-0.03 AMAf(2)-0.035 AMAf(3)*(SizeAf(2)/SizeSA(2)) AMAf(4)]);
AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1) AMAf(2)+0.11 AMAf(3) AMAf(3)*SizeAf(1)/SizeAf(2)]);
AMAf = get(MAf,'Position'); 

%%%%%%%%%%%%%%%%%%%   Australia %%%%%%%%%%%%%%%%%%%%%
% figure
MAu = subplot(4,8,[1 2 9 10]);
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,IncMedianWetForestAu)
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,IncMedianWetNoForestAu)
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% colorbar
% colormap(SMAP_Color_SoilMoisture./300)
set(gca,'color',[0.7 0.7 0.7])
caxis([-5 25])
colormap(cmapRedBu)
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',15)

AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1) AMAu(2) AMSA(3)*(SizeAu(2)/SizeSA(2)) AMAu(4)]);
AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1)-0.08 AMAu(2)+0.18 AMAu(3) AMAu(3)*SizeAu(1)/SizeAu(2)]);
AMAu = get(MAu,'Position'); 


%%%%%%%%%%%%%%%%%%%   IGBP %%%%%%%%%%%%%%%%%%%%%%%
Xticks = 1:16;
% XticksMod = [7,8,9,10,14];
XticksMod = [7,10,14,9,8];
Xticks1Mod = IGBP_Names(XticksMod);
dVWCmatIGBPMod = WetIncIGBP;
dVWCmatIGBPModCount = dVWCmatIGBPMod;
dVWCmatIGBPModCount(isfinite(dVWCmatIGBPModCount))=1;
IGBPBinCount = nansum(dVWCmatIGBPModCount,2);

IGBPBox = subplot(4,8,[17 20]);
boxplot(WetIncIGBP','notch','off')
h=findobj(gca,'tag','Outliers');
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'y','FaceAlpha',0.2);
end
hold on
A = 0:0.01:35;
plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1Mod,'fontsize',15)
xtickangle(25)
ylim([-5 25])
ylabel({'\DeltaVWC/\DeltaSM'; '[(kg m^{-2})/(m^3 m^{-3})]'})
grid on
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines,'linewidth',2);
set(gca,'fontname','arial')

AIGBP = get(IGBPBox,'Position'); 
set(IGBPBox,'Position',[AIGBP(1)-0.04 AIGBP(2)-0.01 AIGBP(3)+0.01 AIGBP(4)+0.14]);
AIGBP = get(IGBPBox,'Position'); 

% for ia = 1:size(IncIGBPConf,1)
%     aNum = [sprintf('%0.0f',IncIGBPConf(ia,1)) '%'];
% annotation('textbox',[0.15+((ia-1)/100)*6.8 0.51 0.1 0.1],'String',aNum,'fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
% end

% IncMedianWetNoForestAll = [IncMedianWetNoForestAf(:); IncMedianWetNoForestAu(:); ...
%                             IncMedianWetNoForestSA(:)];
IncMedianWetBoxAll = [IncMedianWetAf(:); IncMedianWetAu(:); ...
                            IncMedianWetSA(:)];                        
TreeCoverMatAf = TreeCoverMat(rowAf,colAf);
TreeCoverMatSA = TreeCoverMat(rowSA,colSA);
TreeCoverMatAu = TreeCoverMat(rowAu,colAu);
TreeCoverMatAll = [TreeCoverMatAf(:); TreeCoverMatAu(:); TreeCoverMatSA(:)];

IncMedianTreeCover = [TreeCoverMatAll,IncMedianWetBoxAll];
% scatter(IncMedianWetNoForestAll(:),TreeCoverMatAll(:))
IncMedianTreeCover(any(isnan(IncMedianTreeCover),2),:)=[];
[corrTreeFrac corrpTreeFrac]=corrcoef(IncMedianTreeCover(:,1),IncMedianTreeCover(:,2));

%%%%%%%%%%%%%%%% Tree Cover Prevalence %%%%%%%%%%%%%%%%%
edge = [0 5 10 20 30 50 100];
% GPMAll = [GPMAf(:); GPMAu(:); GPMSA(:)];
TreeCoverMatAllCount = [TreeCoverMatAf(:); TreeCoverMatAu(:); TreeCoverMatSA(:)];
% GPMAllCount = [GPMAfCount(:); GPMAuCount(:); GPMSACount(:)];
Boot90 = [Boot90Af(:); Boot90Au(:); Boot90SA(:)];
TotalCount = [WaterMatAf(:); WaterMatAu(:); WaterMatSA(:)];
IncTreeCoverConf = nan(size(edge,2)-1,3);
for i = 1:size(edge,2)-1
   TreeCoverBin = TreeCoverMatAllCount; 
   Inc90Bin = Boot90;
   TotalBin = TotalCount;
   % 1) get rid of values outside of veg areas
   TreeCoverBin(isnan(TotalBin)) = nan;
   Inc90Bin(isnan(TotalBin)) = nan;
   % 2) Get rid of non-pulse response pixels 
   Inc90Bin(Inc90Bin>0) = nan;
   % 3) Find only confident pulse response and total non veg in GPM bin
   Inc90Bin(TreeCoverBin<=edge(i) | TreeCoverBin>edge(i+1)) = nan;
   TreeCoverBin(TreeCoverBin<=edge(i) | TreeCoverBin>edge(i+1)) = nan;
   % 4) Count them up!
   Inc90Bin(isfinite(Inc90Bin))=1;
   TreeCoverBin(isfinite(TreeCoverBin))=1;
   IncTreeCoverConf(i,1) = nansum(Inc90Bin)/nansum(TreeCoverBin);
   IncTreeCoverConf(i,2) = nansum(Inc90Bin);
   IncTreeCoverConf(i,3) = nansum(TreeCoverBin);
   
end
IncTreeCoverConf(:,1) = IncTreeCoverConf(:,1).*100;
% nansum(IncTreeCoverConf(:,3))
% for ia = 1:size(IncTreeCoverConf,1)-1
%     aNum = [sprintf('%0.0f',IncTreeCoverConf(ia,1)) '%'];
% annotation('textbox',[0.63+((ia-1)/100)*6.8 0.51 0.1 0.1],'String',aNum,'fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
% end


Xticks = [1 2 3 4 5 6 7 8];
Xticks1 = {'0' '0-5' '5-10' '10-20' '20-30' '30-50' '50-100'};
% Xticks1 = {'0-5' '5-15' '15-30' '30-50' '50-75'};
TreeCoverBox = subplot(4,8,[21 24]);
edge = [0 5 10 20 30 50 100];
[n,bin] = histc(TreeCoverMatAll(:),edge);
boxplot(IncMedianWetBoxAll(:),bin,'notch','off');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',0.2);
end
h=findobj(gca,'tag','Outliers');
xlim([1.5 6.5])
delete(h)
% hold on
grid on
hold on
ylabel({'\DeltaVWC/\DeltaSM'; '[(kg m^{-2})/(m^3 m^{-3})]'})

A = 0:0.01:35;
plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
xlabel('Tree Cover (%)')
set(gca,'fontsize',15)
set(gca,'fontname','arial')
ylim([-5 25])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines,'linewidth',2);
ATreeCover = get(TreeCoverBox,'Position'); 
set(TreeCoverBox,'Position',[ATreeCover(1)+0.03 AIGBP(2)...
            AIGBP(3) AIGBP(4)]);
      
% Scale up all continents
ScaleFactor = 1;
AMSA = get(MSA,'Position'); % top right fig
AMAf = get(MAf,'Position'); % top right fig
AMAu = get(MAu,'Position'); % top right fig
set(MAf,'Position',[AMAf(1) AMAf(2)+0.03 AMAf(3)*ScaleFactor AMAf(4)*ScaleFactor]);   
set(MSA,'Position',[AMSA(1) AMSA(2)+0.03 AMSA(3)*ScaleFactor AMSA(4)*ScaleFactor]);
set(MAu,'Position',[AMAu(1) AMAu(2)+0.03 AMAu(3)*ScaleFactor AMAu(4)*ScaleFactor]);

edge = [0 5 10 20 30 50 100];
BoxSamplenTreeCover= nan(1,size(edge,2)-1);
for i = 1:size(edge,2)-1
    TreeCoverMatn = TreeCoverMatAll;
    TreeCoverMatn(TreeCoverMatn<edge(i))=nan;
    TreeCoverMatn(TreeCoverMatn>edge(i+1))=nan;    
    TreeCoverMatn(isfinite(TreeCoverMatn))=1;
    IncAfMedianWetn = IncMedianWetBoxAll.*TreeCoverMatn;
    IncAfMedianWetn(isfinite(IncAfMedianWetn))=1;  
    BoxSamplenTreeCover(i) = nansum(IncAfMedianWetn(:));
end

%%%%%%%%%%%%%%%%%%  GPM Prevalence Values %%%%%%%%%%%%%%%%%
% edge = [0 500 1000 1500 2000 2500 3000 10000];
edge = [0 250 500 750 1000 1500 2000 10000];
GPMAll = [GPMAf(:); GPMAu(:); GPMSA(:)];
GPMAllCount = [GPMAfCount(:); GPMAuCount(:); GPMSACount(:)];
Boot90 = [Boot90Af(:); Boot90Au(:); Boot90SA(:)];
% IncAllCount = [IncMedianWetAllAfCount(:); IncMedianWetAllAuCount(:);...
%     IncMedianWetAllSACount(:)];
TotalCount = [WaterMatAf(:); WaterMatAu(:); WaterMatSA(:)];
IncGPMConf = nan(size(edge,2)-1,3);
for i = 1:size(edge,2)-1
   GPMBin = GPMAllCount; 
   Inc90Bin = Boot90;
   TotalBin = TotalCount;
   
   % 1) get rid of values outside of veg areas
   GPMBin(isnan(TotalBin)) = nan;
   Inc90Bin(isnan(TotalBin)) = nan;
   % 2) Get rid of non-pulse response pixels 
   Inc90Bin(Inc90Bin>0) = nan;
   % 3) Find only confident pulse response and total non veg in GPM bin
   Inc90Bin(GPMBin<edge(i) | GPMBin>edge(i+1)) = nan;
   GPMBin(GPMBin<edge(i) | GPMBin>edge(i+1)) = nan;
   % 4) Count them up!
   Inc90Bin(isfinite(Inc90Bin))=1;
   GPMBin(isfinite(GPMBin))=1;
   IncGPMConf(i,1) = nansum(Inc90Bin)/nansum(GPMBin);
   IncGPMConf(i,2) = nansum(Inc90Bin);
   IncGPMConf(i,3) = nansum(GPMBin);  
end
IncGPMConf(:,1) = IncGPMConf(:,1).*100;

BoxSamplenGPM= nan(1,size(edge,2)-1);
for i = 1:size(edge,2)-1
    GPMMatn = GPMAll;
    GPMMatn(GPMMatn<edge(i))=nan;
    GPMMatn(GPMMatn>edge(i+1))=nan;    
    GPMMatn(isfinite(GPMMatn))=1;
    IncAfMedianWetn = IncMedianWetBoxAll.*GPMMatn;
    IncAfMedianWetn(isfinite(IncAfMedianWetn))=1;  
    BoxSamplenGPM(i) = nansum(IncAfMedianWetn(:));
end

% Precipitation
Xticks = [1 2 3 4 5 6 7];
% Xticks1 = {'0' '0-500' '500-1000' '1000-1500' '1500-2000' '2000-2500' '2500-3000'};
Xticks1 = {'0' '0-250' '250-500' '500-750' '750-1000' '1000-1500' '1500-2000'};
% Xticks1 = {'0-5' '5-15' '15-30' '30-50' '50-75'};
PrecipBox = subplot(4,8,[27 30]);
[n,bin] = histc(GPMAll(:),edge);
boxplot(IncMedianWetBoxAll(:),bin,'notch','off');
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'b','FaceAlpha',0.2);
end
h=findobj(gca,'tag','Outliers');
xlim([1.5 7.5])
delete(h)
xtickangle(25)
% hold on
grid on
hold on
ylabel({'\DeltaVWC/\DeltaSM'; '[(kg m^{-2})/(m^3 m^{-3})]'})
A = 0:0.01:35;
plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
xlabel('Annual Precipitation (mm/year)')
set(gca,'fontsize',15)
set(gca,'fontname','arial')
ylim([-5 25])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines,'linewidth',2);
AGPM = get(PrecipBox,'Position'); 
AIGBP = get(IGBPBox,'Position'); 
set(PrecipBox,'Position',[AGPM(1)-0.1 AGPM(2)-0.07 AIGBP(3)+0.15 AIGBP(4)]);

% for ia = 1:size(IncGPMConf,1)-1
%     aNum = [sprintf('%0.0f',IncGPMConf(ia,1)) '%'];
% annotation('textbox',[0.31+((ia-1)/100)*8.05 0.19 0.1 0.1],'String',aNum,'fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
% end

annotation('textbox',[0.02 0.9 0.1 0.1],'String','a','fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.27 0.9 0.1 0.1],'String','b','fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.55 0.9 0.1 0.1],'String','c','fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.03 0.54 0.1 0.1],'String','d','fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.51 0.54 0.1 0.1],'String','e','fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.18 0.22 0.1 0.1],'String','f','fontsize',15,'Fontname','arial','FontWeight','bold','EdgeColor','none')
fig = gcf;
fig.InvertHardcopy = 'off';
cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
% print('FigS5DryResponse','-djpeg','-r300')
% print('FigS5DryResponse_100','-djpeg','-r100')
print('FigS5','-depsc','-r300') % save as high res PNGd    
