clc
clear
load coast   
clat = lat;                                                        
clon = long ; 
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP Data/General Data')
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\General Data')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
load('ForestFlag')
load('SMAPCenterCoordinates9KM')
load('soil_texture_9km')
load('Ancillary_9km_3dB','IGBP_9km_3dB','IGBP_Names')
load('MTDCA2_v13','SMAP_Color_SoilMoisture','SMAP_Color_Scale')
load('MODISTreeCoverEASE2_2016_9km_Africa')
TreeCoverEASE2_2016_9km=TreeCoverEASE2_2016_9km(12:1012,58:858);

cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy')
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\MIT\Manuscripts\WaterExchangeNatureGeoscience')
load('SMThresholdEstimate_NullRemove_Globe')
% load('SMThresholdEstimate_NullRemove_Globe_SMDeseason')
load('IncrementMedian_HalfDeg_4Inc_DryEnd_NullRemove_SMThreshEstNoNull')
SMThresh(IncCount<50)=nan;
load('IncrementMedian_HalfDeg_4Inc_WetEnd_NullRemove_SMThreshEstNoNull')
SMThresh(IncCount<50)=nan;

SMThresh(1:95,446:end)=nan;
SMThresh(1:90,440:end)=nan;
SMThresh(1:80,435:end)=nan;
SMThresh(1:64,432:end)=nan;
SMThresh(1:71,429:end)=nan;
SMThresh(1:65,426:end)=nan;
SMThresh(1:62,425:end)=nan;
SMThresh(1:48,340:360)=nan;
SMThresh(1:48,383:end)=nan;

% IncAfMedianDry = IncAfMedian;
% PixCount(DryDownCount<10)=NaN;
% IncAfMedian(IncCount<50)=NaN; 
% DryDownCount(DryDownCount<10)=NaN;

% IncAfDry = IncAf./0.11;
% % load('IncrementMedian_HalfDeg_4Inc_WetEnd_NullRemove')
% IncAfMedianWet = IncAfMedian;

% PixCount(DryDownCount<10)=NaN;
% IncAfMedianWet(IncCount<50)=NaN;
% DryDownCount(DryDownCount<10)=NaN;
atrans = 0.4;
dll = 0.5; 
AfricaRow = 300:1300;
AfricaCol = 1700:2500;
ForestFlag = ForestFlag(AfricaRow,AfricaCol);
IGBPAfrica = IGBP_9km_3dB(AfricaRow,AfricaCol,:);
ClayAfrica = clay(AfricaRow,AfricaCol,:);
SandAfrica = sand(AfricaRow,AfricaCol,:);
%Remove peninsula
% IncAfMedianDry(1:6,120:end) = NaN;
% IncAfMedianWet(1:6,120:end) = NaN;
% IncAfMedianDry(7:11,123:end) = NaN;
% IncAfMedianWet(7:11,123:end) = NaN;

lon = SMAPCenterLongitudes(AfricaRow,AfricaCol);
lat = SMAPCenterLatitudes(AfricaRow,AfricaCol);  
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


% ForestFlagMat = NaN(nlat,nlon);
% IGBPMat = NaN(nlat,nlon,16);
% ClayMat = NaN(nlat,nlon);
% SandMat = NaN(nlat,nlon);
% TreeCoverMat = NaN(nlat,nlon);

load('GlobeAncillaryHalfDegree')
load('GlobeAncillaryHalfDegreeTwoYearGPM')
load('GlobeAncillaryHalfDegreeFlags')
% for ilat = 1 : nlat     
%     sprintf(['irow = ' sprintf('%0.2f',(ilat/nlat)) ])
%          for ilon = 1 : nlon   
% % Find SMAP EASE2 Rows/Cols
%     [DegRow] = find(latAve>= alat(ilat,ilon)-(dll/2)   & ...
%              latAve<= alat(ilat,ilon)+(dll/2));
%     [DegCol] = find(lonAve>= alon(ilat,ilon)-(dll/2)   & ...
%              lonAve<= alon(ilat,ilon)+(dll/2)   ); 
%            ForestBlock = ForestFlag(DegRow,DegCol);
%            if nansum(ForestBlock(:))==0
%            ForestFlagMat(ilat,ilon) = 1;
%            else
%            ForestFlagMat(ilat,ilon) = nan;           
%            end
%            
%          %Clay fraction
%          A1 = ClayAfrica(DegRow,DegCol);
%           ClayMat(ilat,ilon) = nanmean(A1(:));
%          %Sand fraction
%          A1 = SandAfrica(DegRow,DegCol);
%           SandMat(ilat,ilon) = nanmean(A1(:));          
%          A1 = TreeCoverEASE2_2016_9km(DegRow,DegCol);
%          TreeCoverMat(ilat,ilon) = nanmean(A1(:));         
%         for iG = 1:16
%              A = IGBPAfrica(DegRow,DegCol,iG);
%              IGBPMat(ilat,ilon,iG) = nanmean(A(:));
%         end  
%         end        
% end

IGBPMatAverageGlobe = NaN(nlatGlobe,nlonGlobe,16);
for i = 1:size(IGBPMatAverageGlobe,3)
    A = IGBPMat(:,:,i);
    A(A<0.5)=NaN;
    A(A>=0.5)=1;
    IGBPMatAverageGlobe(:,:,i) = A*i;
end
IGBPMatAverageGlobe = nanmean(IGBPMatAverageGlobe,3);

% Included IGBP Classes and Regions
cIGBP = [  14  10   9  8  7  2  ];                        
IGBP_latlon(14,:) = [   9  15  -9  9 ];                            
IGBP_latlon(10,:) = [  10  17 -10 36 ];                           
IGBP_latlon( 9,:) = [ -22 -15 17 33  ];                            
IGBP_latlon( 8,:) = [ -15  -4 14 30  ];                            
IGBP_latlon( 7,:) = [ -32 -22 17 25  ];                            
IGBP_latlon( 2,:) = [  -4   5 10 30  ];

% IGBP Map
% figure
IGBPx = NaN(size(alatAf));                              
for i = 1 : length(cIGBP)
IGBPz = NaN(size(alatAf));   
%Finding pixels within bounds
dlatlon = find(alatAf>=IGBP_latlon(cIGBP(i),1) & ...
               alatAf<=IGBP_latlon(cIGBP(i),2) & ...
               alonAf>=IGBP_latlon(cIGBP(i),3) & ...
               alonAf<=IGBP_latlon(cIGBP(i),4) & ...
               IGBPMatAverageGlobe(rowAf,colAf)==cIGBP(i));
IGBPx(dlatlon) = cIGBP(i);                                             
IGBPz(dlatlon) = cIGBP(i);                                                     
end

% Matric Potential Conversion
% MatricPotMat = ((SMThreshAfForest)*(9810/10^6))
SatMatSand = -12.1;
SatMatClay = -40.5;
SatMatLoam = -47.8;
bSand = 4.05;
bClay = 11.4;
bLoam = 5.39;
nSand = 0.395;
nClay = 0.482;
nLoam = 0.451;

ForestOnlyMat = ForestFlagMat;
ForestOnlyMat(ForestOnlyMat==1)=0;
ForestOnlyMat(isnan(ForestOnlyMat))=1;
ForestOnlyMat(ForestOnlyMat==0)=NaN;

LoamMat = 1-ClayMat-SandMat;
beffMat = bSand.*SandMat+bClay.*ClayMat+bLoam.*LoamMat;
neffMat = nSand.*SandMat+nClay.*ClayMat+nLoam.*LoamMat;
SateffMat = SatMatSand.*SandMat+SatMatClay.*ClayMat+SatMatLoam.*LoamMat;
MatricMat = ((SateffMat.*(SMThresh./neffMat).^(-beffMat))/100).*(9810/10^6);

MatricMatNoForest = MatricMat.*ForestFlagMat;
MatricMatForest = MatricMat.*ForestOnlyMat;

SMThreshNoForest = SMThresh.*ForestFlagMat;
SMThreshForest= SMThresh.*ForestOnlyMat;

MatricMatAf = MatricMat(rowAf,colAf);
MatricMatSA = MatricMat(rowSA,colSA);
MatricMatAu = MatricMat(rowAu,colAu);

SMThreshAf = SMThresh(rowAf,colAf);
SMThreshSA = SMThresh(rowSA,colSA);
SMThreshAu = SMThresh(rowAu,colAu);

MatricMatNoForestAf = MatricMatNoForest(rowAf,colAf);
MatricMatForestAf = MatricMatForest(rowAf,colAf);
MatricMatNoForestSA = MatricMatNoForest(rowSA,colSA);
MatricMatForestSA = MatricMatForest(rowSA,colSA);
MatricMatNoForestAu = MatricMatNoForest(rowAu,colAu);
MatricMatForestAu = MatricMatForest(rowAu,colAu);

SMThreshNoForestAf = SMThreshNoForest(rowAf,colAf);
SMThreshForestAf = SMThreshForest(rowAf,colAf);
SMThreshNoForestSA = SMThreshNoForest(rowSA,colSA);
SMThreshForestSA = SMThreshForest(rowSA,colSA);
SMThreshNoForestAu = SMThreshNoForest(rowAu,colAu);
SMThreshForestAu = SMThreshForest(rowAu,colAu);

SandMatAf = SandMat(rowAf,colAf);
ClayMatAf = ClayMat(rowAf,colAf);
SandMatSA = SandMat(rowSA,colSA);
ClayMatSA = ClayMat(rowSA,colSA);
SandMatAu = SandMat(rowAu,colAu);
ClayMatAu = ClayMat(rowAu,colAu);

TreeCoverMatAf = TreeCoverMat(rowAf,colAf);
TreeCoverMatSA = TreeCoverMat(rowSA,colSA);
TreeCoverMatAu = TreeCoverMat(rowAu,colAu);
TreeCoverMatAll = [TreeCoverMatAf(:); TreeCoverMatSA(:); TreeCoverMatAu(:)];

GPMAf = PrecipMat(rowAf,colAf);
GPMSA = PrecipMat(rowSA,colSA);
GPMAu = PrecipMat(rowAu,colAu);
GPMAll = [GPMAf(:); GPMSA(:); GPMAu(:)];

SMThreshAll = [SMThreshAf(:); SMThreshSA(:);...
    SMThreshAu(:)];
MatricMatAll = [MatricMatAf(:); MatricMatSA(:);...
    MatricMatAu(:)];
ClayAll = [ClayMatAf(:);ClayMatSA(:);ClayMatAu(:)];
SandAll = [SandMatAf(:);SandMatSA(:);SandMatAu(:)];

SandVect = [SandAll(:),SMThreshAll(:)];
SandVect(any(isnan(SandVect),2),:)=[];
[corrSand corrpSand]=corrcoef(SandVect(:,1),SandVect(:,2));
ClayVect = [ClayAll(:),SMThreshAll(:)];
ClayVect(any(isnan(ClayVect),2),:)=[];
[corrClay corrpClay]=corrcoef(ClayVect(:,1),ClayVect(:,2));

GPMVect = [GPMAll(:),SMThreshAll(:)];
GPMVect(any(isnan(GPMVect),2),:)=[];
[corrGPM corrpGPM]=corrcoef(GPMVect(:,1),GPMVect(:,2));
TreeCoverVect = [TreeCoverMatAll(:),SMThreshAll(:)];
TreeCoverVect(any(isnan(TreeCoverVect),2),:)=[];
[corrTreeCover corrpTreeCover]=corrcoef(TreeCoverVect(:,1),TreeCoverVect(:,2));

MedianSM = nanmedian(SMThreshAll(:));
MedianPotential = nanmedian(MatricMatAll(:));

% MatricMatNoForestAll(MatricMatNoForestAll<-100)=nan;
% [f k] = ksdensity(MatricMatNoForestAll,'bandwidth',0.01)
% plot(k,f)
% hist(MatricMatNoForestAll,10000)
% xlim([-10 0])

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

%% SM Threshold Tropics

SizeSA = size(SMThreshForestSA);
SizeAf = size(SMThreshForestAf);
SizeAu = size(SMThreshForestAu);

figure
set(gcf,'position',[0 0 800 800])
set(gcf,'color','white')

MSA = subplot(1,3,2);
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,SMThreshForestSA)
alpha(atrans)
hold on
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,SMThreshNoForestSA)
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% colorbar
colormap(SMAP_Color_SoilMoisture./300)
set(gca,'color',[0.7 0.7 0.7])
caxis([0 0.6])
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',16)


% 1) Set size of first tile
AMSA = get(MSA,'Position'); % top right fig
set(MSA,'Position',[AMSA(1)-0.06 AMSA(2)+0.51 AMSA(3) AMSA(4)]);
AMSA = get(MSA,'Position'); % top right fig

% 2) Set aspect ratio
% set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(4)*(SizeSA(2)/SizeSA(1)) AMSA(4)]);
set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(3) AMSA(3)*(SizeSA(1)/SizeSA(2))]);
% AMSA = get(MSA,'Position'); % top right fig

% MAf = subplot(3,2,[2 4]);
% figure
MAf = subplot(1,3,3);
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,SMThreshForestAf)
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,SMThreshNoForestAf)
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')

% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
c = colorbar;
% set(get(d,'title'),'string','\DeltaVWC/\DeltaSM')
c.Label.String = 'Estimated Soil Moisture Threshold  [m^3 m^{-3}]';
c.Label.Rotation = 270;
pos = get(c,'Position')
c.Label.Position = [pos(1)+3.5 pos(2)+0.18];
colormap(SMAP_Color_SoilMoisture./300)
set(gca,'color',[0.7 0.7 0.7])
caxis([0 0.6])
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',16)

AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1)-0.1 AMAf(2)+0.5 AMSA(3)*(SizeAf(2)/SizeSA(2)) AMAf(4)]);
AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1) AMAf(2) AMAf(3) AMAf(3)*SizeAf(1)/SizeAf(2)]);
AMAf = get(MAf,'Position'); 

% figure
MAu = subplot(1,3,1);
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,SMThreshForestAu)
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,SMThreshNoForestAu)
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% colorbar
colormap(SMAP_Color_SoilMoisture./300)
set(gca,'color',[0.7 0.7 0.7])
caxis([0 0.6])
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',16)


AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1)+0.005 AMAu(2)+0.6 AMSA(3)*(SizeAu(2)/SizeSA(2)) AMAu(4)]);
AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1) AMAu(2) AMAu(3) AMAu(3)*SizeAu(1)/SizeAu(2)]);

annotation('textbox',[0.1 0.85 0.1 0.1],'String','a','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.32 0.85 0.1 0.1],'String','b','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.565 0.85 0.1 0.1],'String','c','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
fig = gcf;
fig.InvertHardcopy = 'off';
cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
% print('FigS3SMThreshold','-djpeg','-r300') % save as high res PNGd  
% print('FigS3SMThreshold_100','-djpeg','-r100') % save as high res PNGd  
print('FigS3','-depsc','-r300') % save as high res PNGd    

% AMSA = get(MSA,'Position'); % top right fig
% set(MSA,'Position',[AMSA(1) AMSA(2)+0.01 AMSA(3) AMSA(4)]);
% 
% AMAf = get(MAf,'Position'); % bottom right fig
% set(MAf,'Position',[AMAf(1) AMAf(2) AMAf(3) AMAf(4)]);
% 
% AMAu = get(MAu,'Position'); % bottom right fig
% set(MAu,'Position',[AMAu(1)+0.2 AMAu(2)+0.07 AMAu(3) AMAu(4)]);


%% Soil Matric Potential Threshold Tropics

SizeSA = size(SMThreshForestSA);
SizeAf = size(SMThreshForestAf);
SizeAu = size(SMThreshForestAu);

figure
set(gcf,'position',[0 0 800 800])
set(gcf,'color','white')

MSA = subplot(1,3,2);
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,log10(-MatricMatForestSA))
alpha(atrans)
hold on
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,log10(-MatricMatNoForestSA))
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% c1 = colorbar;
% pos = get(c1,'Position')
% c1.Direction = 'reverse';
colormap(flipud(parula))
% c1.Label.String = 'Soil Matric Potential  [MPa]';
% c1.Label.Rotation = 270;
% c1.Label.Position = [pos(1)+3 pos(2)+0.4];
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
% ticks = -1:1:2;
% set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',16)

% 1) Set size of first tile
AMSA = get(MSA,'Position'); % top right fig
set(MSA,'Position',[AMSA(1)-0.06 AMSA(2)+0.51 AMSA(3) AMSA(4)]);
AMSA = get(MSA,'Position'); % top right fig

% 2) Set aspect ratio
% set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(4)*(SizeSA(2)/SizeSA(1)) AMSA(4)]);
set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(3) AMSA(3)*(SizeSA(1)/SizeSA(2))]);
% AMSA = get(MSA,'Position'); % top right fig

MAf = subplot(1,3,3);
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,log10(-MatricMatForestAf))
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,log10(-MatricMatNoForestAf))
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
c1 = colorbar;
pos = get(c1,'Position')
c1.Direction = 'reverse';
colormap(flipud(parula))
c1.Label.String = 'Soil Matric Potential  [MPa]';
c1.Label.Rotation = 270;
c1.Label.Position = [pos(1)+3.4 pos(2)+0.4];
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
ticks = -1:1:2;
set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',16)

AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1)-0.1 AMAf(2)+0.5 AMSA(3)*(SizeAf(2)/SizeSA(2)) AMAf(4)]);
AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1) AMAf(2) AMAf(3) AMAf(3)*SizeAf(1)/SizeAf(2)]);
AMAf = get(MAf,'Position'); 

MAu = subplot(1,3,1);
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,log10(-MatricMatForestAu))
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,log10(-MatricMatNoForestAu))
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
% set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% c1 = colorbar;
% pos = get(c1,'Position')
% c1.Direction = 'reverse';
colormap(flipud(parula))
% c1.Label.String = 'Soil Matric Potential  [MPa]';
% c1.Label.Rotation = 270;
% c1.Label.Position = [pos(1)+3 pos(2)+0.4];
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
% ticks = -1:1:2;
% set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% colorbar
% colormap(SMAP_Color_SoilMoisture./300)
% set(gca,'color',[0.7 0.7 0.7])
% caxis([0 0.6])
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',16)

% set(gcf,'position',[0 0 800 800])
% set(gcf,'color','white')

AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1)+0.005 AMAu(2)+0.6 AMSA(3)*(SizeAu(2)/SizeSA(2)) AMAu(4)]);
AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1) AMAu(2) AMAu(3) AMAu(3)*SizeAu(1)/SizeAu(2)]);


annotation('textbox',[0.1 0.85 0.1 0.1],'String','a','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.32 0.85 0.1 0.1],'String','b','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.565 0.85 0.1 0.1],'String','c','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')

fig = gcf;
fig.InvertHardcopy = 'off';
cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
% print('FigS4MatricThreshold','-djpeg','-r300')
% print('FigS4MatricThreshold_100','-djpeg','-r100')
print('FigS4','-depsc','-r300') % save as high res PNGd    

% AMSA = get(MSA,'Position'); % top right fig
% set(MSA,'Position',[AMSA(1) AMSA(2)+0.01 AMSA(3) AMSA(4)]);
% 
% AMAf = get(MAf,'Position'); % bottom right fig
% set(MAf,'Position',[AMAf(1) AMAf(2) AMAf(3) AMAf(4)]);
% 
% AMAu = get(MAu,'Position'); % bottom right fig
% set(MAu,'Position',[AMAu(1)+0.2 AMAu(2)+0.07 AMAu(3) AMAu(4)]);

%% Boxplots
cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
IGBPMatAf = IGBPMat(rowAf,colAf,:);
IGBPMatSA = IGBPMat(rowSA,colSA,:);
IGBPMatAu = IGBPMat(rowAu,colAu,:);
XticksMod = [7,8,9,10,14];

MatricThreshmatIGBP = nan(5,4000);
for j = 1:size(XticksMod,2)
    iB = XticksMod(j);
    IncDiff1Af = MatricMatNoForestAf;
    IncDiff1SA = MatricMatNoForestSA;
    IncDiff1Au = MatricMatNoForestAu;
    
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
    IncDiff1(isnan(IncDiff1))=[];
    MatricThreshmatIGBP(j,1:size(IncDiff1,1)) = IncDiff1';
end

ClayVec = [ClayMatAf(:); ClayMatAu(:); ClayMatSA(:)];
SMVec = [SMThreshNoForestAf(:); SMThreshNoForestAu(:); SMThreshNoForestSA(:)];
ClayVect = [ClayVec, SMVec];
ClayVect(any(isnan(ClayVect),2),:)=[];

Xticks1Mod = IGBP_Names(XticksMod);
SMThreshmatIGBPModcount = MatricThreshmatIGBP;
SMThreshmatIGBPModcount(isfinite(SMThreshmatIGBPModcount))=1;
SMThreshmatIGBPModcount = nansum(SMThreshmatIGBPModcount,2);

yLog = [0.1:0.1:0.9, 1:1:9, 10:10:100];
xlog = log(yLog)/log(10);

figure
set(gcf,'position',[100 0 1000 500])
set(gcf,'color','white')
Xticks = 1:16;
% IGBP
c = subplot(2,2,[2 4]);
bx = boxplot(log10(-MatricThreshmatIGBP'),'Colors','k');
h=findobj(gca,'tag','Outliers');
set(gca, 'YDir','reverse')
set(gca,'YTick',xlog)
ylim([-1 2])
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'m','FaceAlpha',0.2);
end
hold on
% A = 0:0.01:35;
% plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1Mod,'fontsize',16)
xtickangle(50)
ylabel({'Estimated Soil Matric' ; 'Potential Threshold [MPa]'})
set(gca,'fontname','arial')
set(gca,'fontsize',16)
grid on
set(gca,'fontname','arial')
ticks = (-1:1:2);
set(gca,'YTickLabel',(-10.^xlog))

ax = gca;
labels = string(ax.YAxis.TickLabels); % extract
labels([2:9, 11:18,20:27]) = nan; % remove every other one
ax.YAxis.TickLabels = labels; % set 

b = subplot(2,2,3);
A = log10(-MatricMatAll(:));
A(A<-1)=NaN;
A(A>2)=NaN;
[ks f] = ksdensity(A,'bandwidth',0.2);
plot(f,ks,'LineWidth',2)
set(gca, 'XDir','reverse')
set(gca,'xtick',xlog)
ticks = (-1:1:2);
set(gca,'XTickLabel',(-10.^xlog))
xlabel({'Estimated Soil Matric Potential Threshold [MPa]'})
ylabel({'Probability Density'; 'Function'})
set(gca,'fontname','arial')
set(gca,'fontsize',16)
xlim([-1 2])
grid on

ax = gca;
labels = string(ax.XAxis.TickLabels); % extract
labels([2:9, 11:18,20:27]) = nan; % remove every other one
ax.XAxis.TickLabels = labels; % set 

Xticks = [1 2 3 4 5];
Xticks1 = {'0' '0%-15%' '15%-30%' '30%-45%' '45%-60%'};
a = subplot(2,2,1);
edge = [0 0.15 0.30 0.45 0.6];
[nClay,bin] = histc(ClayAll(:),edge);
boxplot(SMThreshAll(:),bin);
h=findobj(gca,'tag','Outliers');
xlim([1.5 5.5])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1,'fontsize',9)
ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
xlabel('Clay Fraction')
set(gca,'fontname','arial')
set(gca,'fontsize',16)
ylim([0 0.4])
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',0.2);
end
% set(b1,'color','r')
% ba1 = get(get(gca,'children'),'children');   % Get the handles of all the objects
% t1 = get(ba1,'tag');   % List the names of all the objects 
% box1 = ba1(7);   % The 7th object is the first box
% set(box1, 'Color', 'k');   % Set the color of the first box to green
grid on
lines = findobj(gcf, 'type', 'line', 'Tag', 'Median');
set(lines, 'Color', 'k');
set(lines,'linewidth',2);

BoxSamplenClay = nan(1,5);
edge = [0 0.15 0.30 0.45 0.6];
for i = 1:4
    TreeCoverMatn = ClayAll;
    TreeCoverMatn(TreeCoverMatn<edge(i))=nan;
    TreeCoverMatn(TreeCoverMatn>edge(i+1))=nan;    
    TreeCoverMatn(isfinite(TreeCoverMatn))=1;
    IncAfMedianWetn = SMThreshAll.*TreeCoverMatn;
    IncAfMedianWetn(isfinite(IncAfMedianWetn))=1;  
    BoxSamplenClay(i) = nansum(IncAfMedianWetn(:));
end

Aa = get(a,'Position'); % top right fig
set(a,'Position',[Aa(1)-0.03 Aa(2)-0.03 Aa(3)+0.03 Aa(4)+0.03]);
Aa = get(a,'Position'); % top right fig
Bb = get(b,'Position'); % bottom right fig
set(b,'Position',[Aa(1) Bb(2) Aa(3) Aa(4)]);
Cc = get(c,'Position'); % bottom right fig
set(c,'Position',[Cc(1) Cc(2)+0.1 Cc(3) Cc(4)-0.1]);

annotation('textbox',[0.02 0.87 0.1 0.1],'String','a','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.5 0.87 0.1 0.1],'String','c','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.02 0.4 0.1 0.1],'String','b','fontsize',16,'Fontname','arial','FontWeight','bold','EdgeColor','none')
% print('Fig3SMThresholdBox','-djpeg','-r300') % save as high res PNGd
print('Fig3','-depsc','-r1200') % save as high res PNGd    

% % Xticks1 = {'20%-40%' '40%-60%' '60%-80%' '80%-100%'};
% b = subplot(2,2,4)
% edge = [0 0.25 0.50 0.75 1];
% % edge = [0.20 0.40 0.60 0.80 1];
% [nSand,bin] = histc(SandMat(:),edge);
% boxplot(SMThreshAf(:),bin)
% h=findobj(gca,'tag','Outliers');
% xlim([1.5 5.5])
% set(gca,'XTick',Xticks)
% set(gca,'XTickLabel',Xticks1,'fontsize',9)
% ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
% xlabel('Sand Fraction')


MedianSM
MedianPotential
CorrSMClay = corrClay(1,2)
CorrSMSand = corrSand(1,2)
CorrSMGPM = corrGPM(1,2)
corrSMTreeCover = corrTreeCover(1,2)
%%
figure
edge = [0 5 10 30 50 100];
[nClay,bin] = histc(TreeCoverMatAll(:),edge);
boxplot(MatricMatAll(:),bin)
h=findobj(gca,'tag','Outliers');
xlim([1.5 5.5])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1,'fontsize',9)
ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
xlabel('Clay Fraction')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
ylim([-10 0])
delete(h)
h = findobj(gca,'Tag','Box');
for j=1:length(h)
   patch(get(h(j),'XData'),get(h(j),'YData'),'g','FaceAlpha',0.2);
end
grid on
%% Soil Matric Potential Threshold Tropics (OLD CONFIGURATION)

SizeSA = size(SMThreshForestSA);
SizeAf = size(SMThreshForestAf);
SizeAu = size(SMThreshForestAu);

figure
MSA = subplot(3,2,[1 3]);
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,log10(-MatricMatForestSA))
alpha(atrans)
hold on
pcolor(TopLeftCornerLonSA,TopLeftCornerLatSA,log10(-MatricMatNoForestSA))
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% c1 = colorbar;
% pos = get(c1,'Position')
% c1.Direction = 'reverse';
colormap(flipud(parula))
% c1.Label.String = 'Soil Matric Potential  [MPa]';
% c1.Label.Rotation = 270;
% c1.Label.Position = [pos(1)+3 pos(2)+0.4];
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
% ticks = -1:1:2;
% set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',12)

% 1) Set size of first tile
AMSA = get(MSA,'Position'); % top right fig
set(MSA,'Position',[AMSA(1) AMSA(2) AMSA(3) AMSA(4)-0.1]);
AMSA = get(MSA,'Position'); % top right fig

% 2) Set aspect ratio
% pbaspect([1 1*(SizeSA(1)/SizeSA(2)) 1])
AMSA = get(MSA,'Position'); % top right fig
set(MSA,'Position',[AMSA(1)-0.1 AMSA(2) AMSA(4)*(SizeSA(2)/SizeSA(1)) AMSA(4)]);
AMSA = get(MSA,'Position'); % top right fig

MAf = subplot(3,2,[2 4]);
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,log10(-MatricMatForestAf))
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAf,TopLeftCornerLatAf,log10(-MatricMatNoForestAf))
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
c1 = colorbar;
pos = get(c1,'Position')
c1.Direction = 'reverse';
colormap(flipud(parula))
c1.Label.String = 'Soil Matric Potential  [MPa]';
c1.Label.Rotation = 270;
c1.Label.Position = [pos(1)+2.5 pos(2)+0.1];
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
ticks = -1:1:2;
set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',12)

AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1)-0.2 AMAf(2) AMSA(3)*(SizeAf(2)/SizeSA(2)) AMAf(4)]);
AMAf = get(MAf,'Position'); 
set(MAf,'Position',[AMAf(1) AMAf(2) AMAf(3) AMAf(3)*SizeAf(1)/SizeAf(2)]);
AMAf = get(MAf,'Position'); 

MAu = subplot(3,2,[5 6]);
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,log10(-MatricMatForestAu))
alpha(atrans)
hold on
pcolor(TopLeftCornerLonAu,TopLeftCornerLatAu,log10(-MatricMatNoForestAu))
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
set(gca,'ytick',{})
set(gca,'xtick',{})
% c1 = colorbar;
% pos = get(c1,'Position')
% c1.Direction = 'reverse';
colormap(flipud(parula))
% c1.Label.String = 'Soil Matric Potential  [MPa]';
% c1.Label.Rotation = 270;
% c1.Label.Position = [pos(1)+3 pos(2)+0.4];
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
% ticks = -1:1:2;
% set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% colorbar
% colormap(SMAP_Color_SoilMoisture./300)
% set(gca,'color',[0.7 0.7 0.7])
% caxis([0 0.6])
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',12)

set(gcf,'position',[0 0 800 800])
set(gcf,'color','white')

AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1) AMAu(2) AMSA(3)*(SizeAu(2)/SizeSA(2)) AMAu(4)]);
AMAu = get(MAu,'Position'); 
set(MAu,'Position',[AMAu(1) AMAu(2) AMAu(3) AMAu(3)*SizeAu(1)/SizeAu(2)]);

annotation('textbox',[0 0.76 0.1 0.1],'String','a','fontsize',12,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.335 0.76 0.1 0.1],'String','b','fontsize',12,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.3 0.3 0.1 0.1],'String','c','fontsize',12,'Fontname','arial','FontWeight','bold','EdgeColor','none')

AMSA = get(MSA,'Position'); % top right fig
set(MSA,'Position',[AMSA(1) AMSA(2)+0.01 AMSA(3) AMSA(4)]);

AMAf = get(MAf,'Position'); % bottom right fig
set(MAf,'Position',[AMAf(1) AMAf(2) AMAf(3) AMAf(4)]);

AMAu = get(MAu,'Position'); % bottom right fig
set(MAu,'Position',[AMAu(1)+0.2 AMAu(2)+0.07 AMAu(3) AMAu(4)]);

%%
% Histogram
d = subplot(2,1,1);
edge = -2:0.1:10;
c = histogram(log10(-MatricMatAll(:)),edge);
xlim([-1 2])
set(gca, 'XDir','reverse')
set(gca,'xtick',[-1:1:2])
ticks = (-1:1:2);
set(gca,'XTickLabel',(-10.^ticks))
xlabel('Soil Matric Potential  [MPa]')
ylabel('Frequency')
set(gca,'fontname','arial')
set(gca,'fontsize',12)

%%
MatricBlockMedian = nan(6,5);
for j = 1:size(cIGBP,2)
    A = IGBPx;
    A(A~=cIGBP(j))=NaN;
    A(A==cIGBP(j))=1; 
    AMatric = MatricMat.*A;
    AMatric = AMatric(:);
    AMatric(isnan(AMatric))=[];
    MatricBlockMedian(j,1) = cIGBP(j);
    MatricBlockMedian(j,2) = prctile(AMatric(:),25);
    MatricBlockMedian(j,3) = prctile(AMatric(:),50);
    MatricBlockMedian(j,4) = prctile(AMatric(:),75);
    MatricBlockMedian(j,5) = size(AMatric,1);
end



MatricmatIGBPMod =  MatricThreshmatIGBP(XticksMod,:);




%%
%
% figure
% d = subplot(2,2,2);
% [f k] = ksdensity(SMThreshAf(:),'bandwidth',0.01);
% plot(k,f,'LineWidth',2)
% grid on
% xlabel('Estimated SM Threshold [m^3 m^{-3}]')
% ylabel('pdf')
% ylim([0 8])
% set(gca,'fontname','times new roman')
% set(gca,'fontsize',12)

% Xticks = [1 2 3 4 5];
% % Xticks1 = {'0' '0%-25%' '25%-50%' '50%-75%' '75%-100%'};
% Xticks1 = {'0' '0%-15%' '15%-30%' '30%-45%' '45%-60%'};
% d = subplot(2,2,2)
% % edge = [0 0.25 0.50 0.75 1];
% edge = [0 0.15 0.30 0.45 0.6];
% [nClay,bin] = histc(ClayMat(:),edge);
% boxplot(SMThreshAf(:),bin)
% h=findobj(gca,'tag','Outliers');
% xlim([1.5 5.5])
% set(gca,'XTick',Xticks)
% set(gca,'XTickLabel',Xticks1,'fontsize',9)
% ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
% xlabel('Clay Fraction')
% set(gca,'fontname','arial')
% set(gca,'fontsize',11)
% ylim([0 0.4])
% delete(h)
% grid on
% 
% b = subplot(2,2,4);
% boxplot(SMThreshmatIGBPMod')
% h=findobj(gca,'tag','Outliers');
% delete(h)
% set(gca,'XTick',Xticks)
% set(gca,'XTickLabel',Xticks1Mod,'fontsize',9)
% xtickangle(30)
% ylim([0 0.3])
% grid on
% ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
% set(gca,'fontname','times new roman')
% set(gca,'fontsize',12)

% Xticks = [1 2 3 4 5];
% Xticks1 = {'0' '0%-25%' '25%-50%' '50%-75%' '75%-100%'};
% % Xticks1 = {'20%-40%' '40%-60%' '60%-80%' '80%-100%'};
% b = subplot(2,2,4)
% edge = [0 0.25 0.50 0.75 1];
% % edge = [0.20 0.40 0.60 0.80 1];
% [nSand,bin] = histc(SandMat(:),edge);
% boxplot(SMThreshAf(:),bin)
% h=findobj(gca,'tag','Outliers');
% xlim([1.5 5.5])
% set(gca,'XTick',Xticks)
% set(gca,'XTickLabel',Xticks1,'fontsize',9)
% ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
% xlabel('Sand Fraction')
% set(gca,'fontname','arial')
% set(gca,'fontsize',11)
% ylim([0 0.4])
% delete(h)
% grid on



MatricThreshmatIGBP = nan(16,2000);
for iB = 1:size(MatricThreshmatIGBP,1)
    IncDiff1 = MatricMat;
    IGBP1 = IGBPMat(:,:,iB);
    IGBP1(IGBP1<0.5)=NaN;
    IGBP1(IGBP1>=0.5)=1;    
    IncDiff1 = IncDiff1.*IGBP1;
    IncDiff1 = IncDiff1(:);
    IncDiff1(isnan(IncDiff1))=[];
    MatricThreshmatIGBP(iB,1:size(IncDiff1,1)) = IncDiff1';
end

MatricBlockMedian = nan(6,5);
for j = 1:size(cIGBP,2)
    A = IGBPx;
    A(A~=cIGBP(j))=NaN;
    A(A==cIGBP(j))=1; 
    AMatric = MatricMat.*A;
    AMatric = AMatric(:);
    AMatric(isnan(AMatric))=[];
    MatricBlockMedian(j,1) = cIGBP(j);
    MatricBlockMedian(j,2) = prctile(AMatric(:),25);
    MatricBlockMedian(j,3) = prctile(AMatric(:),50);
    MatricBlockMedian(j,4) = prctile(AMatric(:),75);
    MatricBlockMedian(j,5) = size(AMatric,1);
end

Xticks = 1:16;
XticksMod = [8,9,7,10,14];
Xticks1Mod = IGBP_Names(XticksMod);
MatricmatIGBPMod =  MatricThreshmatIGBP(XticksMod,:);
SMThreshmatIGBPModcount = MatricmatIGBPMod;
SMThreshmatIGBPModcount(isfinite(SMThreshmatIGBPModcount))=1;
SMThreshmatIGBPModcount = nansum(SMThreshmatIGBPModcount,2);
%% Map Tropics




%%
% IGBP
b = subplot(2,2,4)
boxplot(log10(-MatricmatIGBPMod'))
h=findobj(gca,'tag','Outliers');
set(gca, 'YDir','reverse')
ylim([-1 2])
delete(h)
hold on
% A = 0:0.01:35;
% plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1Mod,'fontsize',12)
xtickangle(30)
ylabel('Soil Matric Potential  [MPa]')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
grid on
set(gca,'fontname','arial')
ticks = (-1:1:2);
set(gca,'YTickLabel',(-10.^ticks))
% Tree Cover
% Xticks = [1 2 3 4 5 6];
% Xticks1 = {'0' '0-5' '5-10' '10-15' '15-20'};
% d = subplot(2,2,2);
% dge = [0 5 10 15 20];
% [n,bin] = histc(TreeCoverMat(:),edge);
% boxplot(MatricMat(:),bin)
% h=findobj(gca,'tag','Outliers');
% xlim([1.5 5.5])
% delete(h)
% % hold on
% grid on
% hold on
% A = 0:0.01:35;
% plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
% xlabel('Tree Cover (%)')
% set(gca,'fontsize',12)
% set(gca,'fontname','arial')
% ylim([-15 0])
% set(gca,'XTick',Xticks)
% set(gca,'XTickLabel',Xticks1)

% [ks f] = ksdensity(MatricMat(:),'bandwidth',0.01);
% plot(log10(-f),ks)

% Histogram
d = subplot(2,2,2);
edge = -2:0.1:10;
c = histogram(log10(-MatricMat(:)),edge);
xlim([-1 2])
set(gca, 'XDir','reverse')
set(gca,'xtick',[-1:1:2])
ticks = (-1:1:2);
set(gca,'XTickLabel',(-10.^ticks))
xlabel('Soil Matric Potential  [MPa]')
ylabel('Frequency')
set(gca,'fontname','arial')
set(gca,'fontsize',12)
% set(c,'tick',ticks,'TickLabels',-10.^ticks)

% edge = -50:0.01:0;
% histogram(MatricMat(:),edge)
% xlim([-15 0])
% semilogx(log10(-f),ks)

a = subplot(2,2,[1 3]);
pcolor(TopLeftCornerLon,TopLeftCornerLat,log10(-MatricMatForest))
alpha(atrans)
hold on
pcolor(TopLeftCornerLon,TopLeftCornerLat,log10(-MatricMatNoForest))
% title('Estimated SM Threshold [m^3 m^{-3}]')
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.5,'Color','k')
set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
c1 = colorbar;
pos = get(c1,'Position')
c1.Direction = 'reverse';
colormap(flipud(parula))
c1.Label.String = 'Soil Matric Potential  [MPa]';
c1.Label.Rotation = 270;
c1.Label.Position = [pos(1)+3 pos(2)+0.4];
% colormap(SMAP_Color_SoilMoisture./300)
set(gca,'color',[0.7 0.7 0.7])
caxis([-1,2])
ticks = -1:1:2;
set(c1,'Ticks',ticks,'TickLabels',-10.^ticks)
% annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','arial')
set(gca,'fontsize',14)

% annotation('textbox',[0.4 0.23 0.1 0.1],'String','A','fontsize',20,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.89 0.62 0.1 0.1],'String','B','fontsize',20,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.89 0.25 0.1 0.1],'String','C','fontsize',20,'Fontname','Times New Roman','EdgeColor','none')

annotation('textbox',[0 0.85 0.1 0.1],'String','a','fontsize',12,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.53 0.85 0.1 0.1],'String','b','fontsize',12,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.53 0.47 0.1 0.1],'String','c','fontsize',12,'Fontname','arial','FontWeight','bold','EdgeColor','none')

set(gcf,'color','white')
% A = get(a,'Position'); % Left fig
% set(a,'Position',[A(1)-0.1 A(2)+0.15 A(3)+0.07 A(4)-0.15]);
% B = get(b,'Position'); % bottom right fig
% set(b,'Position',[B(1)+0.014 B(2)+0.16 B(3)+0.007 B(4)-0.05]);
% D = get(d,'Position'); % top right fig
% set(d,'Position',[B(1)+0.014 D(2)+0.05 B(3)+0.007 B(4)-0.05]);
set(gcf,'position',[100 100 1100 600])

A = get(a,'Position'); % Left fig
set(a,'Position',[A(1)-0.1 A(2)+0.15 A(3)+0.12 A(4)-0.15]);
D = get(d,'Position'); % top right fig
set(d,'Position',[D(1)+0.04 D(2)+0.05 D(3) D(4)-0.05]);
% B = get(b,'Position'); % bottom right fig
% set(b,'Position',[B(1)+0.015 B(2)+0.02 B(3)+0.007 B(4)+0.06]);
B = get(b,'Position'); % bottom right fig
set(b,'Position',[D(1)+0.04 A(2)+0.15 D(3) D(4)-0.05]);


%%
figure
subplot(1,2,2)
scatter(SandMat(:),SMThresh(:))
xlabel('Sand')
set(gca,'fontname','times new roman')
set(gca,'fontsize',12)
% ylabel('SM Threshold')
subplot(1,2,1)
scatter(ClayMat(:),SMThresh(:))
xlabel('Clay')
ylabel('SM Threshold')
set(gca,'fontname','times new roman')
set(gca,'fontsize',12)

set(gcf,'position',[100 100 1000 400])
set(gcf,'color','white')
%%
Nh = 45  ;        
Nv = 45;
Y1 = SMThresh(:) ;                                        
X1 = SandMat(:);                                      
h1 = hist3([X1 Y1],[Nh Nv]);                                            
xb = linspace(0,max(X1),size(h1,1));               
yb = linspace(0,max(Y1),size(h1,2));               

Nh = 45  ;        
Nv = 45;
Y1 = SMThresh(:) ;                                        
X2 = ClayMat(:);                                      
h2 = hist3([X2 Y1],[Nh Nv]);                                            
xb2 = linspace(0,max(X2),size(h1,1));               
yb2 = linspace(0,max(Y1),size(h1,2));

figure
subplot(1,2,2)
pcolor(xb,yb,h1')  
shading('flat')
xlim([0.2 1])
ylim([0 0.3])
colormap(flipud(gray))
% colorbar
caxis([0 30])
xlabel('Sand Fraction (%)')
ylabel('SM Threshold')
% title('Africa Depolarization')
% set(gcf,'color','white')
set(gca,'fontsize',14)
set(gca,'Fontname','Times New Roman')

subplot(1,2,1)
pcolor(xb2,yb2,h2')  
shading('flat')
xlim([0 0.6])
ylim([0 0.3])
colormap(flipud(gray))
% colorbar
caxis([0 30])
xlabel('Clay Fraction (%)')
ylabel('SM Threshold')
% title('Africa Depolarization')
% set(gcf,'color','white')
set(gca,'fontsize',14)
set(gca,'Fontname','Times New Roman')
set(gcf,'position',[100 100 1000 400])
set(gcf,'color','white')
%%

% Xticks = [1 2 3 4 5];
% Xticks1 = {'0' '0-5' '5-15' '15-30' '30-60'};
% d = subplot(2,2,2)
% edge = [0 0.25 0.50 0.75 1];
% [n,bin] = histc(SandMat(:),edge);
% boxplot(SMThreshAf(:),bin)
% h=findobj(gca,'tag','Outliers');
% xlim([1.5 5.5])
% delete(h)


%%
[~,~,~,~,...
   TopLeftCornerLat,~,~,~] = InterpEdgesCorners(alat);
TopLeftCornerLat(:,1) = TopLeftCornerLat(:,2);
TopLeftCornerLat(1,:) = TopLeftCornerLat(2,:)+0.5;

[~,~,~,~,...
   TopLeftCornerLon,~,~,~] = InterpEdgesCorners(alon);
TopLeftCornerLon(1,:) = TopLeftCornerLon(2,:);
TopLeftCornerLon(:,1) = TopLeftCornerLon(:,2)-0.5;



MatricThreshmatIGBP = nan(16,2000);
for iB = 1:size(MatricThreshmatIGBP,1)
    IncDiff1 = SMThresh;
    IGBP1 = IGBPMat(:,:,iB);
    IGBP1(IGBP1<0.5)=NaN;
    IGBP1(IGBP1>=0.5)=1;    
    IncDiff1 = IncDiff1.*IGBP1;
    IncDiff1 = IncDiff1(:);
    IncDiff1(isnan(IncDiff1))=[];
    MatricThreshmatIGBP(iB,1:size(IncDiff1,1)) = IncDiff1';
end

Xticks = 1:16;
XticksMod = [7,10,14,9,8];
Xticks1Mod = IGBP_Names(XticksMod);
MatricmatIGBPMod =  MatricThreshmatIGBP(XticksMod,:);
SMThreshmatIGBPModcount = MatricmatIGBPMod;
SMThreshmatIGBPModcount(isfinite(SMThreshmatIGBPModcount))=1;
SMThreshmatIGBPModcount = nansum(SMThreshmatIGBPModcount,2);
%
figure
d = subplot(2,2,2);
[f k] = ksdensity(SMThresh(:),'bandwidth',0.01);
plot(k,f,'LineWidth',2)
grid on
xlabel('Estimated SM Threshold [m^3 m^{-3}]')
ylabel('pdf')
ylim([0 8])
set(gca,'fontname','times new roman')
set(gca,'fontsize',12)

b = subplot(2,2,4);
boxplot(MatricmatIGBPMod')
h=findobj(gca,'tag','Outliers');
delete(h)
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1Mod,'fontsize',9)
xtickangle(30)
ylim([0 0.3])
grid on
ylabel({'Estimated SM';'Threshold [m^3 m^{-3}]'})
set(gca,'fontname','times new roman')
set(gca,'fontsize',12)

a = subplot(2,2,[1 3]);
pcolor(TopLeftCornerLon,TopLeftCornerLat,SMThreshAfForest)
alpha(atrans)
hold on
pcolor(TopLeftCornerLon,TopLeftCornerLat,SMThreshAfNoForest)
% title('Estimated SM Threshold [m^3 m^{-3}]')
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.01,'Color','k')
set(gcf,'position',[100 100 1000 500])
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
colorbar
colormap(SMAP_Color_SoilMoisture./300)
set(gca,'color',[0.7 0.7 0.7])
caxis([0 0.6])
annotation('textbox',[0.7 0.83 0.1 0.1],'String','Median = 0.14 m^3/m^3','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
% annotation('textbox',[0.7 0.78 0.1 0.1],'String','Forests Not Included','fontsize',12,'Fontname','Times New Roman','EdgeColor','none')
set(gca,'fontname','times new roman')
set(gca,'fontsize',13)

annotation('textbox',[0.4 0.23 0.1 0.1],'String','A','fontsize',20,'Fontname','Times New Roman','EdgeColor','none')
annotation('textbox',[0.87 0.65 0.1 0.1],'String','B','fontsize',20,'Fontname','Times New Roman','EdgeColor','none')
annotation('textbox',[0.87 0.23 0.1 0.1],'String','C','fontsize',20,'Fontname','Times New Roman','EdgeColor','none')

set(gcf,'color','white')
A = get(a,'Position'); % Left fig
set(a,'Position',[A(1)-0.1 A(2)+0.15 A(3)+0.07 A(4)-0.15]);
B = get(b,'Position'); % bottom right fig
set(b,'Position',[B(1)-0.005 B(2)-0.055 B(3)+0.007 B(4)+0.15]);
D = get(d,'Position'); % top right fig
set(d,'Position',[D(1)+0.014 D(2)+0.05 D(3)-0.014 D(4)-0.05]);
set(gcf,'position',[100 100 1000 600])
%%


% IncAfMedianWetForest = IncAfMedianWet.*ForestOnlyMat;
% IncAfMedianDryForest = IncAfMedianDry.*ForestOnlyMat;
% IncAfMedianWet = IncAfMedianWet.*ForestFlagMat;
% IncAfMedianDry = IncAfMedianDry.*ForestFlagMat;

% IncAfMedianWetPerc = IncAfMedianWet;
% IncAfMedianWetPerc(LIDARAf>10)=NaN;
% IncAfMedianWetPerc(LIDARAf<5)=NaN;
% IncAfMedianWetPerc = IncAfMedianWetPerc(:);
% IncAfMedianWetPerc(isnan(IncAfMedianWetPerc))=[];
% [k f] = ksdensity(IncAfMedianWetPerc,'bandwidth',0.01);
% plot(f,k)


Xticks = [1 2 3 4 5];
Xticks1 = {'0' '2.5' '7.5' '12.5' '17.5'};
dim1 = [.135 .64 .1 .1];
str1 = ['SM<0.12 m^3/m^3'];

% LIDARAf(isnan(LIDARAf))=1;

a = subplot(1,2,2)
set(0, 'DefaultAxesFontSize',14)
set(0,'defaultlinelinewidth',1)
set(0, 'DefaultAxesFontName','Arial')
Ncr = 100;
Ncb = 300;
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBu = [red;flipud(blu)];
pcolor(TopLeftCornerLon,TopLeftCornerLat,IncAfMedianWetForest)
alpha(atrans)
hold on
pcolor(TopLeftCornerLon,TopLeftCornerLat,IncAfMedianWet)
% alpha(1) 

c = colorbar;
% set(get(d,'title'),'string','\DeltaVWC/\DeltaSM')
c.Label.String = '\DeltaVWC/\DeltaSM  [(kg/m^2)/(m^3/m^3)]';
c.Label.Rotation = 270;
pos = get(c,'Position')
c.Label.Position = [pos(1)+2.2 pos(2)+3.8];
shading flat
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
hold on
geoshow(clat,clon,'LineWidth',0.01,'Color','k')
% axis([-18 52 -35 18])
axis([-20 52 -35 18])
% caxis([-5 5])
caxis([-5 15])
colormap(cmapRedBu)
set(gca,'color',[0.7 0.7 0.7])
% xlabel('Longitude')
% ylabel('Latitude')
% title(['Median(\Delta\tau/\Delta\theta) for \theta>' num2str( ThreshInc)])
% nam = 'IncAfr.fig'                                               
set(gcf,'PaperPositionMode','auto')
% savefig(nam)
set(gca,'fontsize',14)
set(gca,'fontname','times new roman')
axes('Position',[.585 .3 .095 .3])
% IncAfWet(IncAfWet==-1000)=NaN;
% IncAfWet(IncAfWet==0)=NaN;
% edge = 5:10:35;
% edge = 0:5:20;
edge = 0:5:20;
[n,bin] = histc(LIDARAf(:),edge);
% figure
% h = hist3([LIDARAf(:), IncAf(:)],[bin1, bin1]);
% xVec = linspace(min(LIDARAf(:)),max(LIDARAf(:)),bin1);
% yVec = linspace(min(IncAf(:)),max(IncAf(:)),bin1);
% pcolor(xVec,yVec,h); shading flat
% scatter(LIDARAf(:),IncAf(:))
boxplot(IncAfMedianWet(:),bin)
h=findobj(gca,'tag','Outliers');
xlim([1.5 5.5])
delete(h)
% hold on
grid on
hold on
A = 0:0.01:35;
plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
% A = 0:50;
% plot(A,zeros(1,size(A,2)),'--k')
% colormap(flipud(gray))
% colorbar
xlabel('LIDAR Height (m)')
% ylabel('E[{\delta \tau}/{\delta \theta}|\theta]')
set(gca,'fontsize',12)
set(gca,'fontname','times new roman')
ylim([-5 5])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
annotation('textbox',dim1,'String',str1,'fontsize',14,'Fontname','Times New Roman','EdgeColor','none')
annotation('textbox',[0.83 0.18 0.1 0.1],'String','B','fontsize',16,'Fontname','Times New Roman','EdgeColor','none')

Xticks = [1 2 3 4 5];
Xticks1 = {'0' '2.5' '7.5' '12.5' '17.5'};
dim1 = [.525 .64 .1 .1];
str1 = ['SM>0.12 m^3/m^3'];

b = subplot(1,2,1)
set(0, 'DefaultAxesFontSize',14)
set(0,'defaultlinelinewidth',1)
set(0, 'DefaultAxesFontName','Arial')
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBu = [red;flipud(blu)];
% Ncb = 100 ;                                                       
% red = [ ones(Ncb,1)    [1:Ncb]'/Ncb  [1:Ncb]'/Ncb ] ;             
% blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
% cmapRedBu = [red;flipud(blu)];
pcolor(TopLeftCornerLon,TopLeftCornerLat,IncAfMedianDryForest)
% set(v,'EdgeColor','k')
alpha(atrans)
hold on
pcolor(TopLeftCornerLon,TopLeftCornerLat,IncAfMedianDry)
% pcolor(TopLeftCornerLon,TopLeftCornerLat,IncAfMedianDry)
% colorbar
set(gca,'yticklabel',{})
set(gca,'xticklabel',{})
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.01,'Color','k')
% axis([-18 52 -35 18])
axis([-20 52 -35 18])
% caxis([-5 5])
caxis([-5 15])
colormap(cmapRedBu)
set(gca,'color',[0.7 0.7 0.7])
% xlabel('Longitude')
% ylabel('Latitude')
% title(['Median(\Delta\tau/\Delta\theta) for \theta>' num2str( ThreshInc)])
% nam = 'IncAfr.fig'                                               
set(gcf,'PaperPositionMode','auto')
% savefig(nam)
set(gca,'fontsize',14.1)
set(gca,'fontname','times new roman')

axes('Position',[.19 .3 .095 .3])
% IncAfDry(IncAfDry==-1000)=NaN;
% IncAfDry(IncAfDry==0)=NaN;
edge = 0:5:20;
[n,bin] = histc(LIDARAf(:),edge);
% figure
% h = hist3([LIDARAf(:), IncAf(:)],[bin1, bin1]);
% xVec = linspace(min(LIDARAf(:)),max(LIDARAf(:)),bin1);
% yVec = linspace(min(IncAf(:)),max(IncAf(:)),bin1);
% pcolor(xVec,yVec,h); shading flat
% scatter(LIDARAf(:),IncAf(:))

% The two commented out commands below show that the left most boxplot is
% plotting all lidar areas that are nan's or don't fit into the predefined
% boxes
% IncAfMedianDry(isfinite(LIDARAf))=NaN;
% boxplot(IncAfMedianDry(:))
boxplot(IncAfMedianDry(:),bin)
grid on
h=findobj(gca,'tag','Outliers');
delete(h)
hold on
A = 0:0.01:35;
plot(A,zeros(size(A,2),1),'--k','LineWidth',1)
xlim([1.5 5.5])
% hold on

% A = 0:50;
% plot(A,zeros(1,size(A,2)),'--k')
% colormap(flipud(gray))
% colorbar
xlabel('LIDAR Height (m)')
% ylabel('E[{\delta \tau}/{\delta \theta}|\theta]')
set(gca,'fontsize',12)
set(gca,'fontname','times new roman')
ylim([-10 20])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
annotation('textbox',dim1,'String',str1,'fontsize',14,'Fontname','Times New Roman','EdgeColor','none')
annotation('textbox',[0.435 0.18 0.1 0.1],'String','A','fontsize',16,'Fontname','Times New Roman','EdgeColor','none')

set(gcf,'position',[100 100 1000 400])
set(gcf,'color','white')
A = get(a,'Position');
set(a,'Position',[A(1)-0.05 A(2)+0.05 A(3)+0.15 A(4)]);
B = get(b,'Position');
set(b,'Position',[B(1)-0.01 B(2)+0.05 B(3)+0.04 B(4)]);