clc
clear
  SmoothDay = 91; %odd number
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\9km Half Year')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\9km Half Year')
cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP (1)\MTDCA_fmincon')

SMmat1 = matfile('MTDCAmincon_SM_201504_201509_9km');
SMmat2 = matfile('MTDCAmincon_SM_201510_201603_9km');
SMmat3 = matfile('MTDCAmincon_SM_201604_201609_9km');
SMmat4 = matfile('MTDCAmincon_SM_201610_201703_9km');
SMmat5 = matfile('MTDCAmincon_SM_201704_201709_9km');
SMmat6 = matfile('MTDCAmincon_SM_201710_201712_9km');
TAUmat1 = matfile('MTDCAmincon_TAU_201504_201509_9km');
TAUmat2 = matfile('MTDCAmincon_TAU_201510_201603_9km');
TAUmat3 = matfile('MTDCAmincon_TAU_201604_201609_9km');
TAUmat4 = matfile('MTDCAmincon_TAU_201610_201703_9km');
TAUmat5 = matfile('MTDCAmincon_TAU_201704_201709_9km');
TAUmat6 = matfile('MTDCAmincon_TAU_201710_201712_9km');

cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\General Data')
% load('MTDCA2_v13','IGBP_36km','water_f','lat','lon','IGBP_Names')
load('Ancillary_9km_3dB','water_9km_3dB','zlidar_9km_3dB')
load('SMAPCenterCoordinates9KM')

cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
load('SMThresholdEstimate_NullRemove_Globe')

lidar = zlidar_9km_3dB;

NPeriod = 1001;
waterf = water_9km_3dB;
waterf(waterf>0.05)=NaN;
waterf(isfinite(waterf))=1;

lon = SMAPCenterLongitudes;
lat = SMAPCenterLatitudes;  
latAve = nanmean(lat,2);
lonAve = nanmean(lon,1);
% Minimum Number of Overpasses in Drydown
dLength = 4;
 ThreshInc = 0.12; 
dll = 0.5;   %aggregate to dll degrees                                        
alat = [ -60 : dll : 60 ];                              
alon = [ -179 : dll : 179 ];                            
nlat = length(alat);                                    
nlon = length(alon); 
alat = flipud(transpose(repmat(alat,nlon,1))); %lat matrix         
alon = repmat(alon,nlat,1);                  %lon matrix     

IncMedian = NaN(nlat,nlon);   %Median increment matrix   
Inc25 = NaN(nlat,nlon);   %Median increment matrix   
Inc75 = NaN(nlat,nlon);   %Median increment matrix   
LIDAR = NaN(nlat,nlon);   %Lidar Africa                         
DryCount = NaN(nlat,nlon); 
IncCount = NaN(nlat,nlon); 
% Do only if SMThreshold is finite

for ilat = 1 : nlat     
    sprintf(['irow = ' sprintf('%0.2f',(ilat/nlat)) ])
         for ilon = 1 : nlon   

% North Africa
% for ilat = 48 : 62     
%     sprintf(['irow = ' sprintf('%0.2f',(ilat/nlat)) ])
%          for ilon = 358 : 429      

% Find SMAP EASE2 Rows/Cols
    [DegRow] = find(latAve>= alat(ilat,ilon)-(dll/2)   & ...
             latAve<= alat(ilat,ilon)+(dll/2));
    [DegCol] = find(lonAve>= alon(ilat,ilon)-(dll/2)   & ...
             lonAve<= alon(ilat,ilon)+(dll/2)   ); 
 
ThreshInc = SMThresh(ilat,ilon);
         
 % Load only SMAP pixels in this subset        
%  cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\9km Half Year')
 cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP (1)\MTDCA_fmincon')
AMat1 = SMmat1.MTDCA_SM_201504_201509(DegRow,DegCol,:);
AMat2 = SMmat2.MTDCA_SM_201510_201603(DegRow,DegCol,:);
AMat3 = SMmat3.MTDCA_SM_201604_201609(DegRow,DegCol,:);
AMat4 = SMmat4.MTDCA_SM_201610_201703(DegRow,DegCol,:);
AMat5 = SMmat5.MTDCA_SM_201704_201709(DegRow,DegCol,:);
AMat6 = SMmat6.MTDCA_SM_201710_201712(DegRow,DegCol,:);
mv_DCA2 = squeeze(cat(3,AMat1,AMat2,AMat3,AMat4,AMat5,AMat6));
AMat1 = TAUmat1.MTDCA_TAU_201504_201509(DegRow,DegCol,:);
AMat2 = TAUmat2.MTDCA_TAU_201510_201603(DegRow,DegCol,:);
AMat3 = TAUmat3.MTDCA_TAU_201604_201609(DegRow,DegCol,:);
AMat4 = TAUmat4.MTDCA_TAU_201610_201703(DegRow,DegCol,:);
AMat5 = TAUmat5.MTDCA_TAU_201704_201709(DegRow,DegCol,:);
AMat6 = TAUmat6.MTDCA_TAU_201710_201712(DegRow,DegCol,:);
tau_DCA2 = squeeze(cat(3,AMat1,AMat2,AMat3,AMat4,AMat5,AMat6));
tau_DCA2 = tau_DCA2./0.11;

NPeriodSmooth = 366;
tau_DCA2_Smooth = nan(size(tau_DCA2,1),size(tau_DCA2,2),NPeriodSmooth,3);
tau_DCA2_Smooth(:,:,:,1) = tau_DCA2(:,:,1:366);
tau_DCA2_Smooth(:,:,:,2) = cat(3,tau_DCA2(:,:,367:700),NaN(size(tau_DCA2,1)...
                            ,size(tau_DCA2,2)),tau_DCA2(:,:,701:731));
tau_DCA2_Smooth(:,:,:,3) = cat(3,tau_DCA2(:,:,732:end),NaN(size(tau_DCA2,1)...
                            ,size(tau_DCA2,2),366-length(tau_DCA2(1,1,732:end))));    
tau_DCA2_Smooth = nanmean(tau_DCA2_Smooth,4);
MTDCAtauSmooth1 = repmat(tau_DCA2_Smooth,[1 1 3]);
MTDCAtauSmooth2 = NaN(size(tau_DCA2,1),size(tau_DCA2,2),NPeriodSmooth);
%     k=1;
Smooth1 = (SmoothDay-1)/2;
for iDay = 1:NPeriodSmooth
    MTDCAtauSmooth2(:,:,iDay) = nanmean(MTDCAtauSmooth1(:,:,NPeriodSmooth+iDay...
        -Smooth1:NPeriodSmooth+iDay+Smooth1),3);
%         k = k+1;
end
tau_DCA2_Smooth = repmat(MTDCAtauSmooth2,[1 1 3]);
tau_DCA2_Smooth = tau_DCA2_Smooth(:,:,1:NPeriod);

mv_DCA2_Smooth = nan(size(mv_DCA2,1),size(mv_DCA2,2),NPeriodSmooth,3);
mv_DCA2_Smooth(:,:,:,1) = mv_DCA2(:,:,1:366);
mv_DCA2_Smooth(:,:,:,2) = cat(3,mv_DCA2(:,:,367:700),NaN(size(mv_DCA2,1)...
                            ,size(mv_DCA2,2)),mv_DCA2(:,:,701:731));
mv_DCA2_Smooth(:,:,:,3) = cat(3,mv_DCA2(:,:,732:end),NaN(size(mv_DCA2,1)...
                            ,size(mv_DCA2,2),366-length(mv_DCA2(1,1,732:end))));    
mv_DCA2_Smooth = nanmean(mv_DCA2_Smooth,4);
MTDCAmvSmooth1 = repmat(mv_DCA2_Smooth,[1 1 3]);
MTDCAmvSmooth2 = NaN(size(mv_DCA2,1),size(mv_DCA2,2),NPeriodSmooth);
%     k=1;
Smooth1 = (SmoothDay-1)/2;
for iDay = 1:NPeriodSmooth
    MTDCAmvSmooth2(:,:,iDay) = nanmean(MTDCAmvSmooth1(:,:,NPeriodSmooth+iDay...
        -Smooth1:NPeriodSmooth+iDay+Smooth1),3);
%         k = k+1;
end
mv_DCA2_Smooth = repmat(MTDCAmvSmooth2,[1 1 3]);
mv_DCA2_Smooth = mv_DCA2_Smooth(:,:,1:NPeriod);

 ji = 0; 
 Np = 40000; 
VecInc = NaN(round(NPeriod/3),Np);
% Do only if SMThreshold is finite
if isfinite(ThreshInc)
% Loop over each pixel within X degree upscaled plot
    for SMAProw = 1:size(DegRow,1)  
         for SMAPcol = 1:size(DegCol,2)                  
         smt  = squeeze(mv_DCA2(SMAProw,SMAPcol,1:NPeriod));                  
         taut = squeeze(tau_DCA2(SMAProw,SMAPcol,1:NPeriod));
         tauts = squeeze(tau_DCA2_Smooth(SMAProw,SMAPcol,1:NPeriod)); 
         smts = squeeze(mv_DCA2_Smooth(SMAProw,SMAPcol,1:NPeriod));
        
     smt(end) = [] ; %remove ends for next part of analysis                                                  
    taut(end) = [];  
    tauts(end) = [];  
    smts(end) = [] ;
    
          taut(isnan(smt)) = [] ;
           tauts(isnan(smt)) = [] ;  
           smts(isnan(smt)) = [];             
           smt(isnan(smt)) = [];                                          
            ipb = find( smt(2:length(smt)-1) > smt(1:length(smt)-2) & ...
                        smt(2:length(smt)-1) > smt(3:length(smt)  ));  
            ipb = ipb + 1;                                                      
            ipe = find( smt(2:length(smt)-1) < smt(1:length(smt)-2) & ...
                        smt(2:length(smt)-1) < smt(3:length(smt)  ));  
            ipe = ipe + 1;                                                      
           for k = 1 : length(ipb) %loop through each dry down
                ke = min(find(ipe>ipb(k))); %find next trough after peak                           
               smv =  smt(ipb(k):ipe(ke)) ; %SM between trough and peak 
               smvs =  smts(ipb(k):ipe(ke)) ; %SM between trough and peak                                         
              tauv = taut(ipb(k):ipe(ke)); %corresponding tau during this time      
              tauvs = tauts(ipb(k):ipe(ke)); %corresponding tau during this time                                                      
             tauv(isnan(smv)) = [];  
             tauvs(isnan(smv)) = []; 
              smvs(isnan(smv)) = [];                                                 
              smv(isnan(smv)) = [];                                    
                                    
             if (length(smv)>=dLength) 
            tauv(smv>ThreshInc) = [] ;   %remove dry soil activity 
          tauvs(smv>ThreshInc) = [] ;  
           smvs(smv>ThreshInc) = [];             
            smv(smv>ThreshInc) = [];
            
            if size(smv,1)<2|size(tauv,1)<2
           else
           SlopeNulltau = diff(tauvs);
            SlopeRegtau = diff(tauv);
            SlopeDifftau = SlopeRegtau-SlopeNulltau;
            tauvC = nan(size(tauv));
            tauvC(1) = tauv(1);
           SlopeNullsm = diff(smvs);
            SlopeRegsm = diff(smv);
            SlopeDiffsm = SlopeRegsm-SlopeNullsm;
            smvC = nan(size(tauv));
            smvC(1) = smv(1);              
            for a = 2:size(SlopeDifftau,1)+1
                tauvC(a) = tauvC(a-1)+SlopeDifftau(a-1);
                smvC(a) = smvC(a-1)+SlopeDiffsm(a-1);                                
            end  
            
                ji = ji + 1;                                           
                  if (ji>Np)
                      disp(['Error: Increase Vector Size for Africa Map'])
                  end
                  
              VecInc(1:length(smv)-1,ji) = (tauvC(2:end)-tauvC(1:end-1))./ ...  
                                      (smvC(2:end)-smvC(1:end-1));  
         
              end
           end
           end % End dry down loop within pixel
        end %pixel column
    end %pixel row
else
end % End SM threshold isfinite check
     VecInc = VecInc(:);
     VecInc(isnan(VecInc))=[];
     VecIncFin = VecInc;
      VecIncFin(VecIncFin<0|VecIncFin>=0)=1;
      lidarBlock = lidar(DegRow,DegCol);
      LIDAR(ilat,ilon) = nanmean(lidarBlock(:));
      DryCount(ilat,ilon) = ji;
      IncCount(ilat,ilon) = nansum(VecIncFin(:));  
      IncMedian(ilat,ilon) = median(VecInc(:));
      Inc25(ilat,ilon) = prctile(VecInc(:),25);
      Inc75(ilat,ilon) = prctile(VecInc(:),75);      
%       if isnan(nansum(VecIncFin(:)))
%       else          
%            
%       end
%       A = VecInc(:);
%       A(isnan(A))=[];
%       A = sort(A);
%       [ki f] = ksdensity(A(:),'bandwidth',0.01);
%       plot(f,ki)
%       xlim([-10 10])
end
end

cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
save -v7.3 IncrementMedian_HalfDeg_4Inc_DryEnd_NullRemove_SMThreshEstNoNull.mat ...
    LIDAR DryCount IncCount IncMedian Inc25 Inc75 alat alon ThreshInc
%%
clc
clear
load coast   
clat = lat;                                                        
clon = long ; 
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\MIT\Manuscripts\WaterExchangeNatureGeoscience')
cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
load('ForestFlag')
load('SMAPCenterCoordinates9KM')
cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\MIT\Manuscripts\WaterExchangeNatureGeoscience')
load('IncrementMedian_HalfDeg_4Inc_DryEnd')
dll = 0.5; 
AfricaRow = 300:1300;
AfricaCol = 1700:2500;
ForestFlag = ForestFlag(AfricaRow,AfricaCol);
lon = SMAPCenterLongitudes(AfricaRow,AfricaCol);
lat = SMAPCenterLatitudes(AfricaRow,AfricaCol);  
latAve = nanmean(lat,2);
lonAve = nanmean(lon,1);
alat = [ -35 : dll : 18 ];                              
alon = [ -18 : dll : 52 ];                             
nlat = length(alat);                                    
nlon = length(alon); 
alat = flipud(transpose(repmat(alat,nlon,1))); %lat matrix         
alon = repmat(alon,nlat,1);                  %lon matrix     
ForestFlagMat = NaN(nlat,nlon);
for ilat = 1 : nlat     
    sprintf(['irow = ' sprintf('%0.2f',(ilat/nlat)) ])
         for ilon = 1 : nlon   
% Find SMAP EASE2 Rows/Cols
    [DegRow] = find(latAve>= alat(ilat,ilon)-(dll/2)   & ...
             latAve<= alat(ilat,ilon)+(dll/2));
    [DegCol] = find(lonAve>= alon(ilat,ilon)-(dll/2)   & ...
             lonAve<= alon(ilat,ilon)+(dll/2)   ); 
           ForestBlock = ForestFlag(DegRow,DegCol);
           if nansum(ForestBlock(:))==0
           ForestFlagMat(ilat,ilon) = 1;
           else
           ForestFlagMat(ilat,ilon) = nan;           
           end
         end
end


DryCount(IncCount<50)=NaN;
IncAf(DryCount<50)=NaN;
IncCount(IncCount<50)=NaN;
% IncAf(isnan(IncAf))=-1000;
  IncAf = IncAf.*ForestFlagMat;
  LIDAR = LIDAR.*ForestFlagMat;
[~,~,~,~,...
   TopLeftCornerLat,~,~,~] = InterpEdgesCorners(alat);
TopLeftCornerLat(:,1) = TopLeftCornerLat(:,2);
TopLeftCornerLat(1,:) = TopLeftCornerLat(2,:)+0.5;

[~,~,~,~,...
   TopLeftCornerLon,~,~,~] = InterpEdgesCorners(alon);
TopLeftCornerLon(1,:) = TopLeftCornerLon(2,:);
TopLeftCornerLon(:,1) = TopLeftCornerLon(:,2)-0.5;

Xticks = [1 2 3 4];
Xticks1 = {'5' '10' '15' '20'};
dim1 = [.15 .54 .1 .1];
str1 = ['SM<0.12 m^3/m^3'];

figure
set(0, 'DefaultAxesFontSize',14)
set(0,'defaultlinelinewidth',1)
set(0, 'DefaultAxesFontName','Arial')
Ncb = 100 ;                                                       
red = [ ones(Ncb,1)    [1:Ncb]'/Ncb  [1:Ncb]'/Ncb ] ;             
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];     
cmapRedBu = [red;flipud(blu)];
pcolor(TopLeftCornerLon,TopLeftCornerLat,IncAf)
colorbar
shading flat
hold on
geoshow(clat,clon,'LineWidth',0.01,'Color','k')
axis([-18 52 -35 18])
caxis([-2.5 2.5])
colormap(cmapRedBu)
set(gca,'color',[0.7 0.7 0.7])
xlabel('Longitude')
ylabel('Latitude')
% title(['Median(\Delta\tau/\Delta\theta) for \theta>' num2str( ThreshInc)])
% nam = 'IncAfr.fig'                                               
set(gcf,'PaperPositionMode','auto')
% savefig(nam)
set(gca,'fontsize',16)
set(gca,'fontname','times new roman')

axes('Position',[.16 .25 .23 .25])
% IncAf(IncAf==-1000)=NaN;
% IncAf(IncAf==0)=NaN;
edge = 5:5:20;
[n,bin] = histc(LIDAR(:),edge);
% colormap(jet)

% figure
% h = hist3([LIDARAf(:), IncAf(:)],[bin1, bin1]);
% xVec = linspace(min(LIDARAf(:)),max(LIDARAf(:)),bin1);
% yVec = linspace(min(IncAf(:)),max(IncAf(:)),bin1);
% pcolor(xVec,yVec,h); shading flat
% scatter(LIDARAf(:),IncAf(:))
boxplot(IncAf(:),bin)
h=findobj(gca,'tag','Outliers');
delete(h)
% hold on
grid on
% A = 0:50;
% plot(A,zeros(1,size(A,2)),'--k')
% colormap(flipud(gray))
% colorbar
xlabel('LIDAR Height (m)')
% ylabel('E[{\delta \tau}/{\delta \theta}|\theta]')
set(gca,'fontsize',12)
set(gca,'fontname','times new roman')
ylim([-0.1 0.5])
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1)
annotation('textbox',dim1,'String',str1,'fontsize',14,'Fontname','Times New Roman','EdgeColor','none')

% grid on
%     f = dir('*.fig')                                             ;
% for i = 1 : length(f)
%    fn = f(i).name                                                ;
%     openfig(fn)                                                  ;
%     fn = fn(1:length(fn)-4)                                      ;
%    print(fn,'-dtiff','-r300')                                    ;
% end
