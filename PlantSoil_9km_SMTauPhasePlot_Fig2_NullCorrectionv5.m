clc
clear
load coast   
clat = lat;                                                        
clon = long ;  
SmoothDay = 91; %odd number
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\9km Half Year')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\9km Half Year')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP (1)/MTDCA_fmincon')
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

cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP Data/General Data')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\General Data')
load('Ancillary_9km_3dB','IGBP_9km_3dB','IGBP_Names')
load('SMAPCenterCoordinates9KM')
for i = 1:size(IGBP_9km_3dB,3)
    A = IGBP_9km_3dB(:,:,i);
    A(A<0.75)=NaN;
    A(A>=0.75)=1;
    IGBP_9km_3dB(:,:,i) = A*i;
end
IGBP_9km_3dB = nanmean(IGBP_9km_3dB,3);

AfricaRow = 300:1300;
AfricaCol = 1700:2500;
NPeriod = 1001;
IGBP_9km = IGBP_9km_3dB(AfricaRow,AfricaCol);
lon = SMAPCenterLongitudes(AfricaRow,AfricaCol);
lat = SMAPCenterLatitudes(AfricaRow,AfricaCol); 

% Minimum Number of Overpasses in Drydown
dLength = 4;

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
IGBPx = NaN(size(IGBP_9km));                              
for i = 1 : length(cIGBP)
IGBPz = NaN(size(IGBP_9km));   
%Finding pixels within bounds
dlatlon = find(lat>=IGBP_latlon(cIGBP(i),1) & ...
               lat<=IGBP_latlon(cIGBP(i),2) & ...
               lon>=IGBP_latlon(cIGBP(i),3) & ...
               lon<=IGBP_latlon(cIGBP(i),4) & ...
               IGBP_9km==cIGBP(i));
IGBPx(dlatlon) = cIGBP(i);                                             
IGBPz(dlatlon) = cIGBP(i);                                                     
end

Np = 600000;                                                                                                        
sminterp = [ 0.01 : 0.01 : 0.60 ];   
Ni = length(sminterp);
IGBPAdjust = [10 8 7 14 9 2];
VecPair = NaN(2,2,Ni,Np); 
% VecPairSeas = NaN(2,2,Ni,Np); 

% for ki = 1 : length( IGBPAdjust) % Loop over each IGBP class
% for ki = 1 : 2 % Loop over each IGBP class
count = 0;

for ki = 3:4
    count = count+1;
     sprintf(['Storing in ' num2str(ki)])
     ji = 0;                                                       
        IGBPz = NaN(size(IGBPx));                                  
        IGBPz(find(IGBPx== IGBPAdjust(ki))) = 1;   %Only pixels in ki IGBP class                      
        [klat,klon] = find(~isnan(IGBPz));  %EASE2 location of pixels                       
kRow = min(klat):max(klat);
kCol = min(klon):max(klon);

% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\9km Half Year')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP (1)/MTDCA_fmincon')

% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\9km Half Year')
AMat1 = SMmat1.MTDCA_SM_201504_201509(AfricaRow(kRow),AfricaCol(kCol),:);
AMat2 = SMmat2.MTDCA_SM_201510_201603(AfricaRow(kRow),AfricaCol(kCol),:);
AMat3 = SMmat3.MTDCA_SM_201604_201609(AfricaRow(kRow),AfricaCol(kCol),:);
AMat4 = SMmat4.MTDCA_SM_201610_201703(AfricaRow(kRow),AfricaCol(kCol),:);
AMat5 = SMmat5.MTDCA_SM_201704_201709(AfricaRow(kRow),AfricaCol(kCol),:);
AMat6 = SMmat6.MTDCA_SM_201710_201712(AfricaRow(kRow),AfricaCol(kCol),:);
mv_DCA2 = squeeze(cat(3,AMat1,AMat2,AMat3,AMat4,AMat5,AMat6));
AMat1 = TAUmat1.MTDCA_TAU_201504_201509(AfricaRow(kRow),AfricaCol(kCol),:);
AMat2 = TAUmat2.MTDCA_TAU_201510_201603(AfricaRow(kRow),AfricaCol(kCol),:);
AMat3 = TAUmat3.MTDCA_TAU_201604_201609(AfricaRow(kRow),AfricaCol(kCol),:);
AMat4 = TAUmat4.MTDCA_TAU_201610_201703(AfricaRow(kRow),AfricaCol(kCol),:);
AMat5 = TAUmat5.MTDCA_TAU_201704_201709(AfricaRow(kRow),AfricaCol(kCol),:);
AMat6 = TAUmat6.MTDCA_TAU_201710_201712(AfricaRow(kRow),AfricaCol(kCol),:);
tau_DCA2 = squeeze(cat(3,AMat1,AMat2,AMat3,AMat4,AMat5,AMat6));  
tau_DCA2 = tau_DCA2./0.11;
tau_DCA2mean = nanmean(tau_DCA2,3);

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

tau_DCA2 = tau_DCA2(:,:,1:NPeriod);
mv_DCA2 = mv_DCA2(:,:,1:NPeriod);
% MTDCAtauSmooth = NaN(size(tau_DCA2));
% Smooth1 = (SmoothDay-1)/2;
% for iDay = Smooth1+1:NPeriod-Smooth1
%     MTDCAtauSmooth(:,:,iDay) = nanmean(tau_DCA2(:,:,iDay-Smooth1:iDay+Smooth1),3);
% end
% tau_DCA2_Smooth = MTDCAtauSmooth;
  klat = klat-min(klat)+1;
  klon = klon-min(klon)+1;
  cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP/Project_PlantSoilVODStudy')
%   cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
%   cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
for i = 1 : length(klat) % Loop over each pixel
     smt  = squeeze( mv_DCA2(klat(i),klon(i),:));                  
     taut = squeeze(tau_DCA2(klat(i),klon(i),:));   
    tauts = squeeze(tau_DCA2_Smooth(klat(i),klon(i),:));
    smts = squeeze(mv_DCA2_Smooth(klat(i),klon(i),:));

%    plot(1:1006,tauts,'-b',1:1006,taut,'ob' )
 smt(end) = [];                                                    
taut(end) = []; 
%  smts(end) = [];
tauts(end) = []; 
tt = [ 1 : length(taut) ]; 
taut(isnan(smt)) = []; 
tauts(isnan(smt)) = [];
% smts(isnan(smt)) = [];
tt(isnan(smt)) = [];                                       
smt(isnan(smt)) = [];

% Original
  [NDry,timevO,smvO,tauvO] = DryDowns(tt,smt,taut,dLength); 
% Seasonal       
  [NDry,timevS,smvS,tauvS] = DryDowns(tt,smt,tauts,dLength);
  
% Populate VPxxx DataBase and Interpolate DryDowns on Common 
% Soil Moisture Scale So Events Can be Average (Padded With NaN)
for j = 1 : NDry
ji = ji+1                                        ;

if (ji>Np)
   disp(['Error: Increase Np Events Count Max'])
end
% Original       
psm  = cell2mat( smvO(j));                                    
ptau = cell2mat(tauvO(j));
psmS  = cell2mat( smvS(j));                                    
ptauS = cell2mat(tauvS(j));
ptt = cell2mat(timevO(j));
ptts = cell2mat(timevS(j)); 
ptauSfinite = ptauS;
ptauSfinite(isfinite(ptauSfinite))=1;
ptauSfinite(isnan(ptauSfinite))=0;
if nansum(ptauSfinite)<=3 % Remove smoothed dry downs on edges of time series
else
psmSeas = smts(ptt);
ptau(isnan(psm)) = []; 
psmSeas(isnan(psm)) = []; 
psm(isnan(psm)) = []; 
ptauS(isnan(psmS)) = [];                                            
psmS(isnan(psmS)) = [];
SlopeNulltau = diff(ptauS);
SlopeRegtau = diff(ptau);
SlopeDifftau = SlopeRegtau-SlopeNulltau;
ptauC = nan(size(ptau));
ptauC(1) = ptau(1);
SlopeNullsm = diff(psmSeas);
SlopeRegsm = diff(psm);
SlopeDiffsm = SlopeRegsm-SlopeNullsm;
psmC = nan(size(psm));
psmC(1) = psm(1);
for a = 2:size(SlopeDifftau,1)+1
    ptauC(a) = ptauC(a-1)+SlopeDifftau(a-1);
    psmC(a) = psmC(a-1)+SlopeDiffsm(a-1);    
end
% plot(psm,ptau,'-b',psm,ptauC,'-g')
% Remove null slope from real slope
tauinterp = interp1(psmC,ptauC,sminterp,'linear',NaN);           
smin = sminterp;                                          
smin(isnan(tauinterp)) = NaN;                                      
VecPair(count,1,:,ji) = smin;                              
VecPair(count,2,:,ji) = tauinterp;                              

% Seasonal                                                                                
% tauinterp = interp1(psm,ptau,sminterp,'linear',NaN);           
% smin = sminterp ;                                         
% smin(isnan(tauinterp)) = NaN ;                                     
% VecPairSeas(count,1,:,ji) = smin;                              
% VecPairSeas(count,2,:,ji) = tauinterp;     
end
end % End Drydown loop within single pixel
end %End single pixel
end %End IGBP

%% Plot
% xtauMax = [ 0.5 0.4 0.5 0.8 0.4 1.3 ];                             
% xtauMin = [ 0.0 0.0 0.1 0.2 0.0 0.5 ];                             
%  xsmMax = [ 0.6 0.6 0.6 0.6 0.6 0.6 ] ;                            
%  xsmMin = [ 0.0 0.0 0.0 0.0 0.0 0.0 ] ;  
 xtauMax = [ 3  7.5 4 3.5 4.5 15 ]; 
% xtauMax = [ 3  4 4 3.5 4.5 15 ];
xtauMin = [ 0.0 1.5 0 0 1 5 ];                             
 xsmMax = [ 0.6 0.5 0.4 0.45 0.5 0.6 ] ;                            
 xsmMin = [ 0.0 0.0 0.0 0.0 0.0 0.0 ] ; 
% Number of Noodle Origins (SM = Nxsm, tau = NXtau)
   Nxsm = 12;                                                     
  Nxtau = 12;                                                    
% Minimum Number of Drydowns to Form a Noode]le
NcaseMin= 100;                                                       
set(0, 'DefaultAxesFontSize' ,22     )
set(0, 'defaultlinelinewidth', 1.0   )
set(0, 'DefaultAxesFontName' ,'arial')
Ncb = 300;
Ncr = 100;
red = [ ones(Ncr,1)    [1:Ncr]'/Ncr  [1:Ncr]'/Ncr ];               
blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];               
slopeDtauDsm = 3;   
figure
count = 0;
for ki = 3:4
    count = count+1;
% for ki = 1 : 1
    if count ==1
       a = subplot(1,2,1)
    else
       b = subplot(1,2,2)
    end

xtau = xtauMin(ki) + ((xtauMax(ki) - xtauMin(ki))*[1:Nxtau ]/(Nxtau )) ;  
xsm =  xsmMin(ki) + (( xsmMax(ki) -  xsmMin(ki))*[1:Nxsm   ]/(Nxsm )) ;  

bsm = squeeze(VecPair(count,1,:,:));                        
[bsm,ibsm] = nanmax(bsm,[],1);                              
bbtau = squeeze(VecPair(count,2,:,:)) ;                       
for i = 1 : length(ibsm)                 
btau(i) = bbtau(ibsm(i),i);                      
end
bbtaufinite = bbtau;
bbtaufinite(isfinite(bbtaufinite))=1;
bbtauMean = nanmedian(bbtau,2);
bbtauCount = nansum(bbtaufinite,2);
bbtauMean(bbtauCount<100)=NaN;
% Step1: Plot lines first: for All Noodle Origin Points              
for ism = 1 : Nxsm-1
    for itu = 1 : Nxtau-1
  
% Find Cases where Drydown Origin (bsm,btau) is Within Box              
            icase = find( (bsm > xsm(ism)) & (bsm < xsm(ism+1)) & ...               
                          (btau>xtau(itu)) & (btau<xtau(itu+1)) ) ;

                      psm = squeeze(VecPair(count,1,:,icase))             ;
                     ptau = squeeze(VecPair(count,2,:,icase))             ;                    

% Count How Many Drydowns at each Uniform Discretization Level                     
                  isbbtau = sum(~isnan(ptau),2)                   ;
          
% Do Not Include In Noodle Statistics if Less than NcaseMin                  
                     ptau(find(isbbtau<NcaseMin),:) = NaN         ;                     
                      psm = nanmedian( psm,2)                     ;    
                     ptau = nanmedian(ptau,2)                     ;   
         psm(isnan(ptau)) = []                                    ; 
        ptau(isnan(ptau)) = []                                    ; 
        
% Plot Noodle if at Least With Given Length         
   if (length(ptau)>3)
             plot(psm,ptau,'-k','linewidth',1.5,'Color',[0.5 0.5 0.5])         
%                   slp = diff(ptau)./diff(psm)                     ;
%                   slp = nanmedian(slp)                            ;
%          scatter(psm(end),ptau(end),65,slp,'filled','MarkerEdgeColor',[0 0 0])
               hold on
              
%               caxis([-slopeDtauDsm slopeDtauDsm*3])

              colormap([red;flipud(blu)])
%               drawnow
              xlabel('Soil Moisture [m^{3} m^{-3}]')
              ylabel('Vegetation Water Content [kg m^{-2}]')
              xlim([  xsmMin(ki)  xsmMax(ki) ]) 
              ylim([ xtauMin(ki) xtauMax(ki) ])             
              grid on
              if ki == 1
                  title(['Sahelo-Soudan ' char(IGBP_Names(IGBPAdjust(ki)))])
              elseif ki == 2
                  title(['South Central African ' char(IGBP_Names(IGBPAdjust(ki)))])
              else
              title(IGBP_Names(IGBPAdjust(ki)))
              end
% if count == 2
% % h = colorbar;
% % set(get(h,'title'),'string','\DeltaVWC/\DeltaSM')
% c1 = colorbar;
% pos = get(c1,'Position')
% c1.Direction = 'reverse';
% % colormap(flipud(parula))
% c1.Label.String = '\DeltaVWC/\DeltaSM [(kg m^{-2})/(m^{3} m^{-3})]';
% c1.Label.Rotation = 270;
% c1.Label.Position = [pos(1)+4.4 pos(2)-1];
% end
   end % Length if length is long enough
  end % itu
end % ism

% Step2: Plot dots on top: for All Noodle Origin Points       
for ism = 1 : Nxsm-1
    for itu = 1 : Nxtau-1
  
% Find Cases where Drydown Origin (bsm,btau) is Within Box              
            icase = find( (bsm > xsm(ism)) & (bsm < xsm(ism+1)) & ...               
                          (btau>xtau(itu)) & (btau<xtau(itu+1)) ) ;

                      psm = squeeze(VecPair(count,1,:,icase))             ;
                     ptau = squeeze(VecPair(count,2,:,icase))             ;                    

% Count How Many Drydowns at each Uniform Discretization Level                     
                  isbbtau = sum(~isnan(ptau),2)                   ;
          
% Do Not Include In Noodle Statistics if Less than NcaseMin                  
                     ptau(find(isbbtau<NcaseMin),:) = NaN         ;                     
                      psm = nanmedian( psm,2)                     ;    
                     ptau = nanmedian(ptau,2)                     ;   
         psm(isnan(ptau)) = []                                    ; 
        ptau(isnan(ptau)) = []                                    ; 
        
% Plot Noodle if at Least With Given Length         
   if (length(ptau)>3)
%              plot(psm,ptau,'-k','linewidth',1,'Color',[0.5 0.5 0.5])         
                  slp = diff(ptau)./diff(psm)                     ;
                  slp = nanmedian(slp)                            ;
         scatter(psm(end),ptau(end),100,slp,'filled','MarkerEdgeColor',[0 0 0])
               hold on
              
              caxis([-slopeDtauDsm slopeDtauDsm*3])

              colormap([red;flipud(blu)])
%               drawnow
              xlabel('Soil Moisture [m^{3} m^{-3}]')
              ylabel('Vegetation Water Content [kg m^{-2}]')
              xlim([  xsmMin(ki)  xsmMax(ki) ]) 
              ylim([ xtauMin(ki) xtauMax(ki) ])             
              grid on
              if ki == 1
                  title(['Sahelo-Soudan ' char(IGBP_Names(IGBPAdjust(ki)))])
              elseif ki == 2
                  title(['South Central African ' char(IGBP_Names(IGBPAdjust(ki)))])
              else
              title(IGBP_Names(IGBPAdjust(ki)))
              end
if count == 2
% h = colorbar;
% set(get(h,'title'),'string','\DeltaVWC/\DeltaSM')
c1 = colorbar;
pos = get(c1,'Position');
c1.Direction = 'reverse';
% colormap(flipud(parula))
% c1.Label.String = '\DeltaVWC/\DeltaSM';
c1.Label.String = '\DeltaVWC/\DeltaSM [(kg m^{-2})/(m^{3} m^{-3})]';
c1.Label.Rotation = 270;
c1.Label.Position = [pos(1)+4.4 pos(2)+3];
end
   end % Length if length is long enough
  end % itu
end % ism

%      hold on
%      plot(sminterp,bbtauMean,'-k','LineWidth',2)
end % ki

if ki ==5
dim1 = [0.035 0.91 0.1 0.1];
str1 = 'c';
    annotation('textbox',dim1,'String',str1,'fontsize',22,'FontWeight','bold',...
        'Fontname','arial','EdgeColor','none')   
    set(gcf,'Position',[100 100 1500 600])
A = get(a,'Position');
set(a,'Position',[A(1)-0.05 A(2)+0.02 A(3)+0.04 A(4)]);
cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
print('FigS2cSMVWCPhase','-djpeg','-r300')
B = get(b,'Position');
set(b,'Position',[B(1)-0.04 B(2)+0.02 B(3)+0.16 B(4)]);
set(gcf,'color','white')

else
dim1 = [0.035 0.91 0.1 0.1];
str1 = 'a';
    annotation('textbox',dim1,'String',str1,'fontsize',22,'FontWeight','bold',...
        'Fontname','arial','EdgeColor','none')
dim2 = [0.495 0.91 0.1 0.1];
str2 = 'b';
    annotation('textbox',dim2,'String',str2,'fontsize',22,'FontWeight','bold',...
        'Fontname','arial','EdgeColor','none')
set(gcf,'Position',[100 100 1500 600])
A = get(a,'Position');
set(a,'Position',[A(1)-0.05 A(2)+0.02 A(3)+0.04 A(4)]);
B = get(b,'Position');
set(b,'Position',[B(1)-0.04 B(2)+0.02 B(3)+0.11 B(4)]);
set(gcf,'color','white')

end
%
% Plot Seasonal Mean (Null)
% % xtauMax = [ 0.5 0.4 0.5 0.8 0.4 1.3 ];                             
% % xtauMin = [ 0.0 0.0 0.1 0.2 0.0 0.5 ];                             
% %  xsmMax = [ 0.6 0.6 0.6 0.6 0.6 0.6 ] ;                            
% %  xsmMin = [ 0.0 0.0 0.0 0.0 0.0 0.0 ] ;  
 xtauMax = [ 3  8 5 4 5 15 ];                             
xtauMin = [ 0.0 1 1 0 1 5 ];                                   
 xsmMax = [ 0.6 0.6 0.6 0.6 0.6 0.6 ] ;                            
 xsmMin = [ 0.0 0.0 0.0 0.0 0.0 0.0 ] ; 
cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
if ki ==2
% print('Fig2SMVWCPhase','-djpeg','-r300') % save as high res PNGd   
print('Fig2','-depsc','-r1200') % save as high res PNGd    
elseif ki ==4
% print('FigS2abSMVWCPhase','-djpeg','-r300') % save as high res PNGd  
print('FigS2ab','-depsc','-r1200') % save as high res PNGd    
elseif ki ==5
print('FigS2c','-depsc','-r1200') % save as high res PNGd    
% print('FigS2cSMVWCPhase','-djpeg','-r300') % save as high res PNGd    
end
% % Number of Noodle Origins (SM = Nxsm, tau = NXtau)
%    Nxsm = 12;                                                     
%   Nxtau = 12;                                                    
% % Minimum Number of Drydowns to Form a Noode]le
% NcaseMin= 100;                                                       
% set(0, 'DefaultAxesFontSize' ,14     )
% set(0, 'defaultlinelinewidth', 1.0   )
% set(0, 'DefaultAxesFontName' ,'Arial')
% Ncb = 100;                                                         
% red = [ ones(Ncb,1)    [1:Ncb]'/Ncb  [1:Ncb]'/Ncb ];               
% blu = [ [1:Ncb]'/Ncb   [1:Ncb]'/Ncb  ones(Ncb,1)  ];               
%               slopeDtauDsm = 5; 
% figure
% for ki = 1 : 2 

% for ki = 1 : 2
%     if ki == 1
%     subplot(1,2,1)
%     axes('Position',[.3 .6 .15 .15])
%     set(gca,'fontname','times new roman')
%     
%     else
%     subplot(1,2,2)
%     axes('Position',[.75 .6 .15 .15])   
%     set(gca,'fontname','times new roman')
%     end
% end
%%
hold on
%
count = 0;
for ki = 5 : 6
    count = count+1;
    if count == 1
%     subplot(1,2,1)
    axes('Position',[.31 .61 .14 .28])
    else
%     subplot(1,2,2)
    axes('Position',[.71 .61 .14 .28])   
%     set(gca,'FontSize',10)        
    end

xtau = xtauMin(ki) + ((xtauMax(ki) - xtauMin(ki))*[1:Nxtau ]/(Nxtau )) ;  
xsm =  xsmMin(ki) + (( xsmMax(ki) -  xsmMin(ki))*[1:Nxsm   ]/(Nxsm )) ;  

bsm = squeeze(VecPairSeas(count,1,:,:));                        
[bsm,ibsm] = nanmax(bsm,[],1);                              
bbtau = squeeze(VecPairSeas(count,2,:,:)) ;                       
for i = 1 : length(ibsm)                 
btau(i) = bbtau(ibsm(i),i);                      
end

% For All Noodle Origin Points              
     for ism = 1 : Nxsm-1
          for itu = 1 : Nxtau-1
  
% Find Cases where Drydown Origin (bsm,btau) is Within Box              
            icase = find( (bsm > xsm(ism)) & (bsm < xsm(ism+1)) & ...               
                          (btau>xtau(itu)) & (btau<xtau(itu+1)) ) ;

                      psm = squeeze(VecPairSeas(count,1,:,icase));             
                     ptau = squeeze(VecPairSeas(count,2,:,icase));                                 

% Count How Many Drydowns at each Uniform Discretization Level                     
                  isbbtau = sum(~isnan(ptau),2)                   ;
          
% Do Not Include In Noodle Statistics if Less than NcaseMin                  
                     ptau(find(isbbtau<NcaseMin),:) = NaN         ;                     
                      psm = nanmedian( psm,2)                     ;    
                     ptau = nanmedian(ptau,2)                     ;   
         psm(isnan(ptau)) = []                                    ; 
        ptau(isnan(ptau)) = []                                    ; 
        
% Plot Noodle if at Least With Given Length         
           if (length(ptau)>3)
             plot(psm,ptau,'-k','linewidth',1,'Color',[0.5 0.5 0.5])         
                  slp = diff(ptau)./diff(psm)                     ;
                  slp = nanmedian(slp)                            ;
         scatter(psm(end),ptau(end),65,slp,'filled','MarkerEdgeColor',[0 0 0])
               hold on
              
              caxis([-slopeDtauDsm slopeDtauDsm])
%                              h = colorbar;
%                set(get(h,'title'),'string','\DeltaVWC/\DeltaSM')
              colormap([red;flipud(blu)])
%               drawnow
%               xlabel('Soil Moisture \theta [m^{3} m^{-3}]')
%               ylabel('Vegetation Water Content [kg m^{-2}]')
                  set(gca,'FontSize',10)
              xlim([  xsmMin(ki)  xsmMax(ki) ]) 
              ylim([ xtauMin(ki) xtauMax(ki) ])             
              grid on
              title('Seasonal Mean')

           end % Length if
          end % itu
     end % ism
end % ki
% dim1 = [0.38 0.09 0.1 0.1];
% str1 = 'A';
%     annotation('textbox',dim1,'String',str1,'fontsize',18,...
%         'Fontname','Times New Roman','EdgeColor','none')
% dim2 = [0.83 0.09 0.1 0.1];
% str2 = 'B';
%     annotation('textbox',dim2,'String',str2,'fontsize',18,...
%         'Fontname','Times New Roman','EdgeColor','none')
% set(gcf,'Position',[100 100 1500 600])
