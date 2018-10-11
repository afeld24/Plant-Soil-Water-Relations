clc
clear
load coast   
clat = lat;                                                        
clon = long ;

startyear = 2015; % First year
startmonth = 4; %First month
startday = 1; %Start day of first month
endyear = 2017; %Final year (if only doing one year, then startyear = endyear)
endmonth = 12;  %Final month
endday = 26; %End day of final month
% pathdata should include folder of raw HDF files in folder names of years
% (e.g. 2016)
t0 = datenum(startyear, startmonth, startday, 00, 00, 00);
tf = datenum(endyear, endmonth, endday, 00, 00, 00);
dt = datenum(0, 0, 1, 00, 00, 00);
date1 = datevec(t0:tf);
DateVector = date1(:,1:3);
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP (1)/MTDCA_fmincon')
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\9km Half Year')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP (1)\MTDCA_fmincon')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\9km Half Year')
% load('AfricaMTDCA_201504_201509','AfricaOmega')
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
load('MTDCAomegaMAT')
tlat =  [9.1];  %[-13.3, 27.8], [13.3, 27.8], [10.9, 27.2], [5 27.2],                                   
tlon =  [27.2];  %[9.1 27.2] %Not leaf out!! Look at [15, 20]
% 9.1, 27.2 is the official figure
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP Data/General Data')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\General Data')
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\General Data')
% load('MTDCA2_v13','IGBP_36km','water_f','lat','lon','IGBP_Names')
load('Ancillary_9km_3dB','IGBP_9km_3dB')
load('soil_texture_9km')
load('SMAPCenterCoordinates9KM')
% IGBP9 = nan(size(IGBP_9km_3dB,1),size(IGBP_9km_3dB,2));
for i = 1:size(IGBP_9km_3dB,3)
    A = IGBP_9km_3dB(:,:,i);
    A(A<0.75)=NaN;
    A(A>=0.75)=1;
    IGBP_9km_3dB(:,:,i) = A*i;
end
IGBP_9km_3dB = nanmean(IGBP_9km_3dB,3);

% AfricaRow = 300:1300;
% AfricaCol = 1700:2500;
IGBP_9km = IGBP_9km_3dB;
lon = SMAPCenterLongitudes;
lat = SMAPCenterLatitudes;                               

% Minimum Number of Overpasses in Drydown
dLength = 4;
 SmoothDay = 91; %odd number
% Included IGBP Classes and Regions
cIGBP = [  14  10   9  8  7  2  ];                        
IGBP_latlon(14,:) = [   9  15  -9  9 ];                            
IGBP_latlon(10,:) = [  10  17 -10 36 ];                           
IGBP_latlon( 9,:) = [ -22 -15 17 33  ];                            
IGBP_latlon( 8,:) = [ -15  -4 14 30  ];                            
IGBP_latlon( 7,:) = [ -32 -22 17 25  ];                            
IGBP_latlon( 2,:) = [  -4   5 10 30  ];                            

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

% Time-Series Example at tlat/tlon
                                       
klat = NaN(length(tlat),1);                              
klon = NaN(length(tlat),1);     
% Find EASE2 pixels closest to tlat tlon
for i = 1 : length(tlat)
[mm,latlon] = min(abs(lat(:)-tlat(i))+abs(lon(:)-tlon(i)));
[klat(i),klon(i)] = ind2sub(size(lat),latlon);                   
end
IGBP_9km(klat,klon)
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP Data\9km Half Year')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP Data\9km Half Year')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP (1)\MTDCA_fmincon')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP (1)/MTDCA_fmincon')
NPeriod = 1001;
% omega = AfricaOmega(klat,klon);
clayAf = clay(klat,klon);
sandAf = sand(klat,klon);
omega = MTDCAomega(klat,klon);

AMat1 = SMmat1.MTDCA_SM_201504_201509(klat,klon,:);
AMat2 = SMmat2.MTDCA_SM_201510_201603(klat,klon,:);
AMat3 = SMmat3.MTDCA_SM_201604_201609(klat,klon,:);
AMat4 = SMmat4.MTDCA_SM_201610_201703(klat,klon,:);
AMat5 = SMmat5.MTDCA_SM_201704_201709(klat,klon,:);
AMat6 = SMmat6.MTDCA_SM_201710_201712(klat,klon,:);
mv_DCA2 = squeeze(cat(3,AMat1,AMat2,AMat3,AMat4,AMat5,AMat6));
mv_DCA2 = mv_DCA2(1:NPeriod,1);
AMat1 = TAUmat1.MTDCA_TAU_201504_201509(klat,klon,:);
AMat2 = TAUmat2.MTDCA_TAU_201510_201603(klat,klon,:);
AMat3 = TAUmat3.MTDCA_TAU_201604_201609(klat,klon,:);
AMat4 = TAUmat4.MTDCA_TAU_201610_201703(klat,klon,:);
AMat5 = TAUmat5.MTDCA_TAU_201704_201709(klat,klon,:);
AMat6 = TAUmat6.MTDCA_TAU_201710_201712(klat,klon,:);
tau_DCA2 = squeeze(cat(3,AMat1,AMat2,AMat3,AMat4,AMat5,AMat6));
tau_DCA2 = tau_DCA2(1:NPeriod,1)./0.11;

% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\SMAP (1)\L1C_TB_E\L1C_TB_E_halfyear_Descending')
cd('/Users/andrewfeldman/Dropbox (MIT)/SMAP (1)/L1C_TB_E/L1C_TB_E_halfyear_Descending')
TBhAFileNam = 'TBh_L1C_9km_201504_201509';
TBhBFileNam = 'TBh_L1C_9km_201510_201603';
TBhCFileNam = 'TBh_L1C_9km_201604_201609';
TBhDFileNam = 'TBh_L1C_9km_201610_201703';
TBhEFileNam = 'TBh_L1C_9km_201704_201709';
TBhFFileNam = 'TBh_L1C_9km_201710_201712';

TBvAFileNam = 'TBv_L1C_9km_201504_201509';
TBvBFileNam = 'TBv_L1C_9km_201510_201603'; 
TBvCFileNam = 'TBv_L1C_9km_201604_201609';
TBvDFileNam = 'TBv_L1C_9km_201610_201703';
TBvEFileNam = 'TBv_L1C_9km_201704_201709';  
TBvFFileNam = 'TBv_L1C_9km_201710_201712'; 
TBhaMatNam = matfile(TBhAFileNam);
TBhbMatNam = matfile(TBhBFileNam); 
TBhcMatNam = matfile(TBhCFileNam);
TBhdMatNam = matfile(TBhDFileNam); 
TBheMatNam = matfile(TBhEFileNam); 
TBhfMatNam = matfile(TBhFFileNam); 
TBvaMatNam = matfile(TBvAFileNam);
TBvbMatNam = matfile(TBvBFileNam); 
TBvcMatNam = matfile(TBvCFileNam);
TBvdMatNam = matfile(TBvDFileNam); 
TBveMatNam = matfile(TBvEFileNam); 
TBvfMatNam = matfile(TBvFFileNam); 
TsaMatNam = matfile('Ts_AM_201504_201509_9km');
TsbMatNam = matfile('Ts_AM_201510_201603_9km');
TscMatNam = matfile('Ts_AM_201604_201609_9km');
TsdMatNam = matfile('Ts_AM_201610_201703_9km'); 
TseMatNam = matfile('Ts_AM_201704_201709_9km');
TsfMatNam = matfile('Ts_AM_201710_201712_9km');
Tsa = TsaMatNam.Ts_AM_201504_201509(klat,klon,:);
Tsb = TsbMatNam.Ts_AM_201510_201603(klat,klon,:);
Tsc = TscMatNam.Ts_AM_201604_201609(klat,klon,:);
Tsd = TsdMatNam.Ts_AM_201610_201703(klat,klon,:);
Tse = TseMatNam.Ts_AM_201704_201709(klat,klon,:);
Tsf = TsfMatNam.Ts_AM_201710_201712(klat,klon,:);

TBha = TBhaMatNam.TBh_L1C_9km_201504_201509(klat,klon,:);
TBhb = TBhbMatNam.TBh_L1C_9km_201510_201603(klat,klon,:);
TBhc = TBhcMatNam.TBh_L1C_9km_201604_201609(klat,klon,:);
TBhd = TBhdMatNam.TBh_L1C_9km_201610_201703(klat,klon,:);
TBhe = TBheMatNam.TBh_L1C_9km_201704_201709(klat,klon,:);
TBhf = TBhfMatNam.TBh_L1C_9km_201710_201712(klat,klon,:);
TBva = TBvaMatNam.TBv_L1C_9km_201504_201509(klat,klon,:);
TBvb = TBvbMatNam.TBv_L1C_9km_201510_201603(klat,klon,:); 
TBvc = TBvcMatNam.TBv_L1C_9km_201604_201609(klat,klon,:);
TBvd = TBvdMatNam.TBv_L1C_9km_201610_201703(klat,klon,:);
TBve = TBveMatNam.TBv_L1C_9km_201704_201709(klat,klon,:);
TBvf = TBvfMatNam.TBv_L1C_9km_201710_201712(klat,klon,:);
Ts = squeeze(cat(3,Tsa,Tsb,Tsc,Tsd,Tse,Tsf));

TBh = squeeze(cat(3,TBha,TBhb,TBhc,TBhd,TBhe,TBhf));
clear TBha TBhb TBhc TBhd TBhe TBhf
TBv = squeeze(cat(3,TBva,TBvb,TBvc,TBvd,TBve,TBvf));
clear TBva TBvb TBvc TBvd TBve TBvf
Ts = squeeze(cat(3,Tsa,Tsb,Tsc,Tsd,Tse,Tsf)); 
clear Tsa Tsb Tsc Tsd Tse Tsf

TBh = TBh(1:NPeriod,:);
TBv = TBv(1:NPeriod,:);
Ts = Ts(1:NPeriod,:);
startyear = 2015; % First year
startmonth = 4; %First month
startday = 1; %Start day of first month
endyear = 2018; %Final year (if only doing one year, then startyear = endyear)
endmonth = 3;  %Final month
endday = 31; %End day of final month
% pathdata should include folder of raw HDF files in folder names of years
% (e.g. 2016)
t0 = datenum(startyear, startmonth, startday, 00, 00, 00);
tf = datenum(endyear, endmonth, endday, 00, 00, 00);
dt = datenum(0, 0, 1, 00, 00, 00);
date1 = datevec(t0:tf);
DateVector = date1(:,1:3);
% Loop over each sample pixel
for i = 1 : length(klat)
    % Create sample tau/sm series and remove NaNs
    smt  = mv_DCA2;                  
    taut = tau_DCA2;   
% Year1Start = 1;
% Year1End = Year1Start+365;
% Year2Start = Year1End+1;
% Year2End = Year2Start+364;
% Year3Start = Year2End+1;
% Year3End = Year3Start+364;
    NPeriodSmooth = 366;
    MTDCAtauSmooth = nan(3,366);
    MTDCAtauSmooth(1,:) = taut(1:366);
    MTDCAtauSmooth(2,:) = [taut(367:700); NaN; taut(701:731)];
    MTDCAtauSmooth(3,:) = [taut(732:end); NaN(366-length(taut(732:end)),1)];    
    MTDCAtauSmooth = nanmean(MTDCAtauSmooth,1)';
    MTDCAtauSmooth1 = repmat(MTDCAtauSmooth,[3 1]);
    MTDCAtauSmooth2 = NaN(NPeriodSmooth,1);
%     k=1;
    Smooth1 = (SmoothDay-1)/2;
    for iDay = 1:NPeriodSmooth
        MTDCAtauSmooth2(iDay,1) = nanmean(MTDCAtauSmooth1(NPeriodSmooth+iDay...
            -Smooth1:NPeriodSmooth+iDay+Smooth1,1),1);
%         k = k+1;
    end
    MTDCAtauSmooth = repmat(MTDCAtauSmooth2,[3 1]);
    MTDCAtauSmooth = MTDCAtauSmooth(1:NPeriod);
    tauts =taut-MTDCAtauSmooth;
    Tbht = TBh;
    Tbvt = TBv;
    Tst = Ts;
     smt(end) = [] ; %remove ends for next part of analysis                                                  
    taut(end) = []; 
    tauts(end) = []; 
    Tbht(end) = [];
    Tbvt(end) = [];
    Tst(end) = [];

    tt = [1:length(taut)];                                 
    taut(isnan(smt)) = [];  
    tauts(isnan(smt)) = [];
    tt(isnan(smt)) = [];
    Tbht(isnan(smt)) = [] ;
    Tbvt(isnan(smt)) = [] ;  
    Tst(isnan(smt)) = [] ;     
    smt(isnan(smt)) = [];   
   taut(isnan(Tbvt)) = [];
   tauts(isnan(Tbvt)) = [];
  Tbht(isnan(Tbvt)) = [];   
   smt(isnan(Tbvt)) = []; 
   Tst(isnan(Tbvt)) = [];    
   Tbvt(isnan(Tbvt)) = [];

%    tautFourier = (taut-min(taut(:)))./(max(taut(:))-min(taut(:)));
%     smtFourier = (smt-min(smt(:)))./(max(smt(:))-min(smt(:)));
%     tautFourier = taut;
%     smtFourier = smt;
   
%         figure 
% set(0,'DefaultAxesFontSize' ,14)
% set(0,'defaultlinelinewidth', 0.5)
% set(0,'DefaultAxesFontName' ,'Arial')
% subplot(2,1,1)
% plot(tt,smt,'bo-',tt,taut,'go-')      
% xlabel('Days')
% title(['[Lat,Lon]=[' num2str(lat(klat(i),klon(i)),3) ',' ...
%          num2str(lon(klat(i),klon(i)),3) ']'])
% 
%   hold on 
    % Find all peaks in soil moisture
    ipb = find( smt(2:length(smt)-1) > smt(1:length(smt)-2) & ...
    smt(2:length(smt)-1) > smt(3:length(smt)));
    ipb = ipb + 1;  

    % Find all troughs in soil moisture
    ipe = find( smt(2:length(smt)-1) < smt(1:length(smt)-2) & ...
    smt(2:length(smt)-1) < smt(3:length(smt)));
    ipe = ipe + 1;  

    jj = 0;
    stauv = NaN(100,100);  
    stausv = NaN(100,100);                                      
    ssmv = NaN(100,100);                                      
    sttv = NaN(100,100);    
    sTBh = NaN(100,100);    
    sTBv = NaN(100,100);   
    sTs = NaN(100,100);     
    for k = 1 : length(ipb) %Cycle through each peak
    ke = min(find(ipe>ipb(k))); %Find the next bottom after the peak    
    TBhDry =  Tbht(ipb(k):ipe(ke)) ; %TBH between trough and peak                          
    TBvDry =  Tbvt(ipb(k):ipe(ke)); %TBV between trough and peak    
    TsDry =  Tst(ipb(k):ipe(ke));
    smv =  smt(ipb(k):ipe(ke)); %SM drydown                            
    tauv = taut(ipb(k):ipe(ke));  %tau drydown  
    tausv = tauts(ipb(k):ipe(ke));  %tau drydown 
    ttv =   tt(ipb(k):ipe(ke));  %tt drydown                         
    tauv(isnan(smv)) = [];   
    tausv(isnan(smv)) = [];  
    ttv(isnan(smv)) = [];                                             
    smv(isnan(smv)) = [];                                    
        if (length(smv)>=dLength)  %Only store if drydown length>=dLength                              

        jj = jj + 1                                         
        stauv(jj,1:length(smv)) = tauv;
        stausv(jj,1:length(smv)) = tausv;
        ssmv(jj,1:length(smv)) =  smv;                      
        sttv(jj,1:length(smv)) =  ttv;
        sTBh(jj,1:length(smv)) =  TBhDry;
        sTBv(jj,1:length(smv)) =  TBvDry;    
        sTs(jj,1:length(smv)) =  TsDry;
        end
    end        

        scase = find(~isnan(ssmv(:,1))) ; %SM Drydown vectors                          

         ncase = length(scase);  
Xticks = [62 154 245 336 428 520 611 701 793 885 976];
Xticks1 = {'6/1/15' '9/1/15' '12/1/15' '3/1/16' '6/1/16' '9/1/16' '12/1/16'...
    '3/1/17' '6/1/17' '9/1/17' '12/1/17'};
figure 
set(0,'DefaultAxesFontSize' ,18)

set(0,'defaultlinelinewidth', 1)
set(0,'DefaultAxesFontName' ,'Arial')

%%%% Top Figure
a = subplot(2,8,1:8)
p1 = plot(tt,smt,'bo-') 
ylabel('SM [m^{3} m^{-3}]')
set(gcf,'position',[100 100 1600 700])    
axis([0 NPeriod 0 max([taut;smt])])  
set(gca,'XTick',Xticks)
set(gca,'XTickLabel',Xticks1,'fontsize',18)
set(gca,'fontname','arial')                              
set(gcf,'PaperPositionMode','auto')
yyaxis right
hold on

p2 = plot(tt,taut,'go-') 
set(gca,'YColor','black')
yyaxis left
hold on 

for j = 1:ncase
    smv = ssmv(j,:);
    ttv = sttv(j,:);
    smv(isnan(smv))=[];
    ttv(isnan(ttv))=[];
    plot(ttv,smv,'r*')      
    hold on
    drawnow    
end
ylim([0 0.45])
yyaxis right
hold on
for j = 1:ncase
    tauv = stauv(j,:);
    ttv = sttv(j,:);
    tauv(isnan(tauv))=[];
    ttv(isnan(ttv))=[];
    plot(ttv,tauv,'k*')      
    hold on
    drawnow    
end

y1 = ylabel('VWC [kg m^{-2}]');
y1h = get(y1,'Position');
set(y1,'Rotation',270,'Position',[y1h(1)+30 y1h(2)])
legend([p1, p2], 'SM', 'VWC','Location','Northwest')     
    Xticks = [100 225 320 420];
    Xticks1 = {'7/9/15' '11/11/15' '2/14/16' '5/24/16'};  
kiVec = [1,4,7,9];

%%%% Bottom Left Figures
        for ki = 1 : 4
            if ki == 1
            
            elseif ki == 2
            elseif ki == 3
            else
            end
           sssttv =  sttv(kiVec(ki),:) ;                                      
            sssmv =  ssmv(kiVec(ki),:);                                       
          ssstauv = stauv(kiVec(ki),:);                                       
         sssttv(isnan(sssmv)) = [] ;                                   
        ssstauv(isnan(sssmv)) = [];                                    
          sssmv(isnan(sssmv)) = [];                                    
        h = subplot(2,8,8+ki)
        A = get(h,'position');
        set(h,'position',[A(1), A(2), A(3), A(4)+0.05])
        yyaxis left
        plot(sssttv,sssmv,'r* ',sssttv,sssmv,'bo-') 
        ylim([min(ssmv(:)) max(ssmv(:))])         
%         plot(sssttv,sssmv,'r* ',sssttv,ssstauv,'k* ', ...
%              sssttv,sssmv,'bo-',sssttv,ssstauv,'go-')       
        hold on
        yyaxis right
        plot(sssttv,ssstauv,'k* ',sssttv,ssstauv,'go-')      
%         ylim([min(stauv(:)) max(stauv(:))])
        ylim([0.35 max(stauv(:))])
        hold on
        axis([sssttv(1)-1 sssttv(end)+1 0  max([taut;smt]) ]) 
        set(gca,'fontname','arial')
        
        if ki == 1
            yyaxis left
            ylabel('SM [m^{3} m^{-3}]') 
            yyaxis right
            set(gca,'yticklabel',{})
        else
            yyaxis left
            set(gca,'yticklabel',{})
            yyaxis right
            set(gca,'yticklabel',{})            
        end
        yyaxis right
        ylim([0.35 max(stauv(:))])   
        set(gca,'YColor','black')
        yyaxis left
        ylim([0 0.45])    
        set(gca,'YColor','black')   
        
        drawnow
        set(gca,'XTick',Xticks(ki))
        set(gca,'XTickLabel',Xticks1(ki),'fontsize',18)
        end

%%% Bottom Right Figure
        
        h = subplot(2,8,13:16)
        for ki = 1 : length(scase)
           sssttv =  sttv(ki,:) ;                                      
            sssmv =  ssmv(ki,:);                                       
          ssstauv = stauv(ki,:);                                       
         sssttv(isnan(sssmv)) = [] ;                                   
        ssstauv(isnan(sssmv)) = [];                                    
          sssmv(isnan(sssmv)) = [];                                    
        
        plot(sssmv,ssstauv,'k.',sssmv,ssstauv,'k-')      
          xlabel('SM [m^{3} m^{-3}]')
          ylabel('VWC [kg m^{-2}]')
        hold on
        drawnow
        end
        ylim([0.35 max(stauv(:))])  
        xlim([min(ssmv(:)) max(ssmv(:))])          
        set(gca,'fontname','arial')
plotBottomRight = gca;
SatMatSand = -12.1;
SatMatClay = -40.5;
SatMatLoam = -47.8;
bSand = 4.05;
bClay = 11.4;
bLoam = 5.39;
nSand = 0.395;
nClay = 0.482;
nLoam = 0.451;

smVec1 = min(ssmv(:)):0.01:max(ssmv(:));
% smVec2 = 0.05:0.01:0.35;

loamAf = 1-clayAf-sandAf;
beffMat = bSand.*sandAf+bClay.*clayAf+bLoam.*loamAf;
neffMat = nSand.*sandAf+nClay.*clayAf+nLoam.*loamAf;
SateffMat = SatMatSand.*sandAf+SatMatClay.*clayAf+SatMatLoam.*loamAf;
MatricVec = ((SateffMat.*(smVec1./neffMat).^(-beffMat))/100).*(9810/10^6);

% x = 0:0.1:1;
x = linspace(0,1,size(smVec1,2));
MatricInterp = [-100 -10 -1 -0.1];
vq = interp1(MatricVec,x,MatricInterp) 


        plot([0.13 0.32],[1 0.6],'r','LineWidth',2)
        hold on
        plot([0.12 0.13],[1 1],'r','LineWidth',2)
        hold on
        plot([0.11 0.12],[0.98 1],'r','LineWidth',2)
        hold on
        plot([0.10 0.11],[0.94 0.98],'r','LineWidth',2)
        hold on
        plot([0.095 0.10],[0.9 0.94],'r','LineWidth',2)
        hold on
        plot([0.095 0.06],[0.9 0.6],'r','LineWidth',2)
        hold on
        plot([0.06 0.07],[0.6 0.75],'r','LineWidth',2)
        hold on
        plot([0.06 0.074],[0.6 0.65],'r','LineWidth',2)
        
% c1 = annotation('textarrow',[0.095 0.06],[0.9 0.6])
% c1.Color = 'red';
% c1.LineWidth = 2;   
c2 = text(0.08,0.6,'Dry Down Progression in Time','Fontsize',18,'FontName','arial')
c2.Color = 'red';

A = get(h,'position');
set(h,'position',[A(1)+0.04, A(2), A(3)-0.04, A(4)+0.05])    
    set(gcf,'color','white')    
 A1 = plotBottomRight.Position;
ax2 = axes('Position',A1,...
    'XAxisLocation','top','Color','none');
set(gca,'YTickLabel','')
set(gca,'YTick','')
set(gca,'XTick',vq)
set(gca,'XTickLabel',{'-100','-10','-1','-0.1'})
xlabel('Soil Matric Potential [MPa]')
set(gca,'fontname','arial')
% set(c2,'Rotation',352)
% c4 = text(0.062,0.87,'VWC Decreasing', 'Fontsize',14)
% c4.Color = 'red';
% set(c4,'Rotation',42)

        % subplot(2,8,13:14)
% hold on
% cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')   
        
%         h = subplot(2,8,15:16)
%         A = get(h,'position');
%         set(h,'position',[A(1)+0.04, A(2), A(3)-0.04, A(4)+0.05])
%         
%         for ki = 1 : length(scase)
%         sssttv =  sttv(ki,:) ;                                      
%         sssmv =  ssmv(ki,:);                                       
%         ssstauv = stausv(ki,:);                                       
%         sssttv(isnan(sssmv)) = [] ;                                   
%         ssstauv(isnan(sssmv)) = [];                                    
%           sssmv(isnan(sssmv)) = [];   
%           
%         plot(sssmv,ssstauv,'k.',sssmv,ssstauv,'k-')
% %         plot(sssmv,ssstauv,'k*',sssmv,ssstauv,'k-')      
%           xlabel('SM [m^{3} m^{-3}]')
%           ylabel('\tau'' [Np]')
%         hold on
% %         axis([sssttv(1)-1 sssttv(end)+1 0  max([taut;smt]) ]) 
%         grid on
%         drawnow
%         ylim([min(stausv(:)) max(stausv(:))])
%         xlim([min(ssmv(:)) max(ssmv(:))])    
%         set(gca,'fontname','times new roman')
%         end
end

%Reposition Tiles
A = get(a,'Position'); % top fig
set(a,'Position',[A(1) A(2)+0.06 A(3) A(4)]);

annotation('textbox',[0.07 0.9 0.1 0.1],'String','a','fontsize',20,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.53 0.47 0.1 0.1],'String','c','fontsize',20,'Fontname','arial','FontWeight','bold','EdgeColor','none')
annotation('textbox',[0.07 0.47 0.1 0.1],'String','b','fontsize',20,'Fontname','arial','FontWeight','bold','EdgeColor','none')

cd('/Users/andrewfeldman/Dropbox (MIT)/MIT/Manuscripts/WaterExchangeNatureGeoscience/FinalSubmission_NaturePlants/PaperFigures')
print('Fig1','-depsc','-r1200') % save as high res PNGd    
%%
D = get(d,'Position'); % bottom left fig
set(d,'Position',[D(1) D(2) D(3) D(4)]);
% B = get(b,'Position'); % bottom right fig
% set(b,'Position',[B(1)+0.015 B(2)+0.02 B(3)+0.007 B(4)+0.06]);
B = get(b,'Position'); % bottom right fig
set(b,'Position',[D(1) A(2) D(3) D(4)]);

% taudiff = diff(stauv,[],2);
% smdiff = diff(ssmv,[],2);
% ADiff = taudiff./smdiff;
% ADiff = ADiff(:);
% ADiff(isnan(ADiff))=[];
% [f x] = ksdensity(ADiff,'bandwidth',0.01);
% plot(x,f)
% median(ADiff)
%%
hold on
cd('C:\Users\afeldmangolf24\Dropbox (MIT)\SMAP\Project_PlantSoilVODStudy')
% cd('C:\Users\Entekhabi-Group\Dropbox (MIT)\MIT\Manuscripts\WaterExchangeNatureGeoscience')
sTBh(ssmv>0.12)=NaN;
sTBv(ssmv>0.12)=NaN;
sTs(ssmv>0.12)=NaN;
TBh1 = nanmean(sTBh(:));
TBv1 = nanmean(sTBv(:));
Ts1 = nanmean(sTs(:));
JRes = 1000;
h = 0.13;
tau = flipud(repmat(linspace(0,1.5,JRes),[JRes 1])');
theta = repmat(linspace(0,0.6,JRes),[JRes 1]);
angle = 40;
TBhMPix = ones(JRes,JRes).*TBh1;
TBvMPix = ones(JRes,JRes).*TBv1;
TsMPix = ones(JRes,JRes).*Ts1;     
ClayPix = ones(JRes,JRes).*clayAf;
omgPix = ones(JRes,JRes).*omega;
y = exp(-tau./cosd(angle));
[k] = real(mironov(1.4,theta,ClayPix));
x = sqrt(secd(angle).*secd(angle).*k -tand(angle).*tand(angle));
rHr = ((1-x).^2)./((1+x).^2).*exp(-h);                                    
rVr = ((k-x).^2)./((k+x).^2).*exp(-h);
TBhMod = (TsMPix.*(1-rHr).*y) + ...
         (TsMPix.*(1-omgPix).*(1-y).*(1+(rHr.*y))) ;
TBvMod = (TsMPix.*(1-rVr).*y) + ...
         (TsMPix.*(1-omgPix).*(1-y).*(1+(rVr.*y))) ;     
JPix1 = ((TBhMPix-TBhMod).^2)+((TBvMPix-TBvMod).^2);

SMTAUMin1 = nan(size(JPix1,1),2);
for i = 1:size(JPix1,1)
%     A = ;
    [a b] = find(JPix1(i,:)==min(JPix1(i,:)));
    SMTAUMin1(i,:) = [tau(i,1) theta(a,b)];
end
plot(SMTAUMin1(:,2),SMTAUMin1(:,1),'-r')
% ylim([0 max(smt(:))])
set(gcf,'color','white')

%%
JPix1(JPix1>300)=NaN;
contour(theta,tau,JPix1)
caxis([0 100])
