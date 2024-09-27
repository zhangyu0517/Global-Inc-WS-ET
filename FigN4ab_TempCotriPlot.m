clc;clear all

%% 1984-2014 moving window
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\AI_Global8214_CRU05.mat'
load Contri_SM_8214MovWin_OBS.mat
load Contri_VPD_8214MovWin_OBS.mat

load Contri_SM_8214MovWin_CMIP_7x7_All.mat
load Contri_VPD_8214MovWin_CMIP_7x7_All.mat

Mod_ind=[2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,20,21,22,23,24];
Num_Mod=length(Mod_ind);%21; %模型数量
OBS_ind=[1,2,3];
Num_OBS=length(OBS_ind);%21; %观测数量

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\mask_ColdRegion.mat'
load Contri_VPD_8214_Products.mat
load Contri_SM_8214_Products.mat

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\Weight_Grid.mat'

lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';
years=[1982:1:2014]';

m=length(lon);n=length(lat);
Contri_SM_8214_mean=nanmean(Contri_SM_8214(:,:,OBS_ind),3);
Contri_VPD_8214_mean=nanmean(Contri_VPD_8214(:,:,OBS_ind),3);

SM_Dominance=zeros(m,n).*nan;
VPD_Dominance=zeros(m,n).*nan;

SM_Dominance(Contri_SM_8214_mean>Contri_VPD_8214_mean)=1; %1:SM主导；2：VPD主导
VPD_Dominance(Contri_SM_8214_mean<=Contri_VPD_8214_mean)=1;

% SM_Dominance(mask_ColdRegion==1)=nan;
% VPD_Dominance(mask_ColdRegion==1)=nan;
SM_Dominance(:,lat>60)=nan;
VPD_Dominance(:,lat>60)=nan;

%不计算cold region
Exclude_ColdRegion=ones(m,n);
% Exclude_ColdRegion(mask_ColdRegion==1)=nan;
Exclude_ColdRegion(:,lat>60)=nan;

%% 11年滑动
years_Mov=[1987:1:2009]';

lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';
m=length(lon);n=length(lat);

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\Weight_Grid.mat'

TemContri_SM_Dryland_OBS=zeros(length(years_Mov),Num_OBS).*nan;
TemContri_VPD_NonDryland_OBS=zeros(length(years_Mov),Num_OBS).*nan;

TemContri_SM_Dryland_MOD=zeros(length(years_Mov),Num_Mod).*nan;
TemContri_VPD_NonDryland_MOD=zeros(length(years_Mov),Num_Mod).*nan;


for i=1:length(years_Mov)
    for k1=1:Num_OBS
        % 面积加权
        d1=Contri_SM_8214MovWin_OBS(:,:,i,OBS_ind(k1)).*SM_Dominance;%.*Weight_Grid;
        TemContri_SM_Dryland_OBS(i,k1)=nansum(d1(:).*Weight_Grid(:))/nansum((isnan(d1(:))==0).*Weight_Grid(:));%nanmean(d1(:));
        

        nd1=Contri_VPD_8214MovWin_OBS(:,:,i,OBS_ind(k1)).*VPD_Dominance;%.*Weight_Grid;
        TemContri_VPD_NonDryland_OBS(i,k1)=nansum(nd1(:).*Weight_Grid(:))/nansum((isnan(nd1(:))==0).*Weight_Grid(:));%nanmean(nd1(:));
    end
    
    for k2=1:Num_Mod
        d3=Contri_SM_8214MovWin_CMIP(:,:,i,Mod_ind(k2)).*SM_Dominance;%.*Weight_Grid;
        TemContri_SM_Dryland_MOD(i,k2)=nansum(d3(:).*Weight_Grid(:))/nansum((isnan(d3(:))==0).*Weight_Grid(:));%nanmean(d3(:));%sum(domi_SM(:))/nansum(Dryland(:))*100;
        

        nd3=Contri_VPD_8214MovWin_CMIP(:,:,i,Mod_ind(k2)).*VPD_Dominance;%.*Weight_Grid;
        TemContri_VPD_NonDryland_MOD(i,k2)=nansum(nd3(:).*Weight_Grid(:))/nansum((isnan(nd3(:))==0).*Weight_Grid(:));%nanmean(nd3(:));%sum(domi_VPD(:))/nansum(NonDryland(:))*100;
    end
    
end

%% 去除空值的行
figure(1);h3=tight_subplot(1,2,[.1 .1],[.2 .2],[.18 .18]);
set(gcf,'position',[80 100 1400 300])
NA1={'(a)','(b)'};
rgbC =othercolor('Spectral10',14);

strm={'Supply-limited regions','Demand-limited regions'};
data_OBS_all=cat(3,TemContri_SM_Dryland_OBS,TemContri_VPD_NonDryland_OBS);
data_MOD_all=cat(3,TemContri_SM_Dryland_MOD,TemContri_VPD_NonDryland_MOD);


for k=1:2
    axes(h3(k)); 

   data_OBS=data_OBS_all(:,:,k);
   data_MOD=data_MOD_all(:,:,k);

    % OBS
    % 标注趋势值
    datain1=[years_Mov,mean(data_OBS,2)];
    [sig1,sen1]=ktaubZY(datain1,0.5);
    %画趋势线
    vv1 = median(datain1(:,2));
    middata1 = datain1(round(length(datain1)/2),1);
    slope1 = vv1 + sen1*(datain1(:,1)-middata1);
        
    f1=plot(datain1(:,1),datain1(:,2),'-o','color',[255 69 0]/255,'LineWidth',1.5);%[32/255,178/255,170/255]);
    hold on;
    f2=plot(datain1(:,1),slope1,'--','color',[255 69 0]/255,'LineWidth',1);%[199/255,21/255,133/255]);  
    hold on;
    
    % MOD
    % 标注趋势值
    datain2=[years_Mov,mean(data_MOD,2)];
    [sig2,sen2]=ktaubZY(datain2,0.5);
    %画趋势线
    vv2 = median(datain2(:,2));
    middata2 = datain2(round(length(datain2)/2),1);
    slope2 = vv2 + sen2*(datain2(:,1)-middata2);
        
    f3=plot(datain2(:,1),datain2(:,2),'-o','color',[120 120 120]/255,'LineWidth',1.5);%[32/255,178/255,170/255]);
    hold on;
    f4=plot(datain2(:,1),slope2,'--','color',[120 120 120]/255,'LineWidth',1);%[199/255,21/255,133/255]);  
    
     mark1= marksig(sig1,0.1,0.01);
     mark2= marksig(sig2,0.1,0.01);
     
     if k==1
      ylim([0.80 1.1])
      ytext=1.07;
    else
      ylim([0.65 0.95])
      ytext=0.92;
     end

         
     text(1986,ytext,strcat('OBS:',{32},num2str(sen1*10,'%.3f'),mark1,{32},'/decade'),'horiz','left','FontSize',14,'color',[255 69 0]/255);
     text(1986,ytext-0.04,strcat('CMIP6-MMM:',{32},num2str(sen2*10,'%.3f'),mark2,{32},'/decade'),'horiz','left','FontSize',14,'color',[120 120 120]/255);


    % xy轴标注
    xlim([1985 2011]);

    %加粗
    ax=gca;
    ax.LineWidth=1.2;
    set(gca,'TickDir','in')
    
    
    set(gca,'FontSize',14)
%     xlabel('Year','FontSize',16)
    if k==1 
        ylabel('\DeltaET(\it{SM}|\it{VPD})','FontSize',16)
        text(1981,1.15,NA1(k),'horiz','center','FontSize',20);
        
    else
        ylabel('\DeltaET(\it{VPD}|\it{SM})','FontSize',16)     
        text(1981,1,NA1(k),'horiz','center','FontSize',20);
    end
    
    title(strm(k),'FontSize',22)
    yt = get(gca,'YTick');    % 获取纵坐标轴刻度句柄
    
    ax=gca;
    ax.LineWidth=1.0;
  
end
print(gcf, '-dpng', '-r600','.\Figure4\Fig4ab.png');


%%
function res= marksig(sig,num1,num2)
if sig<=num1 & sig>num2
    res='*';
elseif sig<=num2
    res='**';
else
    res=' ';
end
    
end

% function res= marksig(sig,num1)
% if sig<=num1 
%     res='*';
% else
%     res=' ';
% end
%     
% end