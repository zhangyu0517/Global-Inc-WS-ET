%%
clc;clear all

%%
lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';

years=[1982:1:2014]';

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\ET_Ecology05_8214.mat'
load 'F:\ZhangYu\Global_ET_SMVPD\Code\ET_Global8216_GLEAM05.mat'
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\ET_BESS05_8214.mat'
% load 'F:\ZhangYu\Global_ET_SMVPD\Code\ET_PMLv2_8214.mat'
% load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\ET_REA05_8214.mat'

% load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\SM_Merge05_8214.mat' 
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\SM100_Global8214_ERA5Land.mat'
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\VPD_Global8214_CRU05.mat'

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN2\mask_Global05.mat'

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\GS_start_gobal05.mat' 
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\GS_end_gobal05.mat' 

load GrowSeason_gobal05.mat
%%
% 
m=length(lon);n=length(lat);

mask=isnan(mean(ET_Ecology05_8214,3))==0;
% mask=mask_Global05.*(AI_Gloablv3_05>=0.05);

Z_Season_et= zeros(m,n,length(years)*12,3).*nan;
Z_Season_sm= zeros(m,n,length(years)*12).*nan;
Z_Season_vpd= zeros(m,n,length(years)*12).*nan;

for i=1:m
    i
    for j=1:n
        if mask(i,j)==1
            % 标注生长季数据
             % 非生长季 设定为NaN
            dg= repmat(reshape(GrowSeason_gobal05(i,j,:),[],1),length(years),1);

            det1=reshape(ET_Ecology05_8214(i,j,:),[],1);
            det2=reshape(ET_Global8216_GLEAM05(i,j,1:end-24),[],1);
            det3=reshape(ET_BESS05_8214(i,j,:),[],1);
%             det4=reshape(ET_PMLv2_8214(i,j,:),[],1);
%             det4=reshape(ET_REA05_8214(i,j,:),[],1);


            
%             dsm=reshape(SM100_Global8214_ERA5Land(i,j,:),[],1);
            dsm=reshape(SM100_Global8214_ERA5Land(i,j,:),[],1);
            dvpd=reshape(VPD_Global8214_CRU05(i,j,:),[],1);

            Z_Season_et(i,j,:,1)= Zscore_month(det1).*dg;
            Z_Season_et(i,j,:,2)= Zscore_month(det2).*dg;
            Z_Season_et(i,j,:,3)= Zscore_month(det3).*dg;
%             Z_Season_et(i,j,:,4)= Zscore_month(det4).*dg;
            
            Z_Season_sm(i,j,:)= Zscore_month(dsm).*dg;
            Z_Season_vpd(i,j,:)= Zscore_month(dvpd).*dg;
%             Z_Season_sm(i,j,:)= dsm.*dg;
%             Z_Season_vpd(i,j,:)= dvpd.*dg;
            
        end
    end
end
%%
AA=Z_Season_et(:,:,6,2);
% BB=std(AA,0,3,'omitnan');

pcolor(lon,lat,AA');shading flat
%% 7x7 滑动
% MovWin=1;
Cor_ET_SM_8214=zeros(m,n,3).*nan;
Cor_ET_VPD_8214=zeros(m,n,3).*nan;
Cor_VPD_SM_8214=zeros(m,n).*nan;
sig_ET_SM_8214=zeros(m,n,3).*nan;
sig_ET_VPD_8214=zeros(m,n,3).*nan;
sig_VPD_SM_8214=zeros(m,n).*nan;

for i=1:m
    i
    for j=1:n
        if mask(i,j)==1
            data_sm=Z_Season_sm(i,j,:);
            data_vpd=Z_Season_vpd(i,j,:);
            
            sample_sm=reshape(data_sm,[],1);
            sample_vpd=reshape(data_vpd,[],1);

            [R_vs,p_vs]=corr(sample_vpd,sample_sm,'Rows','complete');%partialcorr(S_et,S_vpd,S_sm);%
            Cor_VPD_SM_8214(i,j)=R_vs;
            sig_VPD_SM_8214(i,j)=double(p_vs<=0.05);
                   
               
            for k=1:3
                data_et=Z_Season_et(i,j,:,k);
                % 空间窗口的所有样本
                sample_et=reshape(data_et,[],1);
            
            % 删去空值
                if isempty(sample_et)==0            
                    [R_es,p_es]=corr(sample_sm,sample_et,'Rows','complete');%partialcorr(S_et,S_sm,S_vpd);%
                    Cor_ET_SM_8214(i,j,k)=R_es;
                    sig_ET_SM_8214(i,j,k)=double(p_es<=0.01);

                    [R_ev,p_ev]=corr(sample_vpd,sample_et,'Rows','complete');%partialcorr(S_et,S_vpd,S_sm);%
                    Cor_ET_VPD_8214(i,j,k)=R_ev;
                    sig_ET_VPD_8214(i,j,k)=double(p_ev<=0.01);


                end
            end
        end
    end
end
%% Fig1 abc correlation
Cor_ET_SM_Mean=nanmean(Cor_ET_SM_8214,3);
sig_ET_SM_Mean=mean(sig_ET_SM_8214,3); %所有产品显著相关才能认为显著

Cor_ET_VPD_Mean=nanmean(Cor_ET_VPD_8214,3);
sig_ET_VPD_Mean=mean(sig_ET_VPD_8214,3);

Cor_VPD_SM_Mean=Cor_VPD_SM_8214;
sig_VPD_SM_Mean=sig_VPD_SM_8214;
% load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\mask_ColdRegion.mat'

% Cor_ET_SM_8214(mask_ColdRegion==1)=nan;
% sig_ET_SM_8214(mask_ColdRegion==1)=nan;
% 
% Cor_ET_VPD_8214(mask_ColdRegion==1)=nan;
% sig_ET_VPD_8214(mask_ColdRegion==1)=nan;
% 
% Cor_VPD_SM_8214(mask_ColdRegion==1)=nan;
% sig_VPD_SM_8214(mask_ColdRegion==1)=nan;

Cor_ET_SM_Mean(:,lat>=60)=nan;
sig_ET_SM_Mean(:,lat>=60)=nan;

Cor_ET_VPD_Mean(:,lat>=60)=nan;
sig_ET_VPD_Mean(:,lat>=60)=nan;

Cor_VPD_SM_Mean(:,lat>=60)=nan;
sig_VPD_SM_Mean(:,lat>=60)=nan;

NA1={'(a)','(b)','(c)'};
figure(3);h3=tight_subplot(2,1,[.1 .1],[.05 .08],[.05 .05]);
set(gcf,'position',[80 100 650 500]) 
hold on
data_all=cat(3,Cor_ET_SM_Mean,Cor_ET_VPD_Mean,Cor_VPD_SM_Mean);
data_sig=cat(3,sig_ET_SM_Mean,sig_ET_VPD_Mean,sig_VPD_SM_Mean);
rgb1=flipud((othercolor('PiYG10',10))); 
% rgb1=flipud((othercolor('RdYlGn10',10))); 

for k=1:2
    axes(h3(k)); 

    data=data_all(:,:,k);
    sig=data_sig(:,:,k);

     m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-59.75 60]);
    m_pcolor(lon,lat,data'); 
%     shading interp
    m_coast('linewidth',0.5,'color','k');  %% 娴峰哺绾?
    m_grid('linestyle','none','tickdir','on','linewi',1,'color','k','fontsize',12);  %%  澶栧洿杈圭嚎

     caxis(gca,[-1 1]);
     colormap(gca,rgb1);
     if k==2
        m_contfbar(0.97,[.6 1.6],[-1:0.2:1],[-1:0.2:1],'axfrac',.02,'fontsize',12);
     end
 %显著性打点
     for i1=1:m
         for j1=1:n
              if sig(i1,j1)>=0.6
                 x=lon(i1,1);
                 y=lat(j1,1);
                 hold on;
                 if mod(i1,4)==0 & mod(j1,4)==0
                     m_plot(x,y,'k.','markersize',2.5)
                        %m_plot([x x+0.25],[y y+0.25],'k-')
                 end
               end
         end
    end


    if k==1 
        title('r (SM,ET)','fontsize',18);
    elseif k==2 
        title('r (VPD,ET)','fontsize',18);
    else
        title('r (SM,VPD)','fontsize',18);
    end

    hold on
%     m_plot([-180 180],[0 0 ],'--','color',[130 130 130]/255)
%     m_text(-178,8,'JJA','vertical','middle','color','k','fontsize',10) ;
%     m_text(-178,-8,'DJF','vertical','middle','color','k','fontsize',10) ;
    m_text(-180,72,NA1(k),'fontsize',22);

end
print(gcf, '-dpng', '-r600','.\Figure4\FigN1ab.png');
%% density plot
load cmocean_thermal.csv
dataX=Cor_ET_SM_Mean(:);
dataY=Cor_ET_VPD_Mean(:);

ind_NAN=find(isnan(dataX) | isnan(dataY));
dataX(ind_NAN)=[];
dataY(ind_NAN)=[];


figure(1)
set(gcf,'position',[200 200 500 500]) 
axis equal

p=polyfit(dataX,dataY,1);
yfit=polyval(p,dataX);%求拟合后的y值;
% plot(dataX1,dataY1,'b.',dataX1,yfit,'r-','Markersize',3.5);%画图;这是散点的大小颜色样式

rgbD=cmocean_thermal/255;%othercolor('BuPu6',60);
CData=density2C(dataX,dataY,-1:0.1:1,-1:0.1:1,rgbD);
scatter(dataX,dataY,5,CData,'filled');
colormap;%(rgbD)

hold on
%趋势线
plot(dataX,yfit,'k--','LineWidth',0.5);hold on
mdl = fitlm(dataX,dataY);
R2=mdl.Rsquared.Adjusted;
% hold on; plot(mdl)
axis equal
xlim([-1 1]);xticks([-1:0.2:1])
ylim([-1 1]);yticks([-1:0.2:1])
set(gca,'FontSize',12)
text(0.6,0.8,strcat('R^{2} = ',num2str(R2,'%.2f')),'horiz','center','FontSize',18);

% 坐标轴加粗
ax=gca;
ax.LineWidth=1.2;
ax.FontSize=12;
set(gca,'TickDir','in')
text(-1.2,1.2,'(c)','horiz','center','FontSize',22);

ylabel('r(VPD, ET)','FontSize',18);
xlabel('r(SM, ET)','FontSize',18);

box on   
% axis equal


print(gcf, '-dpng', '-r600','.\Figure4\FigN1d.png');
%% 二维密度图
function [CData,h,XMesh,YMesh,ZMesh,colorList]=density2C(X,Y,XList,YList,colorList)
[XMesh,YMesh]=meshgrid(XList,YList);
XYi=[XMesh(:) YMesh(:)];
F=ksdensity([X,Y],XYi);
ZMesh=zeros(size(XMesh));
ZMesh(1:length(F))=F;

h=interp2(XMesh,YMesh,ZMesh,X,Y);
if nargin<5
colorList=[0.2700         0    0.3300
    0.2700    0.2300    0.5100
    0.1900    0.4100    0.5600
    0.1200    0.5600    0.5500
    0.2100    0.7200    0.4700
    0.5600    0.8400    0.2700
    0.9900    0.9100    0.1300];
end
colorFunc=colorFuncFactory(colorList);
CData=colorFunc((h-min(h))./(max(h)-min(h)));
colorList=colorFunc(linspace(0,1,100)');

function colorFunc=colorFuncFactory(colorList)
x=(0:size(colorList,1)-1)./(size(colorList,1)-1);
y1=colorList(:,1);y2=colorList(:,2);y3=colorList(:,3);
colorFunc=@(X)[interp1(x,y1,X,'pchip'),interp1(x,y2,X,'pchip'),interp1(x,y3,X,'pchip')];
end
end

%% 逐月计算Zscore的函数
function y = Zscore_month(x)

y=zeros(length(x),1);

for i=1:12
    x1=x(i:12:end);
    x2=zscore(x1);
    
    y(i:12:end,1)=x2;
end

end
