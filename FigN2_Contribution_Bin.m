clc;clear all
%% Figure1 abc
lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';

% load Contri_SM_8214_Products.mat 
% load Contri_VPD_8214_Products.mat 
load Contri_SM_8214_Products.mat 
load Contri_VPD_8214_Products.mat 

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\LAI_GIMMS05_8214.mat'

m=length(lon);n=length(lat);
Data_Dominance=zeros(m,n).*nan;
% 空间上 所有产品取均值
Contri_SM_8214_mean=nanmean(Contri_SM_8214(:,:,1:3),3);
Contri_VPD_8214_mean=nanmean(Contri_VPD_8214(:,:,1:3),3);
% 根据主导因素分区
Data_Dominance(Contri_SM_8214_mean>Contri_VPD_8214_mean)=1; %1:SM主导；2：VPD主导
Data_Dominance(Contri_SM_8214_mean<=Contri_VPD_8214_mean)=2;
% 只保留北纬60以下
Contri_SM_8214_mean(:,lat>=60)=nan;
Contri_VPD_8214_mean(:,lat>=60)=nan;
Data_Dominance(:,lat>=60)=nan;

% Fig1ab  右侧附上 纬度平均
Contri_VPD_Lat=zeros(length(lat),1);
Contri_SM_Lat=zeros(length(lat),1);

Contri_VPD_Lat_upper=zeros(length(lat),1);
Contri_SM_Lat_upper=zeros(length(lat),1);

Contri_VPD_Lat_lower=zeros(length(lat),1);
Contri_SM_Lat_lower=zeros(length(lat),1);

for i=1:length(lat)
    Contri_VPD_Lat(i,1)=nanmean(Contri_VPD_8214_mean(:,i));
    Contri_VPD_Lat_upper(i,1)=nanmean(Contri_VPD_8214_mean(:,i))+nanstd(Contri_VPD_8214_mean(:,i));%prctile(Contri_VPD_8214_mean(:,i),75,1);
    Contri_VPD_Lat_lower(i,1)=nanmean(Contri_VPD_8214_mean(:,i))-nanstd(Contri_VPD_8214_mean(:,i));%prctile(Contri_VPD_8214_mean(:,i),25,1);
    
    Contri_SM_Lat(i,1)=nanmean(Contri_SM_8214_mean(:,i));
    Contri_SM_Lat_upper(i,1)=nanmean(Contri_SM_8214_mean(:,i))+nanstd(Contri_SM_8214_mean(:,i));%prctile(Contri_SM_8214_mean(:,i),75,1);
    Contri_SM_Lat_lower(i,1)=nanmean(Contri_SM_8214_mean(:,i))-nanstd(Contri_SM_8214_mean(:,i));%prctile(Contri_SM_8214_mean(:,i),25,1);
    
end
%去除空值的行
ind_NAN=find(isnan(Contri_VPD_Lat) | isnan(Contri_SM_Lat));
latN=lat;
latN(ind_NAN)=[];
Contri_VPD_Lat(ind_NAN)=[];
Contri_VPD_Lat_upper(ind_NAN)=[];
Contri_VPD_Lat_lower(ind_NAN)=[];
Contri_SM_Lat(ind_NAN)=[];
Contri_SM_Lat_upper(ind_NAN)=[];
Contri_SM_Lat_lower(ind_NAN)=[];

%  计算面积比例
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\AI_Global8214_CRU05.mat'
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\Weight_Grid.mat'
mask_Dom=Data_Dominance.*AI_Global8214_CRU05.*Weight_Grid;
ind=find(isnan(mask_Dom)==1);
%
dataAI=AI_Global8214_CRU05(:);
dataAI(ind)=[]; %删去空值

dataDom=Data_Dominance(:);
dataDom(ind)=[]; %删去空值

dataGrid=Weight_Grid(:);
dataGrid(ind)=[]; %删去空值
%
%格点面积加权
Dom_SM_Global= sum((dataDom==1).*dataGrid)/sum((isnan(dataDom)==0).*dataGrid)*100;
Dom_VPD_Global=sum((dataDom==2).*dataGrid)/sum((isnan(dataDom)==0).*dataGrid)*100;

% ×面积占比
Dom_SM_Dryland= sum((dataDom==1 & dataAI<=0.65).*dataGrid)/sum((dataAI<=0.65).*dataGrid)*100;
Dom_VPD_Dryland=sum((dataDom==2 & dataAI<=0.65).*dataGrid)/sum((dataAI<=0.65).*dataGrid)*100;

Dom_SM_NonDryland= sum((dataDom==1 & dataAI>0.65).*dataGrid)/sum((dataAI>0.65).*dataGrid)*100;
Dom_VPD_NonDryland=sum((dataDom==2 & dataAI>0.65).*dataGrid)/sum((dataAI>0.65).*dataGrid)*100;

%%
NA1={'(a)','(b)','(e)'};
figure(1);h3=tight_subplot(3,1,[.05 .05],[.05 .06],[.15 .2]);
set(gcf,'position',[80 100 800 800]) 
hold on
data_all=cat(3,Contri_SM_8214_mean,Contri_VPD_8214_mean,Data_Dominance);
% data_all=cat(3,Contri_SM_8214(:,:,3),Contri_VPD_8214(:,:,3),Data_Dominance);


rgb0=flipud((othercolor('RdBu9',18))); 
rgb01=flipud((othercolor('RdBu11',10))); 
rgb1=[rgb01(1:4,:);rgb0(9:10,:);rgb01(7:end,:)]; 

rgb2=[250 134 000;
    033 158 188]/255;

for k=1:3
    axes(h3(k)); 

    data=data_all(:,:,k);
    % 投影 
    m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-60 60]);
    %m_proj('robinson','lon',[-179.75 179.75],'lat',[-59.75 89.75]);
    m_pcolor(lon,lat,data');%shading interp
    m_coast('linewidth',0.8,'color','k');  %% 娴峰哺绾?
    m_grid('box','off','linestyle','none','tickdir','on','linewi',1,'color','k','xtick',[],'ytick',[],'fontsize',10);  %%  澶栧洿杈圭嚎
    % bar设置
    if k<=2
       caxis(gca,[-2 2]);
       colormap(gca,rgb1);
       if k==2

%         ch = colorbar('horiz');

       end
    else
       caxis(gca,[0.5 2.5]);
       colormap(gca,rgb2);
        h1=plot(-9999, 0, 's','MarkerSize',30, 'Color',rgb2(1,:), 'MarkerFaceColor', rgb2(1, :));hold on
        h2=plot(-9999, 0, 's','MarkerSize',30, 'Color',rgb2(2,:), 'MarkerFaceColor', rgb2(2, :));hold on
        lg3=legend([h1,h2],{'Supply-limited','Demand-limited'},'location','Best','FontSize',14)
        set(lg3,'Box','off');
    end

    if k==1 
        title('\DeltaET(\it{SM}|\it{VPD})','fontsize',18);
    elseif k==2       
        title('\DeltaET(\it{VPD}|\it{SM})','fontsize',18);
    else
%         title('Dominant ','fontsize',18);
    end

    hold on
    m_text(-180,72,NA1(k),'fontsize',22);
    
    % 叠加图层
    pos=[0.80 0.683 0.08 0.25;
        0.80 0.37 0.08 0.25;
        0.072 0.09 0.16 0.15];
    h4= axes('Position',pos(k,:));
    
    if k<=2
    
        p0=plot([0 0],[-60 80],'color',[130 130 130]/255,'LineWidth',1.5); %辅助线
        hold on

        % shading 25% 75%
        if k==1
            Contri_Lat_mean=Contri_SM_Lat;
            Contri_Lat_lower=Contri_SM_Lat_upper;
            Contri_Lat_upper=Contri_SM_Lat_lower;
        else
            Contri_Lat_mean=Contri_VPD_Lat;
            Contri_Lat_lower=Contri_VPD_Lat_upper;
            Contri_Lat_upper=Contri_VPD_Lat_lower;
        end

        p1=plot(Contri_Lat_mean,latN,'color',rgb1(2,:),'LineWidth',1.5);
        
        y_vector = [latN', fliplr(latN')];
        pa1=fill( [Contri_Lat_lower',fliplr(Contri_Lat_upper')], y_vector,rgb1(3,:));
        set(pa1, 'edgecolor', 'none');
        set(pa1, 'FaceAlpha', 0.3);

        % y轴标注
        ax = gca;
     %    ax.XAxisLocation = 'origin';
        ax.YAxisLocation = 'origin';

        set(gca,'FontSize',9)

        xlim([-1.5 1.8])
        ylim([-60 60])
        yticks([-60:20:60])
        yticklabels({'60\circS','40\circS','20\circS','0\circ','20\circN','40\circN','60\circN'})
    else % 第三幅图的附图
        dataD=[[Dom_SM_Global,Dom_VPD_Global];[Dom_SM_Dryland,Dom_VPD_Dryland];[Dom_SM_NonDryland,Dom_VPD_NonDryland]];
        b=bar(dataD,0.8);
        b(1).FaceColor = [250 134 000]/255;
        b(2).FaceColor = [033 158 188]/255;

        ax=gca;
        
%         set(gca,'yaxislocation','right');
        
        % 横坐标 斜体
    ax=gca;
    set(gca,'xlim',[0.6 3.4],'xtick',[1:1:3]);
    ax.XTickLabel={'Global','Drylands','Humid lands'};
    ax.LineWidth=1;
    ax.FontSize=12;
    
    ylim([0 100])
   ylabel('Area fraction (%)' ,'FontSize',14)
        
    xtb = get(gca,'XTickLabel');
    xt = get(gca,'XTick');   % 获取横坐标轴刻度句柄
    yt = get(gca,'YTick');    % 获取纵坐标轴刻度句柄
    ytextp=yt(1)*ones(1,length(xt)); 
    xtextp=xt; 
    text(xtextp,ytextp-0.28,xtb,'HorizontalAlignment','right','rotation',40,'fontsize',12); 
    set(gca,'xticklabel','');% 将原有的标签隐去

        set(gca,'TickDir','out')
    end
 
    box off
    set(h4,'color','none') %透明

end

% print(gcf, '-dpng', '-r600','.\Figure4\Fig2abd.png');
%% colorbar
tick=-2:0.4:2;
mode='h';
set(gcf,'position',[80 50 700 200]) 
colorbarn(tick,rgb1,mode);
print(gcf, '-dpng', '-r600','.\Figure1\Fig2bar.png');

%% Fig 2
% 散点密度图
LAI_mean=mean(LAI_GIMMS05_8214,3);
mask_corr=Contri_VPD_8214_mean.*Contri_SM_8214_mean.*AI_Global8214_CRU05;
ind=find(isnan(mask_corr(:))==1);
%
dataAI=reshape(AI_Global8214_CRU05,[],1);
dataAI(ind)=[]; %删去空值

dataCorSM=reshape(Contri_SM_8214_mean,[],1);
dataCorSM(ind)=[]; %删去空值

dataCorVPD=reshape(Contri_VPD_8214_mean,[],1);
dataCorVPD(ind)=[]; %删去空值

%
% figure(4);h3=tight_subplot(2,1,[.06 .1],[.1 .05],[.1 .1]);
% set(gcf,'position',[80 100 350 580]) 
figure(2);h3=tight_subplot(3,1,[.05 .05],[.05 .06],[.15 .2]);
set(gcf,'position',[80 100 800 800]) 

% 散点图
% scatter(dataX,dataY,20,'b','.');hold on;
strm1={'\DeltaET(\it{SM}|\it{VPD})','\DeltaET(\it{VPD}|\it{SM})'};
strm2={'(c)','(d)'};

cm=(othercolor('RdPu9',60)); 

for k=1:2
    
    axes(h3(k));    
    %
    dataX=dataAI;
%     dataY=dataLAI;

    if k==1
        dataY=dataCorSM;
%         dataZ=dataCorSM;
    else
        dataY=dataCorVPD;
%         dataZ=dataCorVPD;
    end
    
    % bin centers (integers), let the number of xbins equals ybins
    xbins = [0:0.1:ceil(max(dataX))];%linspace(floor(min(dataX)), ceil(max(dataX)), 50);%
    ybins = [floor(min(dataY)):0.2:ceil(max(dataY))];%linspace(floor(min(dataY)), ceil(max(dataY)), 50);%

    xNumBins = numel(xbins); 
    yNumBins = numel(ybins);

    % map X/Y values to bin indices
    Xi = floor( interp1(xbins, 1:xNumBins, dataX, 'linear', 'extrap') );
    Yi = floor( interp1(ybins, 1:yNumBins, dataY, 'linear', 'extrap') );

    % count number of elements in each bin
    H = accumarray([Yi(:) Xi(:)], 1, [yNumBins xNumBins]);
    DS=H./sum(H(:))*100;

    % plot 2D histogram
    pcolor(xbins, ybins, DS); shading flat
    cm(1, :) = [1 1 1];  % set the 0 value as blank

   
    axis square
    % colorbar设置
     colormap(gca,cm)
     if k==1
        ch = colorbar('horiz','Position',[0.45 0.75 0.08 0.010],'FontSize',12,'Ticks',[])
        ch.LineWidth=1;
        ch.AxisLocation='in';
        ch.TickDirection='out';
        ch.Title.String='Increasing Density \rightarrow';
        box on
     end
%     caxis([0.1 1])

    hold on
    % 标记AI=0.65 旱地分界线
    plot([0.65 0.65],[-10 10],'--','color',[156 156 156]/255,'LineWidth',1.2); hold on
    text(0.012,-1.8,'Drylands','color','k','FontSize',13);    
    text(1,-1.8,'Humid lands','color','k','FontSize',13);
    
    text(-0.25,2.3,strm2(k),'horiz','center','FontSize',22);
    
    
 % 坐标轴加粗
    ax=gca;
    ax.LineWidth=1.2;
    ax.FontSize=12;
    set(gca,'TickDir','out')
    box on

    xlim([0 2])
    xticks([0:0.5:2])
    ylim([-2 2])
    yticks([-2:1:2])
    set(gca,'yaxislocation','right');
    ylabel(strm1(k),'FontSize',18);
    if k==2
        xlabel('Aridity index','FontSize',14);
    end

    %
end


% print(gcf, '-dpng', '-r300','.\Figure1\Fig2ab.png');
print(gcf, '-dpng', '-r300','.\Figure4\Fig2cd.png');



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






