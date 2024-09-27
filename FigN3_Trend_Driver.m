clc;clear all

%%
lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';
m=length(lon);n=length(lat);
years=[1982:1:2014]';
load Contri_VPD_8214_Products_SMmerge.mat
load Contri_SM_8214_Products_SMmerge.mat
% load Contri_VPD_8214_Products_SMmerge.mat
% load Contri_SM_8214_Products_SMmerge.mat
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\Weight_Grid.mat'

Contri_SM_8214_mean=nanmean(Contri_SM_8214(:,:,1:3),3);
Contri_VPD_8214_mean=nanmean(Contri_VPD_8214(:,:,1:3),3);

Data_Dominance=zeros(m,n).*nan;
Data_Dominance(Contri_SM_8214_mean>Contri_VPD_8214_mean)=1; %1:SM主导；2：VPD主导
Data_Dominance(Contri_SM_8214_mean<=Contri_VPD_8214_mean)=2;


Data_Dominance(:,lat>=60)=nan;
Weight_Grid(:,lat>=60)=nan;

%% Fig3 ab
load Contri_SM_8214Trend_OBS_BefMean_SMmerge.mat
load Contri_VPD_8214Trend_OBS_BefMean_SMmerge.mat
% load Contri_SM_8214Trend_OBS_BefMean_SMmerge.mat
% load Contri_VPD_8214Trend_OBS_BefMean_SMmerge.mat

rgb1=flipud((othercolor('BuDOr_18',11))); 
rgb2=flipud(rgb1([3,5,6,10],:));

% 高纬度
%不计算cold region
Contri_SM_8214Trend_OBS_BefMean(:,lat>=60)=nan;
Contri_VPD_8214Trend_OBS_BefMean(:,lat>=60)=nan;
Contri_SM_8214Sig_OBS_BefMean(:,lat>=60)=nan;
Contri_VPD_8214Sig_OBS_BefMean(:,lat>=60)=nan;

data_all=cat(3,Contri_SM_8214Trend_OBS_BefMean.*double(Data_Dominance==1),Contri_VPD_8214Trend_OBS_BefMean.*double(Data_Dominance==2));
sig_all=cat(3,Contri_SM_8214Sig_OBS_BefMean.*double(Data_Dominance==1),Contri_VPD_8214Sig_OBS_BefMean.*double(Data_Dominance==2));
data_all(data_all==0)=-9999; 
sig_all(sig_all==0)=-9999;

leg1={'Dec*','Dec','Inc','Inc*'};

rgb01=(othercolor('BuDOr_12',20));
rgb02=(othercolor('BuDOr_12',14));
rgb3=[rgb02(1:6,:);rgb01(10:11,:);rgb02(9:14,:)];
rgb4=[rgb01(7,:);rgb02(2,:);rgb02(13,:);rgb01(12,:)];
NA1={'(a)','(b)'};

figure(4);h3=tight_subplot(1,2,[.05 .05],[.05 .06],[.05 .1]);
set(gcf,'position',[80 100 1800 400]) 
hold on


for k=1:2
    axes(h3(k)); 
 % 图层1
    data=data_all(:,:,k);   
    sig=sig_all(:,:,k);  
% 非限制区 画灰色
    hold on
    data_mask=double(data==-9999);
    data_mask(data_mask==0)=nan;
    m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-60 60]);

    m_pcolor(lon,lat,data_mask');   
    colormap([0.9 0.9 0.9]);freezeColors
    caxis([0.5 1.5])  
    m_coast('linewidth',1,'color','k');  %% 海岸线
  
    % 趋势上色
    data(data==-9999)=nan;% 非限制区
    sig(sig==-9999)=nan;% 非限制区
    m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-60 60]);   
    m_grid('box','off','linestyle','none','tickdir','on','linewi',1,'color','k','xtick',[],'ytick',[],'fontsize',10);
    
    m_pcolor(lon,lat,data'*10);%shading interp % /10a
    caxis(gca,[-0.7 0.7]);
    m_coast('linewidth',0.8,'color','k');  %% 海岸线
    colormap(gca,rgb3);freezeColors
    
     % 显著性打点
    for i1=1:m
        for j1=1:n
            if sig(i1,j1)<=0.1
               x=lon(i1,1);
               y=lat(j1,1);
               hold on;
               if mod(i1,4)==0 & mod(j1,4)==0
                  m_plot(x,y,'k.','markersize',3)
               end
             end
        end
    end
     
    if k==1 
        title('Trends in \DeltaET(\it{SM}|\it{VPD})','fontsize',24);
    else 
        title('Trends in \DeltaET(\it{VPD}|\it{SM})','fontsize',24);
    end

    hold on
    m_text(-190,68,NA1(k),'fontsize',24);
    
  %colorbar
    ch2 = colorbar('v','Position',[0.91 0.32 0.012 0.4],'FontSize',14)
    ch2.LineWidth=1;
    ch2.AxisLocation='out';
    ch2.TickDirection='out';
    ch2.YTick=[-0.6:0.3:0.6];

    % 叠加图层
    pos=[0.07 0.21 0.08 0.31;
        0.52 0.21 0.08 0.31;];
    h4= axes('Position',pos(k,:));
%    set(h4,'color','none') %透明

    data1=zeros(m,n).*nan;
    data1(data>0 & sig>0.1)=4; % 不显著增加
    data1(data>0 & sig<=0.1)=3; % 显著增加
    data1(data<=0 & sig<=0.1)=2;% 显著减少
    data1(data<=0 & sig>0.1)=1;% 不显著减小
    % 非限制区
    for ii=1:4
        
         db=nansum((data1(:)==ii).*Weight_Grid(:))/nansum((isnan(data1(:))==0).*Weight_Grid(:))*100;
         
         b(ii)=bar(ii,db,0.8);%
         b(ii).LineWidth = 0.5;
        hold on
        b(ii).EdgeColor = rgb4(ii,:);b(ii).FaceColor = rgb4(ii,:);
        
        % 标记数字    
        if k==1 & ii==3
            b(ii).EdgeColor = [0 0 0];
            b(ii).LineWidth = 3;
            text(ii-0.38,db+8,strcat(num2str(db,2),'%'),'fontsize',18,'Color',[0 0 0]);
        end
    end

    % 横坐标 斜体
    ax=gca;
    set(gca,'xlim',[0.4 4.6],'xtick',[1:1:4]);
    ax.XTickLabel={'Dec','Dec*','Inc*','Inc'};
    ax.LineWidth=1.2;
    ax.FontSize=14;
        
    xtb = get(gca,'XTickLabel');
    xt = get(gca,'XTick');   % 获取横坐标轴刻度句柄
    yt = get(gca,'YTick');    % 获取纵坐标轴刻度句柄
    ytextp=yt(1)*ones(1,length(xt)); 
    xtextp=xt; 
    text(xtextp,ytextp-2,xtb,'HorizontalAlignment','right','rotation',40,'fontsize',18); 
    set(gca,'xticklabel','');% 将原有的标签隐去
    
    
    ylim([0 50])
    yticks([0:25:50]);
    ylabel('Area fraction (%)','fontsize',18)
    
    
   box off
   set(h4,'color','none') %透明
 end
     

print(gcf, '-dpng', '-r600','.\FigureR\FigR3ab_Mean_Trend.png');
%% Bar
tick=-0.7:0.1:0.7;
mode='v';
figure(2);
set(gcf,'position',[80 50 800 300]) 
label={'-0.6',' ',' ', '-0.3',' ',' ', '0',' ',' ', '0.3',' ',' ', '0.6'};

colorbarn(tick,rgb3,mode,label);
print(gcf, '-dpng', '-r600','.\Figure1\FigNN3bar.png');
%% Attribution 
load Attri_SM_pMK.mat; %AttriSMpMK4;
load Attri_VPD_pMK.mat; 


Attribution_SM=Attri_SM_pMK;
Attribution_SM(Data_Dominance==2)=-9999;
Attribution_VPD=Attri_VPD_pMK;
Attribution_VPD(Data_Dominance==1)=-9999;

Attribution_SM(:,lat>60)=nan;
Attribution_VPD(:,lat>60)=nan;
NA2={'(c)','(d)'};
hold on

rgb5=[150 150 150;
%     42 145 199;
    254 183 005;
    251 132 002;
    42,152,69;
    24 116 205]/255;

leg1={'Not sig','SM','VPD','LAI','Other drivers'};
% leg1={'Not sig','P','T','SM','VPD','LAI'};
data_all=cat(3,Attribution_SM,Attribution_VPD);

figure(4);h3=tight_subplot(1,2,[.05 .05],[.05 .06],[.05 .1]);
set(gcf,'position',[80 100 1800 400]) 
Driver_All=zeros(4,1).*nan; %统计显著变化区域上各因子主导的比例
Driver_Inc=zeros(4,1).*nan;
Driver_Dec=zeros(4,1).*nan;

for k=1:2
    axes(h3(k)); 
 % 图层1
    data=data_all(:,:,k);   
% 非限制区 画灰色
    hold on
    data_mask=double(data==-9999);
    data_mask(data_mask==0)=nan;
    m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-60 60]); 
    m_pcolor(lon,lat,data_mask');      
    colormap([0.9 0.9 0.9]);freezeColors
    caxis([0.5 1.5])  
    m_coast('linewidth',0.8,'color','k');  %% 海岸线
  
    % 主导因素上色
    data(data==-9999)=nan;% 非限制区
    m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-60 60]); 
    m_grid('box','off','linestyle','none','tickdir','on','linewi',1,'color','k','xtick',[],'ytick',[],'fontsize',10);
    
    m_pcolor(lon,lat,data');%shading interp
    caxis(gca,[-0.5 4.5]); colorbar
    m_coast('linewidth',0.8,'color','k');  %% 海岸线
   
     
    if k==1 
        title('Dominant drivers for trends in \DeltaET(\it{SM}|\it{VPD})','fontsize',24);
    else 
        title('Dominant drivers for trends in \DeltaET(\it{VPD}|\it{SM})','fontsize',24);
    end

    hold on
    m_text(-190,68,NA2(k),'fontsize',24);
    
     %bar
     colormap(gca,rgb5);freezeColors
   
    ch2 = colorbar('v','Position',[0.905 0.3 0.012 0.4],'FontSize',12,'Ticks',[0:1:4],'Ticklabels',leg1,'fontsize',18)
    ch2.LineWidth=1;
    ch2.AxisLocation='out';
    ch2.TickDirection='out';
        

% 叠加图层
pos=[0.06 0.29 0.09 0.21;
        0.51 0.29 0.09 0.21;];
h4= axes('Position',pos(k,:));
%    set(h4,'color','none') %透明

    data(data==-9999)=nan;
    % 考虑四个因子  统计在显著变化的区域
    if k==1
        dSig=double(Contri_SM_8214Sig_OBS_BefMean<0.1);
        dSig(dSig==0)=nan;
        dInc=(Contri_SM_8214Trend_OBS_BefMean>0).*(Contri_SM_8214Sig_OBS_BefMean<0.1);
        dInc(dInc==0)=nan;
        dDec=(Contri_SM_8214Trend_OBS_BefMean<=0).*(Contri_SM_8214Sig_OBS_BefMean<0.1);
        dDec(dDec==0)=nan;
    else
        dSig=double(Contri_VPD_8214Sig_OBS_BefMean<0.1);
        dSig(dSig==0)=nan;
        dInc=(Contri_VPD_8214Trend_OBS_BefMean>0).*(Contri_VPD_8214Sig_OBS_BefMean<0.1);
        dInc(dInc==0)=nan;
        dDec=(Contri_VPD_8214Trend_OBS_BefMean<=0).*(Contri_VPD_8214Sig_OBS_BefMean<0.1);
        dDec(dDec==0)=nan;
        
    end
    

    for ii=1:4
        db0=nansum((data(:)==ii).*dSig(:).*Weight_Grid(:))/nansum((isnan(data(:))==0).*Weight_Grid(:))*100;
        Driver_All(ii,1)=db0;
        
        db1=nansum((data(:)==ii).*dInc(:).*Weight_Grid(:))/nansum((isnan(data(:))==0).*Weight_Grid(:))*100;
        Driver_Inc(ii,1)=db1;
        
        db2=nansum((data(:)==ii).*dDec(:).*Weight_Grid(:))/nansum((isnan(data(:))==0).*Weight_Grid(:))*100;
        Driver_Dec(ii,1)=db2;
    end
    
    % 堆叠图
    
        x=[1:1:2]';
        Drivers=[Driver_Dec,Driver_Inc];
        b=barh(x,Drivers,0.8,'stacked');
        for ii=1:4
             set(b(ii),'FaceColor',rgb5(ii+1,:),'EdgeColor',rgb5(ii+1,:));
        end

    ax=gca;
    ax.LineWidth=1.2;
    ax.FontSize=14;
    set(gca,'ylim',[0.4 2.6],'ytick',[1:1:2],'YTickLabel',{'Dec*','Inc*'},'fontsize',18);
    
    xlim([0 50]);xticks([0:25:50]);
    xlabel('Area fraction (%)','fontsize',18)
    
    
   box off
   set(h4,'color','none') %透明
 end
     

print(gcf, '-dpng', '-r600','.\FigureR\FigNN3cd_Mean_Trend.png');
%% 相对贡献
load RelContr_SM_pMK.mat; 
load RelContr_VPD_pMK.mat; 

RelContr_SM=zeros(3,1).*nan;
RelContr_VPD=zeros(3,1).*nan;

for k=1:3
    % SM Driver
    d1=RelContr_SM_pMK(:,:,k);
    d1(Data_Dominance==2)=nan;
    d1(:,lat>60)=nan;
    d1M=nansum(d1(:).*Weight_Grid(:))/nansum((isnan(d1(:))==0).*Weight_Grid(:))*100;
%     
    RelContr_SM(k,1)=nanmean(reshape(d1*100,[],1));

    % VPD Driver
    d2=RelContr_VPD_pMK(:,:,k);
    d2(Data_Dominance==1)=nan;
    d2(:,lat>60)=nan;
    d2M=nansum(d2(:).*Weight_Grid(:))/nansum((isnan(d2(:))==0).*Weight_Grid(:))*100;

    RelContr_VPD(k,1)=nanmean(reshape(d2*100,[],1));
    
end

figure(4);h3=tight_subplot(1,2,[.1 .1],[.15 .2],[.1 .1]);
set(gcf,'position',[80 100 800 500]) 
data_all=[RelContr_SM,RelContr_VPD];
for k=1:2
    axes(h3(k)); 
 % 图层1
    data=data_all(:,k);
    
    for ii=1:3
        
         db=data(ii);
 
         b(ii)=bar(ii,db,0.7);%
         b(ii).LineWidth = 0.5;
        hold on
        b(ii).EdgeColor = rgb5(ii+1,:);b(ii).FaceColor = rgb5(ii+1,:);
%         text(ii-0.37,db+5,strcat(num2str(db,2),'%'),'fontsize',16,'Color',[0 0 0]);
    end
     xticks([1:1:3]);xlim([0.4 3.6]);
    ax=gca;
    ax.XTickLabel={'SM','VPD','LAI'};
    ax.FontSize=14;
    
    ylim([0 50])
    ylabel('Relative contribution (%)','fontsize',16)
    
%     text(0,62,NA1(k),'fontsize',20);
        if k==1 
        title(['(a)   Drivers for trends in' sprintf('\n') '\DeltaET(\it{SM}|\it{VPD})'],'fontsize',18);
    else 
        title(['(b)   Drivers for trends in' sprintf('\n') '\DeltaET(\it{VPD}|\it{SM})'],'fontsize',18);
    end
end

print(gcf, '-dpng', '-r300','.\Figure4\FigS4.png');
%%
function [ax1,c] = colorbarn(tick,color,mode,label)
[m,~]=size(color);
len_tick=length(tick);
if len_tick-m ~= 1
    return
end
% figure;

set(gcf,'color','w');
if isequal(mode,'h')
    c=axes('position',[0.1 0.1 0.8 0.03]);
    cdata=zeros(1,m);
    for i=1:m
        cdata(i)=(tick(i)+tick(i+1))/2;
    end
    axis([0 m 0 1])
    line([0 1],[0 0],'linewidth',2,'color','w','parent',c)
    line([m-1 m],[0 0],'linewidth',2,'color','w','parent',c)
    colormap(color);
    %------------------------
    pat_v1=[1 0;1 1;2 1;2 0];pat_f1=reshape(1:(m-2)*4,4,m-2)';
    for j=2:m-2
        pat_v1=[pat_v1;[j 0;j 1;j+1 1;j+1 0]];
    end
    pat_col1=[cdata(2:end-1)]';
    patch('Faces',pat_f1,'Vertices',pat_v1,'FaceVertexCData',pat_col1,'FaceColor','flat');
    %--------------------------------------------------   
    pat_v2=[0 0.5;1 0;1 1;m-1 1;m-1 0;m 0.5];
    pat_f2=[1 2 3;4 5 6];
    pat_col2=[cdata(1);cdata(end)];
    patch('Faces',pat_f2,'Vertices',pat_v2,'FaceVertexCData',pat_col2,'FaceColor','flat');
    %---------------------------------------------------------
    set(c,'color','none','xcolor','k','ycolor','none');
    box off
%     line([1 m-1],[0 0],'color','k')
%     line([1 m-1],[1 1],'color','k')
%     for i=2:m-2
%         line([i i],[0 1],'color','k')
%     end
 %   set(c,'xtick',1:m-1,'xticklabel',num2cell(tick(2:end-1)),'ytick',[])
    set(c,'xtick',1:m-1,'xticklabel',label,'ytick',[])
%     ax1=axes('position',[0.1 0.2 0.8 0.7]);
elseif isequal(mode,'v')
    c=axes('position',[0.87 0.11 0.02 0.81]);
    set(c,'yaxislocation','right');
    cdata=zeros(1,m);
    for i=1:m
        cdata(i)=(tick(i)+tick(i+1))/2;
    end
    axis([0 1 0 m])
    line([1 1],[0 1],'linewidth',2,'color','w','parent',c)
    line([1 1],[m-1 m],'linewidth',2,'color','w','parent',c)
    colormap(color);
    %--------------------------------------------
    pat_v1=[0 1;1 1;1 2;0 2];pat_f1=reshape(1:(m-2)*4,4,m-2)';
    for j=2:m-2
        pat_v1=[pat_v1;[0 j;1 j;1 j+1;0 j+1]];
    end
    pat_col1=[cdata(2:end-1)]';
    patch('Faces',pat_f1,'Vertices',pat_v1,'FaceVertexCData',pat_col1,'FaceColor','flat','lineWidth',1.0);
    %--------------------------------------------------   
    pat_v2=[0.5 0;1 1;0 1;1 m-1;0.5 m;0 m-1];
    pat_f2=[1 2 3;4 5 6];
    pat_col2=[cdata(1);cdata(end)];
    patch('Faces',pat_f2,'Vertices',pat_v2,'FaceVertexCData',pat_col2,'FaceColor','flat','lineWidth',1.0);
    %------------------------------------------------------
    set(c,'color','none','xcolor','none','ycolor','k');
    box off
    set(c,'ytick',1:m-1,'yticklabel',label,'xtick',[],'fontsize',14)
%     ax1=axes('position',[0.11 0.11 0.74 0.81]);
else
    disp('colorbarn格式出错')
    return
end
end


