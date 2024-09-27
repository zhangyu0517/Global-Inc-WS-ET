%%
clc;clear all

%%
Mod_ind=[2,3,4,5,6,7,8,9,10,11,12,13,14,17,18,20,21,22,23,24];
Num_Mod=length(Mod_ind);%21; %模型数量
OBS_ind=[1,2,3];
Num_OBS=length(OBS_ind);%21; %观测数量

lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';
m=length(lon);n=length(lat);

years=[1982:1:2014]';
 addpath('F:\ZhangYu\Global_ET_SMVPD\Code')
load GS_start_gobal05.mat 
load GS_end_gobal05.mat 
load mask_Global05.mat
%% 保留 观测显著变化的区域
load Contri_SM_8214Trend_OBS_BefMean.mat
load Contri_VPD_8214Trend_OBS_BefMean.mat

% Contri_SM_8214Sig_OBS_mean=nanmean(Contri_SM_8214Sig_OBS<=0.1,3);
% Contri_VPD_8214Sig_OBS_mean=nanmean(Contri_VPD_8214Sig_OBS<=0.1,3);
load Contri_VPD_8214_Products.mat
load Contri_SM_8214_Products.mat
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\Weight_Grid.mat'

Contri_SM_8214_mean=nanmean(Contri_SM_8214(:,:,OBS_ind),3);
Contri_VPD_8214_mean=nanmean(Contri_VPD_8214(:,:,OBS_ind),3);

Data_Dominance=zeros(m,n).*nan;
Data_Dominance(Contri_SM_8214_mean>Contri_VPD_8214_mean)=1; %1:SM主导；2：VPD主导
Data_Dominance(Contri_SM_8214_mean<=Contri_VPD_8214_mean)=2;

%% Figure
Exclude_ColdRegion=ones(m,n);
Exclude_ColdRegion(:,lat>60)=nan;
%Mod
load Contri_SM_8214Trend05_CMIP_7x7_All.mat
load Contri_VPD_8214Trend05_CMIP_7x7_All.mat

Contri_SM_8214Trend_MOD=Contri_SM_8214Trend.*Exclude_ColdRegion;
Contri_VPD_8214Trend_MOD=Contri_VPD_8214Trend.*Exclude_ColdRegion;

SM_Dominance=zeros(m,n).*nan;
VPD_Dominance=zeros(m,n).*nan;
SM_Dominance(Contri_SM_8214_mean>Contri_VPD_8214_mean)=1; %1:SM主导；2：VPD主导
VPD_Dominance(Contri_SM_8214_mean<=Contri_VPD_8214_mean)=1;

%% PDF
NA1={'(c)','(d)'};
figure(1);h3=tight_subplot(1,2,[.1 .1],[.2 .2],[.18 .18]);
set(gcf,'position',[80 100 1400 300])

for k=1:2
    
   axes(h3(k)); 
    
    if k==1
        dOBS=Contri_SM_8214Trend_OBS_BefMean.*SM_Dominance*10; %/10a
        
        dMOD=Contri_SM_8214Trend_MOD.*SM_Dominance*10;   %/10a
    else
        dOBS=Contri_VPD_8214Trend_OBS_BefMean.*VPD_Dominance*10;   %/10a
        
        dMOD=Contri_VPD_8214Trend.*VPD_Dominance*10;   %/10a
    end
    % 辅助线
    p0=plot([0 0],[0 5],'-','Color',[0 0 0]/255,'linewidth',1); hold on

     for i2=1:Num_Mod %13个simulation
            da=dMOD(:,:,Mod_ind(i2));
            Dis_sim=da(isnan(da)==0 & isinf(da)==0);
            
             % 绘制PDF
            [f,xi] = ksdensity(Dis_sim);
            p4=plot(xi,f,'Color',[150 150 150]/255,'linewidth',0.8);
          
            hold on
            
        end


        do=dOBS;
        Dis_obs=do(isnan(do)==0 & isinf(do)==0);
            
             % 绘制PDF
        [f,xi] = ksdensity(Dis_obs);
        p1=plot(xi,f,'-','Color',[255 69 0]/255,'linewidth',1.5); hold on


         ax=gca;
        ax.LineWidth=1.0;
        set(ax,'FontSize',14)

       text(-1.1,4.0,NA1(k),'fontsize',20);
       xlim([-0.9 0.9])
       xticks([-0.8:0.2:0.8])
       
       
    if k==1 
        xlabel('Trends in \DeltaET(\it{SM}|\it{VPD})','FontSize',16)       
    else
        xlabel('Trends in \DeltaET(\it{VPD}|\it{SM})','FontSize',16)       
    end
    
   ylabel('PDF','FontSize',16)
        
   
    
    ylim([0 3.5])
    yticks([0:1:3])
    

    lg1=legend([p1],{'OBS'},'location','NorthWest','fontsize',14,'TextColor',[255 69 0]/255);    
    set(lg1,'Orientation','horizon')
    legend boxoff
    ah2 = axes('position',get(gca,'position'),'visible','off');
    lg2=legend(ah2,[p4],{'CMIP6 members'},'location','NorthWest','fontsize',14,'TextColor',[120 120 120]/255);
    legend boxoff
end


print(gcf, '-dpng', '-r600','.\Figure4\Fig4cd.png');
% data_sig=cat(3,nanmean(Contri_SM_8214Sig<=0.05,3),nanmean(Contri_VPD_8214Sig<=0.05,3));
%% 统计 变化方向  先窗口平均再求趋势
% 观测平均
Contri_SM_8214Trend_OBS_mean=Contri_SM_8214Trend_OBS_BefMean.*SM_Dominance;%nanmean(Contri_SM_8214Trend_OBS,3);
Contri_VPD_8214Trend_OBS_mean=Contri_VPD_8214Trend_OBS_BefMean.*VPD_Dominance;%nanmean(Contri_VPD_8214Trend_OBS,3);

% 模拟平均
load Contri_SM_8214MovWin_CMIP_7x7_All.mat
load Contri_VPD_8214MovWin_CMIP_7x7_All.mat
Contri_SM_8214MovWin_MOD_mean=nanmean(Contri_SM_8214MovWin_CMIP(:,:,:,Mod_ind),4);
Contri_VPD_8214MovWin_MOD_mean=nanmean(Contri_VPD_8214MovWin_CMIP(:,:,:,Mod_ind),4);

Contri_SM_8214Trend_MOD_mean=zeros(m,n).*nan;
Contri_VPD_8214Trend_MOD_mean=zeros(m,n).*nan;
Contri_SM_8214Sig_MOD_mean=zeros(m,n).*nan;
Contri_VPD_8214Sig_MOD_mean=zeros(m,n).*nan;
MovWin=7; %5
TimeWin=11; %11年滑动
Num_TemWin=length(years)-TimeWin+1;

for i=1:m    
     i
    for j=1:n
            % 趋势
         dy1=reshape(Contri_SM_8214MovWin_MOD_mean(i,j,:),[],1);
         dy2=reshape(Contri_VPD_8214MovWin_MOD_mean(i,j,:),[],1);
         if isnan(dy1(1))==0 & isnan(dy2(1))==0 %不存在空值
             num=[1:1:Num_TemWin]';
             data1=[num dy1];
             [sig1,sen1]=ktaubZY(data1,0.5);
             Contri_SM_8214Trend_MOD_mean(i,j)=sen1;
             Contri_SM_8214Sig_MOD_mean(i,j)=sig1;

             data2=[num dy2];
             [sig2,sen2]=ktaubZY(data2,0.5);
             Contri_VPD_8214Trend_MOD_mean(i,j)=sen2;
             Contri_VPD_8214Sig_MOD_mean(i,j)=sig2;
         end
      end
end


% 1：Both Increas 2：Both decrease 3：Diff
Compare_SM=zeros(m,n).*nan;
Compare_SM(isnan(Contri_SM_8214Trend_MOD_mean)==0)=1; %Diff
Compare_SM(Contri_SM_8214Trend_MOD_mean<=0 & Contri_SM_8214Trend_OBS_mean<=0)=2;
Compare_SM(Contri_SM_8214Trend_MOD_mean>0 & Contri_SM_8214Trend_OBS_mean>0)=3;
Compare_SM(Data_Dominance==2)=-9999;


Compare_VPD=zeros(m,n).*nan;
Compare_VPD(isnan(Contri_VPD_8214Trend_MOD_mean)==0)=1; %Diff
Compare_VPD(Contri_VPD_8214Trend_MOD_mean<=0 & Contri_VPD_8214Trend_OBS_mean<=0)=2;
Compare_VPD(Contri_VPD_8214Trend_MOD_mean>0 & Contri_VPD_8214Trend_OBS_mean>0)=3;
Compare_VPD(Data_Dominance==1)=-9999;
%%
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN3\Weight_Grid.mat'
% 0 变化不显著
% data_all=cat(3,Compare_SM.*(Contri_SM_8214Sig_OBS_BefMean<=0.05),Compare_VPD.*(Contri_VPD_8214Sig_OBS_BefMean<=0.05));
data_all=cat(3,Compare_SM,Compare_VPD);

NA1={'(e)','(f)'};
figure(1);h3=tight_subplot(1,2,[.05 .05],[.1 .1],[.1 .1]);
set(gcf,'position',[80 100 1200 400]) 

% set(gcf,'position',[80 100 1500 300]) 
hold on

rgb3=[
%     222 222 222;
    150 150 150;
    0 160 135;
    230 75 43;]/255;

for k=1:2
    axes(h3(k)); 
    
    data=data_all(:,:,k); 
    
    % 非限制区 画灰色
    data_mask=double(data==-9999);
    data_mask(data_mask==0)=nan;
    m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-59.75 60]);
    m_pcolor(lon,lat,data_mask');      hold on
    colormap([0.9 0.9 0.9]);freezeColors
    caxis([0.5 1.5])  
    m_coast('linewidth',0.8,'color','k');  %% 海岸线

    % 上色
    data(data==-9999)=nan;
     m_proj('Miller Cylindrical','lon',[-179.75 179.75],'lat',[-59.75 60]);
%m_proj('Miller','lon',[-179.75 179.75],'lat',[-89.75 89.75]);
    m_pcolor(lon,lat,data');% shading interp
    m_coast('linewidth',0.8,'color','k');  %% 娴峰哺绾?
   m_grid('box','off','linestyle','none','tickdir','on','linewi',1,'color','k','xtick',[],'ytick',[],'fontsize',10);

     caxis(gca,[0.5 3.5]);
     colormap(rgb3);freezeColors
%      colorbar

    hold on

    m_text(-195,72,NA1(k),'fontsize',20);
%      colormap(gca,rgb3);
     % 要加权！
%      d0=sum((data(:)==0).*Weight_Grid(:))/sum((isnan(data(:))==0).*Weight_Grid(:))*100;
     d1=nansum((data(:)==1).*Weight_Grid(:))/nansum((isnan(data(:))==0).*Weight_Grid(:))*100;
     d2=nansum((data(:)==2).*Weight_Grid(:))/nansum((isnan(data(:))==0).*Weight_Grid(:))*100;
     d3=nansum((data(:)==3).*Weight_Grid(:))/nansum((isnan(data(:))==0).*Weight_Grid(:))*100;

    if k==2
         f1=plot(-9999, 0, 's','MarkerSize',30, 'Color',rgb3(3,:), 'MarkerFaceColor', rgb3(3, :));hold on
         f2=plot(-9999, 0, 's','MarkerSize',30, 'Color',rgb3(2,:), 'MarkerFaceColor', rgb3(2, :));hold on
         f3=plot(-9999, 0, 's','MarkerSize',30, 'Color',rgb3(1,:), 'MarkerFaceColor', rgb3(1, :));hold on
         lg3=legend([f1,f2,f3],{'+ +';'- -';'Inconsistent'},'Orientation','horizon','location','Best','FontSize',16)
         set(lg3,'Box','off');
    end
 
      % 叠加图层
    pos=[0.100 0.26 0.08 0.16;
        0.525 0.26 0.08 0.16;];
    h4= axes('Position',pos(k,:));
    
    p=pie([d3 d2 d1]);
    p(1).FaceColor = rgb3(3,:);p(2).FontSize = 14;
    p(3).FaceColor = rgb3(2,:);p(4).FontSize = 14;
    p(5).FaceColor = rgb3(1,:);p(6).FontSize = 14;
    
    
    box off
    set(h4,'color','none') %透明
    
    

end
print(gcf, '-dpng', '-r600','.\Figure4\FigN4Diff.png');
%%



