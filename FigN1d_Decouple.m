%%
clc;clear all

%%
lat=[89.75:-0.5:-89.75]';
lon=[-179.75:0.5:179.75]';

years=[1982:1:2014]';
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\SM100_Global8214_ERA5Land.mat'
load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\VPD_Global8214_CRU05.mat'
% 
% load GS_start_gobal05.mat
% load GS_end_gobal05.mat

load 'F:\ZhangYu\Global_ET_SMVPD\CodeN2\mask_Global05.mat'
% load 'F:\ZhangYu\Global_ET_SMVPD\CodeN\mask_ColdRegion.mat'
load GrowSeason_gobal05.mat

mask_Global05(:,lat>60)=nan; %不要寒区
%%
bin=10;
bins=[0:bin:100];
num_bin=length(bins)-1;

m=length(lon);n=length(lat);

S_sm=zeros(length(years),1).*nan;
S_vpd=zeros(length(years),1).*nan;

Z_Season_sm= zeros(m,n,length(years)*12).*nan;
Z_Season_vpd= zeros(m,n,length(years)*12).*nan;
for i=1:m
    i
    for j=1:n
        if mask_Global05(i,j)==1
             % 非生长季 设定为NaN
            dg= repmat(reshape(GrowSeason_gobal05(i,j,:),[],1),length(years),1);
            dsm=reshape(SM100_Global8214_ERA5Land(i,j,:),[],1);
            dvpd=reshape(VPD_Global8214_CRU05(i,j,:),[],1);

       
            
            Z_Season_sm(i,j,:)= Zscore_month(dsm).*dg;
            Z_Season_vpd(i,j,:)= Zscore_month(dvpd).*dg;
        end
    end
end

%% 空间滑动
MovWin=7; %5
Cor_VPD_SM_all=zeros(m,n).*nan;
sig_VPD_SM_all=zeros(m,n).*nan;
Cor_VPD_SM_SMbin=zeros(m,n,num_bin).*nan;
sig_VPD_SM_SMbin=zeros(m,n,num_bin).*nan;
Cor_VPD_SM_VPDbin=zeros(m,n,num_bin).*nan;
sig_VPD_SM_VPDbin=zeros(m,n,num_bin).*nan;

for i=1:m    
     i
    for j=1:n
        if (mask_Global05(i,j)==1 & i>(MovWin-1)/2 & j>(MovWin-1)/2 & m-i>(MovWin-1)/2 & n-j>(MovWin-1)/2)
                data_sm=Z_Season_sm(i-(MovWin-1)/2:i+(MovWin-1)/2,j-(MovWin-1)/2:j+(MovWin-1)/2,:);
                data_vpd=Z_Season_vpd(i-(MovWin-1)/2:i+(MovWin-1)/2,j-(MovWin-1)/2:j+(MovWin-1)/2,:);
                
                % 空间窗口的所有样本
                sample_sm=reshape(data_sm,[],1);
                sample_vpd=reshape(data_vpd,[],1);
                
                num_sample=length(sample_sm);
                % 删去空值 ET为负的值
                indNA=find(isnan(sample_sm)==1 |isnan(sample_vpd)==1);
                sample_sm(indNA)=[];
                sample_vpd(indNA)=[];
                %保证样本量
            if length(sample_sm)>num_sample/2 & length(sample_vpd)>num_sample/2 %& length(unique(sample_sm))>length(sample_vpd)/2%不为空集 不重复值小于50%
                x1=sample_sm;
                x2=sample_vpd;
                
                [R_vs0,p_vs0]=corr(x1,x2);
                Cor_VPD_SM_all(i,j)=R_vs0;
                sig_VPD_SM_all(i,j)=p_vs0;
                
                for k=1:num_bin 
                    % SM bin
                    ind1=find(x1>= prctile(x1,bins(k)) & x1< prctile(x1,bins(k+1)));    
                    if isempty(ind1)==0
                        [R_vs1,p_vs1]=corr(x2(ind1,1),x1(ind1,1));%partialcorr(S_et,S_vpd,S_sm);%
                        Cor_VPD_SM_SMbin(i,j,k)=R_vs1;
                        sig_VPD_SM_SMbin(i,j,k)=p_vs1;                   
                    end


                    %VPD bin
%                     k=10                   
                    ind2=find(x2>= prctile(x2,bins(k)) & x2< prctile(x2,bins(k+1)));
                    if isempty(ind2)==0
                        [R_vs2,p_vs2]=corr(x2(ind2,1),x1(ind2,1));%partialcorr(S_et,S_vpd,S_sm);%
                        Cor_VPD_SM_VPDbin(i,j,k)=R_vs2;
                        sig_VPD_SM_VPDbin(i,j,k)=p_vs2;
%                     R_vs2
                    
                    end              
                end
                
                
            end
        end
    end
end
%% FigS2 a correlation boxplot
Corr_Bin=zeros(length(Cor_VPD_SM_all(:)),num_bin*2 +1);
for p=1:num_bin*2 +1
    if p==1
       Corr_Bin(:,1)= reshape(Cor_VPD_SM_all,[],1);
    elseif p>1 & p<=num_bin+1
        Corr_Bin(:,p)= reshape(Cor_VPD_SM_SMbin(:,:,p-1),[],1);
    else
        Corr_Bin(:,p)= reshape(Cor_VPD_SM_VPDbin(:,:,p-(1+num_bin)),[],1);
    end
end

Corr_Bin=rmmissing(Corr_Bin); %删除有缺失数据的行

figure(1);
set(gcf,'position',[80 100 1120 300]) 
h3=tight_subplot(1,1,[.05 .05],[.25 .15],[.05 .05]);

b1=boxplot(Corr_Bin(:,1),'Colors','k','symbol','','Position',[1],'Widths',0.15);%'PlotStyle','compact'
% b1=violinChart(gca, [1],Corr_Bin(:,1),'Colors',[180 180 180]/255);
hold on
b2=boxplot(Corr_Bin(:,2:11),'Colors','k','symbol','','Position',[1.9:1:10.9],'Widths',0.15); 
hold on
b3=boxplot(Corr_Bin(:,12:21),'Colors','k','symbol','','Position',[2.1:1:11.1],'Widths',0.15); 
hold on

boxobj = findobj(gca,'Tag','Box');
h = findobj(gca,'Tag','Box');
LW = findobj(gca,'Tag','Lower Whisker');
UW = findobj(gca,'Tag','Upper Whisker');
Uav = findobj(gca,'Tag','Upper Adjacent Value');
Lav = findobj(gca,'Tag','Lower Adjacent Value');
M = findobj(gca,'Tag','Median');

color=flipud([[180 180 180]/255;repmat([255 140 0]/255,10,1);repmat([78 98 171]/255,10,1)]);

for j=1:length(h)
    LW(j).LineStyle='-';
    UW(j).LineStyle='-';
    Uav(j).LineStyle='none';
    Lav(j).LineStyle='none';
    M(j).Color='k';
    pc(j)=patch(get(h(j),'XData'),get(h(j),'YData'),color(j,:),'FaceAlpha',0.6);
end

set(b1,'LineWidth',1.2);
set(b2,'LineWidth',1.2);
set(b3,'LineWidth',1.2);

p1=plot([0 12],[0 0],'k--');
% 横坐标 斜体
ax=gca;
set(gca,'xlim',[0.5 11.5],'xtick',[1:1:11]);
ax.XTickLabel={'ALL','0-10 th','10-20 th','20-30 th','30-40 th',...
    '40-50 th','50-60 th','60-70 th','70-80 th','80-90 th','90-100 th'};
xtb = get(gca,'XTickLabel');
xt = get(gca,'XTick');   % 获取横坐标轴刻度句柄
yt = get(gca,'YTick');    % 获取纵坐标轴刻度句柄
ytextp=yt(1)*ones(1,length(xt)); 
xtextp=xt; 
text(xtextp,ytextp-0.23,xtb,'HorizontalAlignment','right','rotation',30,'fontsize',14); 
set(gca,'xticklabel','');% 将原有的标签隐去
        
% xlabel('Percentile','fontsize',18)

% findall(gca,'Tag','Box');
legend(pc([12,1]),'SM bins','VPD bins','fontsize',16,'Orientation','horizon')
% legend boxoff
% 坐标轴加粗
ax=gca;
ax.LineWidth=1.2;
ax.FontSize=13;

text(0.5,0.78,'(d)','fontsize',22);
set(gca,'ylim',[-0.8 0.6],'ytick',[-0.5:0.5:0.5]);

ylabel('r(SM,VPD)','fontsize',18)

print(gcf, '-dpng', '-r600','.\Figure4\FigN1d.png');

%% 逐月计算Zscore的函数
function y = Zscore_month(x)

y=zeros(length(x),1);

for i=1:12
    x1=x(i:12:end);
    x2=zscore(x1);
    
    y(i:12:end,1)=x2;
end

end