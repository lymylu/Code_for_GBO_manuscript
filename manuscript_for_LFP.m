%% main function of the S1 silicon LFP data
% the silicon electrode is ranging from 100 to 1600 under the dura.
global datamat % the data was stored in the datamat
basepath='/media/huli/My Book/YLP/siliconS1/';
subjectnum={'rat1silicon','rat2silicon','rat3silicon','rat4silicon','rat5silicon2','rat6silicon','rat7silicon','rat8silicon','rat9silicon','rat10silicon','rat11silicon','rat12silicon2','rat13silicon','rat14silicon','rat15silicon','rat16silicon'}; % electrode 1-3 in rat9 are broken.  delete 1 and 10 because the electrode modality are bad
lfpt=linspace(-2,4,6001);lfpf=0:100; depth=linspace(100,1600,16);
condition_laser={'Left_3_5_S1_3.mat','Right_3_5_S1_3.mat','Left_3_S1_3.mat','Right_3_S1_3.mat','sound_S1_3.mat'};
condition_elec={'Left_1_S1_3.mat','Right_1_S1_3.mat','Left_75_S1_3.mat','Right_75_S1_3.mat'};
datamat=matfile([basepath,'workspaceofLFP_for_manuscript.mat']);
%% generate plot
%% For Figure S2,S3 and Figure 4B-G Figure 5B-G
varsuffix='_rearrange_correct_blcorrect';
figsavesuffix='rearrange';
tmpvar=eval(['datamat.ERPsummary_laser',varsuffix]);
figure;
plotERP_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,[-0.2,0.5],[-200,50],[1:16],'groupdepthaverage','overlap',[1,10]);
% savefig(gcf,fullfile(basepath,'Figurenew',['ERP_laser_',figsavesuffix,'_averagedepth.fig'])); % Figure 4B
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['ERP_laser_',figsavesuffix,'_averagedepth.eps']),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
% plot the rearranged data averaged from all channels and trials for each condition.
channelindex=[2,5,8,11,14];
figure;
for i=1:5
subplot(1,5,i);plotERP_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,[-0.2,0.5],[-300,100],channelindex(i),'groupdepthaverage','overlap',[1,10]);
legend({'Left High','Righ High','Left Low','Right Low'});
end
% savefig(gcf,fullfile(basepath,'Figurenew',['ERP_laser_',figsavesuffix,'_amongdepth.fig'])); % Figure S2A
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['ERP_laser_',figsavesuffix,'_amongdepth.eps']),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
tmpvar=eval(['datamat.ERPsummary_elec',varsuffix]);
figure;
plotERP_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,[-0.1,0.3],[-300,100],[1:16],'groupdepthaverage','overlap',[1,10]);
% savefig(gcf,fullfile(basepath,'Figurenew',['ERP_elec_',figsavesuffix,'_averagedepth.fig'])); % Figure 5E
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['ERP_elec_',figsavesuffix,'_averagedepth.eps']),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
figure;
for i=1:5
subplot(1,5,i)
plotERP_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,[-0.1,0.3],[-400,100],channelindex(i),'groupdepthaverage','overlap',[1,10]);
legend({'Left High','Righ High','Left Low','Right Low'});
end
% savefig(gcf,fullfile(basepath,'Figurenew',['ERP_elec_',figsavesuffix,'_amongdepth.fig'])); % Figure S3A
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['ERP_elec_',figsavesuffix,'_amongdepth.eps']),'ContentType','Vector');

tmpvar=eval(['datamat.Specsummary_laser',varsuffix]);
plotSpec_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,0:100,depth,[-0.2,0.5],[0,50],[-40,120],[1:16],'groupdepthaverage',[1,10]);
colormap jet;
% savefig(gcf,fullfile(basepath,'Figurenew',['Spec_laser_',figsavesuffix,'_averagedepth.fig'])); %Figure 4E
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['Spec_laser_',figsavesuffix,'_averagedepth.eps']),'ContentType','Vector');
figure;
for i=1:5
plotSpec_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,0:100,depth,[-0.2,0.5],[0,50],[-40,120],channelindex(i),'groupdepthaverage',[1,10]);
colormap jet;
% savefig(gcf,fullfile(basepath,'Figurenew',['Spec_laser_',figsavesuffix,'_depth_',num2str(channelindex(i)*100),'.fig']));%Figure S3B
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['Spec_laser_',figsavesuffix,'_depth_',num2str(channelindex(i)*100),'.eps']),'ContentType','Vector');
end
plotSpec_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,0:100,depth,[-0.2,0.5],[50,100],[-1,3],[1:16],'groupdepthaverage',[1,10]);
colormap jet;
% savefig(gcf,fullfile(basepath,'Figurenew',['GBO_laser_',figsavesuffix,'_averagedepth.fig'])); %Figure 4E
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['GBO_laser_',figsavesuffix,'_averagedepth.eps']),'ContentType','Vector');

figure;
for i=1:5
plotSpec_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,0:100,depth,[-0.2,0.5],[50,100],[-1,3],channelindex(i),'groupdepthaverage',[1,10]);
colormap jet;
% savefig(gcf,fullfile(basepath,'Figurenew',['GBO_laser_rearrange_depth_',num2str(channelindex(i)*100),'.fig']));%Figure S3B
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['GBO_laser_rearrange_depth_',num2str(channelindex(i)*100),'.eps']),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
end
% generate plot Specsummary_elec 
tmpvar=eval(['datamat.Specsummary_elec',varsuffix]);
plotSpec_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,0:100,depth,[-0.1,0.3],[0,100],[-200,600],[1:16],'groupdepthaverage',[1,10]);%%
colormap jet;
% savefig(gcf,fullfile(basepath,'Figurenew',['Spec_elec_',figsavesuffix,'_averagedepth.fig'])); %Figure 5E
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['Spec_elec_',figsavesuffix,'_averagedepth.eps']),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
channelindex=[2,5,8,11,14];
%
figure;
for i=1:5
plotSpec_amongdepth(tmpvar,{'Left High','Right High','Left Low','Right Low'},subjectnum,lfpt,0:100,depth,[-0.1,0.2],[0,100],[-300,900],channelindex(i),'groupdepthaverage',[1,10]);%%
colormap jet;
% savefig(gcf,fullfile(basepath,'Figurenew',['Spec_elec_',figsavesuffix,'_depth_',num2str(channelindex(i)*100),'.fig']));%Figure S4B
% exportgraphics(gcf,fullfile(basepath,'Figurenew',['Spec_elec_',figsavesuffix,'_depth_',num2str(channelindex(i)*100),'.eps']),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
end
%% timepoint anova to find the location and intensity effects on the ERP and SPEC data
   mypool = parpool;
   addAttachedFiles(mypool,'media/huli/My Book/YLP/functions/ylp/neuroview/stat/simple_mixed_anova.m');
[pERPsummary_laser_correct,fERPsummary_laser_correct]=timevarying_anova(datamat.ERPsummary_laser_rearrange_correct_blcorrect,[-0.5,1],[]);
[pERPsummary_elec_correct,fERPsummary_elec_correct]=timevarying_anova(datamat.ERPsummary_elec_rearrange_correct_blcorrect,[-0.5,1],[]);
[pSpecsummary_laser_correct,fSpecsummary_laser_correct]=timevarying_anova(cellfun(@(x) squeeze(nanmean(x(:,60:90,:,:),2)),datamat.Specsummary_laser_correct_blcorrect,'UniformOutput',0),[-0.5,1],[1,10]); 
[pSpecsummary_elec_correct,fSpecsummary_elec_correct]=timevarying_anova(cellfun(@(x) squeeze(nanmean(x(:,60:90,:,:),2)),datamat.Specsummary_elec_correct_blcorrect,'UniformOutput',0),[-0.5,1],[1,10]);
savedatamat(datamat,{'pERPsummary_laser_rearrange_correct_all','pERPsummary_elec_rearrange_correct_all','pSpecsummary_laser_rearrange_correct_all','pSpecsummary_elec_rearrange_correct_all'},{pERPsummary_laser_correct,pERPsummary_elec_correct,pSpecsummary_laser_correct,pSpecsummary_elec_correct}); 
savedatamat(datamat,{'fERPsummary_laser_rearrange_correct_all','fERPsummary_elec_rearrange_correct_all','fSpecsummary_laser_rearrange_correct_all','fSpecsummary_elec_rearrange_correct_all'},{fERPsummary_laser_correct,fERPsummary_elec_correct,fSpecsummary_laser_correct,fSpecsummary_elec_correct});
%%
fdr_correct(datamat.pERPsummary_laser_rearrange_correct_all,[],lfpt,[-0.2,0.15;0.15,0.5],1); % Figure 4D
colormap gray; 
%savefig(gcf,fullfile(basepath,'Figurenew','pERPsummary_laser_rearrange_correct_all.fig'));
%exportgraphics(gcf,fullfile(basepath,'Figurenew','pERPsummary_laser_rearrange_correct_all.eps'),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
fdr_correct(datamat.pSpecsummary_laser_rearrange_correct,[],lfpt,[-0.2,0.15;0.15,0.5],1);% Figure 4G
colormap gray;
%savefig(gcf,fullfile(basepath,'Figurenew','pGBOsummary_laser_rearrange_correct_all.fig'));
%exportgraphics(gcf,fullfile(basepath,'Figurenew','pGBOsummary_laser_rearrange_correct_all.eps'),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
%%
fdr_correct(datamat.pERPsummary_elec_rearrange_correct_all,[],lfpt,[-0.1,0.01;0.01,0.03;0.03,0.3],1);% Figure 5D
colormap gray;
%savefig(gcf,fullfile(basepath,'Figurenew','pERPsummary_elec_rearrange_correct_all.fig'));
%exportgraphics(gcf,fullfile(basepath,'Figurenew','pERPsummary_elec_rearrange_correct_all.eps'),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
fdr_correct(datamat.pSpecsummary_elec_rearrange_correct_all,[],lfpt,[-0.1,0.01;0.01,0.03;0.03,0.3],1);% Figure 5G
colormap gray;
%savefig(gcf,fullfile(basepath,'Figurenew','pGBOsummary_elec_rearrange_correct_all.fig'));
%exportgraphics(gcf,fullfile(basepath,'Figurenew','pGBOsummary_elec_rearrange_correct_all.eps'),'ContentType','Vector');%% plot the rearranged data averaged from all channels and trials for each condition.
%%
 function [datahigh,datalow]=Rearrange(filename,dataup,datadown,eventup,eventdown,subject)
 try   
 [a,b]=xlsread(fullfile(subject,filename));
    indexhigh=a(eventup,3);
    indexlow=a(eventdown,3);
    nanhigh=find(isnan(indexhigh)==1);
    nanlow=find(isnan(indelow)==1);
    [~,index]=sort(cat(1,indexhigh,indexlow));
    orderhigh=zeros(length(indexhigh));
    orderlow=zeros(length(indexlow));
    orderlow=index([1:(length(indexlow)-length(nanlow)),nanlow]);
    orderhigh=index([end-(length(indexhigh)-length(nanhigh))+1:end,nanhigh]);
    if size(dataup,4)==1
     dataall=cat(3,dataup,datadown);
    datalow=dataall(:,:,orderlow);
    datahigh=dataall(:,:,orderhigh);
    else
    dataall=cat(4,dataup,datadown);
    datalow=dataall(:,:,:,orderlow);
    datahigh=dataall(:,:,:,orderhigh);
    end
 catch
     datahigh=dataup;
     datalow=datadown;
 end
end
 function [data,eventindex]=Loaddata(path,file,name,prop)
 % to reduce
switch prop
    case 'origin'
        try
        tmp=matfile(fullfile(path,file));
        tmp=eval(['tmp.',name]);
        data=getfield(tmp,prop);
        catch
            data=nan([6001,16]);
        end
    case 'Spec'
          try
            tmp=matfile(fullfile(path,file));
            tmp=eval(['tmp.',name]);
            data=getfield(tmp,prop);
          catch
                data=nan([6001,101,16]);
          end
end
try 
    index=getfield(tmp.Chooseinfo,'Eventindex');
    eventindex=cellfun(@(x) str2num(x),index,'UniformOutput',1);
catch
    eventindex=[];
end
data=single(data);
 end
function savedatamat(datamat,variablenames,variable)
for i=1:length(variablenames)
    eval(['datamat.',variablenames{i},'=variable{i}']);
end
end
function CSDoutput=csd_cal(varargin)
p=inputParser;
addRequired(p,'data');
addRequired(p,'level');
addParameter(p,'timerange','',@(x) isnumeric(x));
addParameter(p,'type',[],@(x) ischar(x));
addParameter(p,'filter',[],@(x) isnumeric(x));
addParameter(p,'caxis',[],@(x) isnumeric(x));
parse(p,varargin{:});
t=linspace(-2,4,6001);
t2=t(find(t>=p.Results.timerange(1)&t<=p.Results.timerange(2)));
data=p.Results.data;
if ~isempty(p.Results.filter)&&~strcmp(p.Results.type,'spectral')
    for i=1:size(data,3)
        for j=1:size(data,2)
        try
           datafilt(j,:,i)=eegfilt(data(:,j,i)',1000,p.Results.filter(1),p.Results.filter(2));
        catch
           datafilt(j,:,i)=nan(size(data,1),1)';
        end
    end
    end
    data=permute(datafilt,[2,1,3]);
end
if ~isempty(p.Results.timerange)
    data=data(find(t>=p.Results.timerange(1)&t<=p.Results.timerange(2)),:,:,:);
end
% data=basecorrect(data,t,p.Results.timerange(1),0,'subtract');
switch p.Results.type
    case 'origin'
    case 'amplitude'
        for i=1:size(data,3)
            data(:,:,i)=abs(hilbert(data(:,:,i)));
        end
    case 'spectral'
         data=squeeze(nanmean(data(:,p.Results.filter(1):p.Results.filter(2),:,:),2));
end
%      data(isnan(data))=[];
% using iCSD method in CSDplotter
el_pos=(0.1:0.1:1.6)*1e-3;
Fc = F_const(el_pos,2.0e-3,0.3,0.3);
switch p.Results.level
        case 'subjectlevel'
            CSDoutput=nan(6001,16,size(data,3));
        for i=1:size(data,3)
 %          figure;
%             f = fspecial('gaussian',[3 3],0.2);
%             data(:,:,i)=imfilter(data(:,:,i),f,'corr','same');
%          CSDoutput(:,:,i)=CSD(data(:,:,i)./1E6,1000,1E-4,'unitsLength','mm','unitsCurrent','uA','timeaxis',p.Results.timerange,'inverse',1E-3);
 %           close(gcf);
            
            try nandepth=logical(mean(isnan(data(:,:,i))));
             el_pos=(0.1:0.1:1.6)*1e-3;
             el_pos(nandepth)=[];
             Fc = F_const(el_pos,2.0e-3,0.3,0.3);      
             tmpCSD= (Fc^-1*(data(:,~nandepth,i)'/1E6))';
             CSDoutput(:,~nandepth,i)=tmpCSD; 
            end
        end
        case 'grouplevel'
            data=nanmean(data,3);
%             f = fspecial('gaussian',[3 3],0.2);
%             data=imfilter(data,f,'corr','same');
         h=figure;
         CSDoutput=CSD(data./1E6,1000,1E-4,'unitsLength','mm','unitsCurrent','uA','timeaxis',p.Results.timerange,'inverse',1E-3);
         close(h);  
         % CSDoutput=(Fc^-1*data'/1E6)';
end
CSDoutput=basecorrect(CSDoutput,t,-0.1,0,'Subtract');
end
function [pCSDoutput]=pCSD_cal(CSDoutput,timerange,lfpt)
            tmpCSDoutput=permute(CSDoutput,[2,1,3]);
            [pCSDoutput]=sub_tfd_bootstrp(tmpCSDoutput(:,find(lfpt>=timerange(1)&lfpt<=0),:),tmpCSDoutput,200,10);
            pCSDoutput=reshape(mafdr(pCSDoutput(:)),size(pCSDoutput,1),size(pCSDoutput,2));
end
function [CSDsmooth,depth]=CSDplot(varargin) %CSDplot(CSDoutput_laser,condition,subjectnum,lfpt,[-0.1,0.2],'grouplevel','caxis',[-2,2]);
    %% #1 CSDoutput, #2 p, #3,t #4 pCSDoutput
    p=inputParser;
    addRequired(p,'CSDoutput');
    addRequired(p,'condition');
    addRequired(p,'subjectnum');
    addRequired(p,'t');
    addRequired(p,'timerange');
    addRequired(p,'level');
    addParameter(p,'cmap',[],@(x) isnumeric(x));
    addParameter(p,'pCSDoutput',[],@(x) isnumeric(x));
    addParameter(p,'invalid',[],@(x) isnumeric(x));
    parse(p,varargin{:});
    t=p.Results.t;
    cmap=p.Results.cmap;
    timerange=p.Results.timerange;
    CSDoutputall=p.Results.CSDoutput;   
    subjectnum=p.Results.subjectnum;
    subjectnum(p.Results.invalid)=[];
    figure;
    for i=1:length(CSDoutputall)
%         
    switch p.Results.level
        case 'grouplevel'
            CSDoutputall{i}(:,:,p.Results.invalid)=[];
            CSDoutputall{i}=basecorrect(CSDoutputall{i},t,-0.1,0,'subtract');
            CSDoutput=nanmean(CSDoutputall{i},3);
%             subplot(1,2,1);
%             lagging=max(abs(CSDoutput));
%             lagging=cumsum(repmat(max(lagging),[1,size(CSDoutput,2)]));
%             plot(t,bsxfun(@minus,CSDoutput,lagging));
%             axis tight;
%             xlim([timerange(1),timerange(2)]); 
            subplot(1,length(p.Results.condition),i);
            [CSDsmooth{i},depth]=csd_plot(CSDoutput,t,cmap,timerange);
            xlim([timerange(1),timerange(2)]); title(p.Results.condition{i});
        case 'subjectlevel'
            figure;suptitle(p.Results.condition{i});
            CSDoutputall{i}(:,:,p.Results.invalid)=[];
            CSDoutput=CSDoutputall{i};
            for j=1:size(CSDoutputall{i},3)
            subplot(round(sqrt(size(CSDoutput,3))),ceil(size(CSDoutput,3)/round(sqrt(size(CSDoutput,3)))),j);
            [CSDsmooth{i}(:,:,j),depth]=csd_plot(CSDoutput(:,:,j),t,cmap,timerange);
            xlim([timerange(1),timerange(2)]);
            try
                title(subjectnum{j});
            end
            end
    end
    end
end
function [CSDsmooth,zs]=csd_plot(CSDoutput,t,cmap,timerange)
% plot CSD by CSDplotter
gauss_sigma=1e-4;
filter_range=5*gauss_sigma;
CSDoutput=CSDoutput';
el_pos=(0.1:0.1:1.6)*1e-3;
le = length(el_pos);
h = el_pos(2)-el_pos(1);
first_z = el_pos(1)-h/2; %plot starts at z1-h/2;
mfactor = ceil(200/le);
minizs = 0:h/mfactor:(mfactor-1 )*h/mfactor;
for i=1:size(CSDoutput,1) % all rows
zs((1:mfactor)+mfactor*(i-1)) = first_z+(i-1)*h+minizs;
new_CSD_matrix((1:mfactor)+mfactor*(i-1),:) = repmat(CSDoutput(i,:),mfactor,1);
end
[zs,new_CSD_matrix]=gaussian_filtering(zs,new_CSD_matrix,gauss_sigma,filter_range);
[x2,y2]=meshgrid(t,zs.*1e3);
imagesc(gca,t,zs*1e3,new_CSD_matrix); colormap jet;
CSDsmooth=new_CSD_matrix';
% try
% pCSDoutputsmooth=pCSD_cal(new_CSD_matrix,timerange,t);
% % %plot CSD by neuroview
% %     [x,y]=meshgrid(1:size(CSDoutput,1),1:size(CSDoutput,2));
% %     [x2,y2]=meshgrid(1:size(CSDoutput,1),1:0.2:size(CSDoutput,2));
% %      CSDoutput=(interp2(x,y,CSDoutput',x2,y2))';
% %      f = fspecial('gaussian',[3 3],0.2);
% %      CSDoutput=imfilter(CSDoutput,f,'corr','full');
% %      if ~isempty(pCSDoutput)
% %          pCSDoutputsmooth=(interp2(x,y,pCSDoutput,x2,y2))';
% %      end
% %     imagesc(gca,t,y2(:,1),(CSDoutput')); colormap jet;
% %     imagesc(gca,t,y(:,1),CSDoutput');colormap jet;
%         
%         hold on;
%         contour(gca,x2,y2,(pCSDoutputsmooth<0.05),[1,1],'black','LineWidth',1);
%     end
%     set(gca,'ydir','reverse');   
    caxis(cmap); 
    separate=[2,6,8.5,10.5,13]/10;
    hold on;line(gca,[timerange(1),timerange(2)],[separate(1),separate(1)],'LineWidth',1); % separate layer I and II/III
    line(gca,[timerange(1),timerange(2)],[separate(2),separate(2)],'LineWidth',1); % separate layer II/III and layer IV
      line(gca,[timerange(1),timerange(2)],[separate(3),separate(3)],'LineWidth',1); % separate layer IV and layer Va
       line(gca,[timerange(1),timerange(2)],[separate(4),separate(4)],'LineWidth',1,'LineStyle','--'); % separate layer Va and layer Vb
       line(gca,[timerange(1),timerange(2)],[separate(5),separate(5)],'LineWidth',1); % separate layer Vb and layer VI
       set(gca,'YTick',[separate(1)/2,separate(1)+(separate(2)-separate(1))/2, separate(2)+(separate(3)-separate(2))/2,separate(3)+(separate(4)-separate(3))/2,separate(4)+(separate(5)-separate(4))/2,separate(5)+(16-separate(5))/2]);
       set(gca,'YTickLabel',{'Layer I','Layer II/III','Layer IV','Layer Va','Layer Vb','Layer VI'});
       xlabel('Time(s)');
end
function adddepthline(timerange)
  separate=[1.5,5.5,8.5,11,13.5];
    hold on;line(gca,[timerange(1),timerange(2)],[separate(1),separate(1)],'LineWidth',1); % separate layer I and II/III
    line(gca,[timerange(1),timerange(2)],[separate(2),separate(2)],'LineWidth',1); % separate layer II/III and layer IV
      line(gca,[timerange(1),timerange(2)],[separate(3),separate(3)],'LineWidth',1); % separate layer IV and layer Va
       line(gca,[timerange(1),timerange(2)],[separate(4),separate(4)],'LineWidth',1,'LineStyle','--'); % separate layer Va and layer Vb
       line(gca,[timerange(1),timerange(2)],[separate(5),separate(5)],'LineWidth',1); % separate layer Vb and layer VI
       set(gca,'YTick',[separate(1)/2,separate(1)+(separate(2)-separate(1))/2, separate(2)+(separate(3)-separate(2))/2,separate(3)+(separate(4)-separate(3))/2,separate(4)+(separate(5)-separate(4))/2,separate(5)+(16-separate(5))/2]);
       set(gca,'YTickLabel',{'Layer I','Layer II/III','Layer IV','Layer Va','Layer Vb','Layer VI'});
       xlabel('Time (s)');
end
function plotERP_amongdepth(data,condition,subjectnum,lfpt,xrange,yrange,depthindex,option1,option2,invalid)
% option1 grouplevel/subjectlevel/conditionlevel option2 overlap show/separate show  
subjectnum(invalid)=[];
for i=1:length(condition) 
    data{i}(:,:,invalid)=[];
    switch option1
        case 'grouplevel'   
        ERPplot=nanmean(data{i},3);
        lagging=max(abs(ERPplot));
        lagging=cumsum(repmat(max(lagging),[1,16]));
        subplot(1,length(condition),i);
        switch option2
            case 'overlap'
            plot(lfpt,ERPplot);
            case 'separate'
            plot(lfpt,bsxfun(@minus,ERPplot,lagging));
        end
        xlim(xrange);ylim(yrange);title(condition{i})
        case 'subjectlevel'
        ERPplot=data{i};
        figure; suptitle(condition{i})
        for j=1:size(ERPplot,3)
            try
            lagging=max(abs(nanmean(ERPplot(:,:,j),2)));
            lagging=cumsum(repmat(max(lagging),[1,size(ERPplot,2)]));
            subplot(round(sqrt(size(ERPplot,3))),ceil(size(ERPplot,3)/round(sqrt(size(ERPplot,3)))),j);
            switch option2
                case 'separate'
                    plot(lfpt,bsxfun(@minus,zeros(size(ERPplot,1),size(ERPplot,2)),lagging),'black','LineWidth',0.2,'LineStyle',':');
                    hold on;
                plot(lfpt,bsxfun(@minus,ERPplot(:,:,j),lagging));
                case 'overlap'
                plot(lfpt,ERPplot(:,:,j));
            end
            xlim(xrange);ylim(yrange); title(subjectnum{j});
            end
        end
        case 'groupdepthaverage'
             ERPplot=nanmean(nanmean(data{i}(:,depthindex,:),2),3);
              plot(lfpt,ERPplot); hold on;
              xlim(xrange);ylim(yrange);
    end
    end

switch option1
     case 'conditionlevel'
       for i=1:length(data)
        data{i}(:,:,invalid)=[];
       end
        for i=1:16
            figure; suptitle(['depth=',num2str(100*i),'um']);
            for j=1:size(data{1},3)
                subplot(round(sqrt(size(data{1},3))),ceil(size(data{1},3)/round(sqrt(size(data{1},3)))),j);
                hold on; cellfun(@(x) plot(lfpt,x(:,i,j)),data,'UniformOutput',0);
                legend(condition); xlim(xrange);
                title(subjectnum{j});
            end
        end        
end
end
function plotSpec_amongdepth(data,condition,subjectnum,lfpt,lfpf,lfpz,xrange,yrange,crange,depthindex,option1,invalid)
subjectnum(invalid)=[];
    switch option1
        case 'grouplevel' 
         for i=1:length(condition) 
            data{i}(:,:,:,invalid)=[];   
            Specplot(:,:,:,i)=nanmean(data{i},4);
            Specplot=Specplot(:,:,depthindex,i);
         end
            obj=NeuroPlot.imagesc3D;
            obj.plot(lfpt,lfpf,lfpz,Specplot,crange,'xlim',xrange,'ylim',yrange,'subplottitle',condition);
       
        case 'subjectlevel'
            for i=1:length(data)
                data{i}(:,:,:,invalid)=[];   
            Specplot=data{i}(:,:,depthindex,:);
            figure; suptitle(condition{i})
            obj=NeuroPlot.imagesc3D;
            obj.plot(lfpt,lfpf,lfpz,Specplot,crange,'xlim',xrange,'ylim',yrange,'subplottitle',subjectnum,'parent',gcf);
            end
     case 'conditionlevel'
         data{1}(:,:,invalid)=[];
        data{2}(:,:,invalid)=[];
        data{3}(:,:,invalid)=[];
        data{4}(:,:,invalid)=[];
        for i=1:length(depthindex)
            figure; suptitle(['depth=',num2str(100*depthindex(i)),'um']);
            for j=1:size(data{1},3)
                subplot(round(sqrt(size(data{1},3))),ceil(size(data{1},3)/round(sqrt(size(data{1},3)))),j);
                hold on; cellfun(@(x) plot(lfpt,x(:,depthindex(i),j)),data,'UniformOutput',0);
                legend(condition); xlim(xrange);
                title(subjectnum{j});
            end
        end
        case 'groupdepthaverage'
         for i=1:length(condition) 
            data{i}(:,:,:,invalid)=[];   
            Specplot=nanmean(nanmean(data{i}(:,:,depthindex,:),4),3);
%          obj=NeuroPlot.imagesc3D;
%          obj.plot(lfpt,lfpf,1:size(Specplot,3),Specplot,crange,'subplottitle',condition,'xlim',xrange,'ylim',yrange,'parent',gcf);
          subplot(round(sqrt(length(condition))),ceil(length(condition))/round(sqrt(length(condition))),i);
          imagesc(lfpt,lfpf,Specplot'); axis xy; xlim(xrange);ylim(yrange); caxis(crange);
            title(condition{i});
         end
end
end
function [pERP,fERP]=timevarying_anova(summarydata,timerange,invalid)
%% calculate the two way anova among depth of ERP/Spec data intensity&side 2*2 in each time point
if size(summarydata{1},4)~=1
    summarydata{1}=permute(summarydata{1},[1,3,4,2]);
    summarydata{2}=permute(summarydata{2},[1,3,4,2]);
    summarydata{3}=permute(summarydata{3},[1,3,4,2]);
    summarydata{4}=permute(summarydata{4},[1,3,4,2]);
end
    summarydata{1}(:,:,invalid,:)=[];
    summarydata{2}(:,:,invalid,:)=[];
    summarydata{3}(:,:,invalid,:)=[];
    summarydata{4}(:,:,invalid,:)=[];
F_intensity=nan(size(summarydata{1},1),size(summarydata{1},2),size(summarydata{1},4));
P_intensity=nan(size(summarydata{1},1),size(summarydata{1},2),size(summarydata{1},4));
F_Side=nan(size(summarydata{1},1),size(summarydata{1},2),size(summarydata{1},4));
P_Side=nan(size(summarydata{1},1),size(summarydata{1},2),size(summarydata{1},4));
F_interaction=nan(size(summarydata{1},1),size(summarydata{1},2),size(summarydata{1},4));
P_interaction=nan(size(summarydata{1},1),size(summarydata{1},2),size(summarydata{1},4));
lfpt=linspace(-2,4,6001);
multiWaitbar('Processing',0);
 for j=1:size(summarydata{1},2) % for the channel
     for n=1:size(summarydata{1},4) % for the frequency domain if exist
        datamat(:,:,1,1)=summarydata{1}(:,j,:,n);
        datamat(:,:,1,2)=summarydata{2}(:,j,:,n);
        datamat(:,:,2,1)=summarydata{3}(:,j,:,n);
        datamat(:,:,2,2)=summarydata{4}(:,j,:,n);
     parfor i=1:size(summarydata{1},1) % for the time varying.
        if lfpt(i)>=timerange(1)&&lfpt(i)<=timerange(2)
        try
            [tbl,~] = simple_mixed_anova(squeeze(datamat(i,:,:,:)), [],{'Intensity','Side'});
            F_intensity(i,j,n)=tbl{3,4};
            P_intensity(i,j,n)=tbl{3,5};
            F_Side(i,j,n)=tbl{5,4};
            P_Side(i,j,n)=tbl{5,5};
            F_interaction(i,j,n)=tbl{7,4};
            P_interaction(i,j,n)=tbl{7,5}; 
        catch
            F_intensity(i,j,n)=NaN;
            P_intensity(i,j,n)=NaN;
            F_Side(i,j,n)=NaN;
            P_Side(i,j,n)=NaN;
            F_interaction(i,j,n)=NaN;
            P_interaction(i,j,n)=NaN; 
        end
        end
     end
     end
      multiWaitbar('Processing',j/size(summarydata{1},2));
 end
fERP={F_intensity,F_Side,F_interaction};
pERP={P_intensity,P_Side,P_interaction};
end
function pcorrect=sub_fdr_correct(p)
    [m,n]=size(p);
    for i=1:n
    pcorrect(:,i)=fdr_BH(p(:,i),0.05);
    end
end
function fdr_correct(pERPsummary,fERPsummary,lfpt,timeband,option)
if option
 P_intensity=[];P_side=[];P_interaction=[];
    for k=1:size(timeband,1)
        P_intensity=cat(1,P_intensity,sub_fdr_correct(pERPsummary{1}(lfpt>timeband(k,1)&lfpt<timeband(k,2),:)));
        P_side=cat(1,P_side,sub_fdr_correct(pERPsummary{2}(lfpt>timeband(k,1)&lfpt<timeband(k,2),:)));
        P_interaction=cat(1,P_interaction,sub_fdr_correct(pERPsummary{3}(lfpt>timeband(k,1)&lfpt<timeband(k,2),:)));
    end
    index=zeros(1,length(lfpt));
    for k=1:size(timeband,1)
        index=index|lfpt>timeband(k,1)&lfpt<timeband(k,2);
    end
    lfpt=lfpt(index);
    timeband=reshape(timeband,[],1);
    timeband=[min(timeband),max(timeband)];
    figure; subplot(2,3,1);imagesc(lfpt,1:16,P_intensity');xlim(timeband);caxis([0,0.05]);title('Intensity P');set(gca,'FontSize',14);
    adddepthline(timeband);
    subplot(2,3,2);imagesc(lfpt,1:16,P_side');xlim(timeband);caxis([0,0.05]);title('Side P'); adddepthline(timeband);set(gca,'FontSize',14);
    subplot(2,3,3);imagesc(lfpt,1:16,P_interaction');xlim(timeband);caxis([0,0.05]);title('Side*Intensity P'); adddepthline(timeband);set(gca,'FontSize',14);
else
    timeband=reshape(timeband,[],1);
    timeband=[min(timeband),max(timeband)];
    figure; subplot(1,3,1);imagesc(lfpt,1:16,pERPsummary{1}');xlim(timeband);caxis([0,0.05]);title('Intensity P'); adddepthline(timeband);set(gca,'FontSize',14);
    subplot(1,3,2);imagesc(lfpt,1:16,pERPsummary{2}');xlim(timeband);caxis([0,0.05]);title('Side P'); adddepthline(timeband);set(gca,'FontSize',14);
    subplot(1,3,3);imagesc(lfpt,1:16,pERPsummary{3}');xlim(timeband);caxis([0,0.05]);title('Side*Intensity P'); adddepthline(timeband);set(gca,'FontSize',14);
end
    timeband=reshape(timeband,[],1);
    timeband=[min(timeband),max(timeband)];
    end
function len=channellength(x)
if ischar(x)
    len=length(str2num(x));
elseif isnumeric(x)
    if isnan(x)
        len=0;
    else
        len=length(x);
    end
end
end