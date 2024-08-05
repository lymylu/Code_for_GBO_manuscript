%% main function of the S1 silicon SPK data
datamat=matfile('workspaceofSPK_for_manuscript.mat');
depth=linspace(100,1600,16);
condition_laser={'sound_S1_3.mat','Left_3_5_S1_3.mat','Right_3_5_S1_3.mat','Left_3_S1_3.mat','Right_3_S1_3.mat'};
condition_elec={'Left_1_S1_3.mat','Right_1_S1_3.mat','Left_75_S1_3.mat','Right_75_S1_3.mat'};
% spkt=linspace(-2,4,6/0.01+1); % for the bin size of 10 ms
%% using the SPKsummary in workspaceofSPK_rawraster.mat to plot the raster show;
%using the SPKsummary in workspaceofSPK.mat to plot the gaussian firing
%rate show 
% Laser using 4,4,16, electric using  22, 1 , 16 (Intensity, Side, Both)
index=ismember(SPKsummary.Laserresponse,'Side');
index1=ismember(SPKsummary.putativeCellType,'Interneuron');
index=index&index1;
index=find(index==1);
firingRate=SPKsummary.firingRate(index);
rastermat=matfile('workspaceofSPK_rawraster.mat');
SPKsummary_raster=getfield(rastermat,'SPKsummary');
binmat=matfile('workspaceofSPK_binspike.mat');
SPKsummary_bin=getfield(binmat,'SPKsummary');
%%
figure; 
spikepresent_show(SPKsummary_bin,spkt_bin,condition_laser(2:5),index,[-0.1,0.3]);
%%
figure;
spikepresent_show(SPKsummary_raster,rastermat.spkt',condition_laser(2:5),index,[-0.1,0.3]);
%% %% For Figure S3 C-D Figure S4 C-D Figure 4I Figure 5I
piename={'SPKsummary.Laserresponse','SPKsummary.Electricresponse'};
figure;
for i=1:2
subplot(1,2,i);
x=eval(['tabulate(',piename{i},');']);
Intensity=x{ismember(x(:,1),'Intensity'),2};
Side=x{ismember(x(:,1),'Side'),2};
Both=x{ismember(x(:,1),'Both'),2};
Interaction=x{ismember(x(:,1),'Interaction'),2};
none=x{ismember(x(:,1),'none'),2};
pie([Intensity,Side,Both,Interaction,none],{['Intensity=',num2str(Intensity)],['Side=',num2str(Side)],['Both=',num2str(Both)],['Interaction=',num2str(Interaction)],['none=',num2str(none)]}); 
end
figure;
for i=1:2
subplot(1,2,i);
x=eval(['tabulate(',piename{i},'(ismember(SPKsummary.putativeCellType,''Pyramidal Cell'')));']);
Intensity=x{ismember(x(:,1),'Intensity'),2};
Side=x{ismember(x(:,1),'Side'),2};
Both=x{ismember(x(:,1),'Both'),2};
% Interaction=x{ismember(x(:,1),'Interaction'),2};
% none=x{ismember(x(:,1),'none'),2};
pie([Intensity,Side,Both],{['Intensity=',num2str(Intensity)],['Side=',num2str(Side)],['Both=',num2str(Both)]}); 
end
figure;
for i=1:2
subplot(1,2,i);
x=eval(['tabulate(',piename{i},'(ismember(SPKsummary.putativeCellType,''Interneuron'')));']);
Intensity=x{ismember(x(:,1),'Intensity'),2};
Side=x{ismember(x(:,1),'Side'),2};
Both=x{ismember(x(:,1),'Both'),2};
% Interaction=x{ismember(x(:,1),'Interaction'),2};
% none=x{ismember(x(:,1),'none'),2};
pie([Intensity,Side,Both],{['Intensity=',num2str(Intensity)],['Side=',num2str(Side)],['Both=',num2str(Both)]}); 
end
%%
Celldist=tabulate(SPKsummary.putativeCellType);
m=1;
figure;
for i=1:2
    for j=1:2
    subplot(2,3,m)
    x=eval(piename{i});
    if j==1
        x=strrep(x,'Both','Intensity');
    elseif j==2
        x=strrep(x,'Both','Side');
    end
    y=tabulate(SPKsummary.putativeCellType(ismember(x,type{j})));
    for k=1:2
        typey(k)=y{ismember(y(:,1),Celldist(k,1)),2}/Celldist{k,2};
    end
    bar(typey);
    set(gca,'xticklabel',Celldist(:,1));
    m=m+1;
    title([piename{i},type{j}]);
    end
end
%% sankey plot between function, layer and neuron type For Figure 4H 5H
functiontype={'Intensity','Side','Both'};
layer={'LayerII/III','LayerIV','LayerVa','LayerVb','LayerVI'};
type={'Pyramidal Cell','Interneuron'};
List=[];k=1;
for i=1:length(layer)
    for j=1:length(functiontype)
        a=ismember(SPKsummary.Layer,layer{i});
        b=ismember(SPKsummary.Laserresponse,functiontype{j}); 
        if sum(a&b)>0
        List=cat(1,List,[layer(i),{sum(a&b)},functiontype(j)]);
        k=k+1;
        end
        
    end
end
for i=1:length(functiontype)
    for j=1:length(type)
           a=ismember(SPKsummary.Laserresponse,functiontype{i});
        b=ismember(SPKsummary.putativeCellType,type{j});
        if sum(a&b)>0
        List=cat(1,List,[functiontype(i),{sum(a&b)},type(j)]);
        k=k+1;
        end
    end
end
figure;
sankey2([],'List',List,'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.05);
for i=1:length(layer)
    for j=1:length(functiontype)
        a=ismember(SPKsummary.Layer,layer{i});
        b=ismember(SPKsummary.Electricresponse,functiontype{j}); 
        if sum(a&b)>0
        List=cat(1,List,[layer(i),{sum(a&b)},functiontype(j)]);
        k=k+1;
        end
        
    end
end
for i=1:length(functiontype)
    for j=1:length(type)
           a=ismember(SPKsummary.Electricresponse,functiontype{i});
        b=ismember(SPKsummary.putativeCellType,type{j});
        if sum(a&b)>0
        List=cat(1,List,[functiontype(i),{sum(a&b)},type(j)]);
        k=k+1;
        end
    end
end
figure;
sankey2([],'List',List,'XLim',[0,2],'YLim',[0,1],'PieceWidth',0.05);
%% plot the troughtopeak distribution among spk units for Figure S1 1F
figure; celltype={'Pyramidal Cell','Interneuron'};color={'r','b'};
for i=1:2
    index=ismember(SPKsummary.putativeCellType,celltype{i})&SPKsummary.firingRate>0;
subplot(1,2,1)
histogram(SPKsummary.troughToPeak(index),'BinWidth',0.04); title('VallytoPeak distribution'); hold on;
subplot(1,2,2)
scatter(SPKsummary.troughToPeak(index),SPKsummary.norm_depth(index)); hold on;
adddepthline([0,1]); ylim([0,16]);set(gca,'YDir','reverse');
end


%% main function of the S1 silicon Phase locking data
SpikeTriggeredsummary_elec_ERP=datamat.SpikeTriggeredsummary_elec_Gamma;
SpikeTriggeredsummary_laser_ERP=datamat.SpikeTriggeredsummary_laser_Gamma;
%% for Figure S3 F-H S4 F-H Figure 4 J-L Figure 5 J-L
showt=linspace(-0.1,0.1,201); % showt=linspace(-0.05,0.05,101) for Gamma;
type={'Intensity','Side','Both','none'};
celltype={'Interneuron','Pyramidal Cell'};
for j=1:2
    figure; hold on;
for k=1:3
    tmpdata=eval(['SpikeTriggeredsummary_laser_ERP;']);
    tmpcondition=tmpdata.spikelfp(ismember(SPKsummary.Laserresponse,type{k})&ismember(SPKsummary.putativeCellType,celltype{j}));
    tmpshowall=[];
    tmpshow=tmpcondition;
    for c=1:length(tmpshow)
        tmpshowall=cat(3,tmpshowall,tmpshow{c});
    end
    invalidindex=find(squeeze(nanmean(nanmean(tmpshowall,2),1))==0);
    tmpshowall(:,:,invalidindex)=[];
      amp=squeeze(nanmean(tmpshowall,2));
      amp=squeeze(nanmean(tmpshowall,2));
      amp_all{k}=[];amp_lat{k}=[];
      [amp_all{k},amp_lat{k}]=min(amp(showt>0&showt<0.1,:));
    [~,invalidindex3]=deleteoutliers(amp_all{k},0.05);
    amp_all{k}(invalidindex3)=[];
       showttmp=showt(showt>0&showt<0.025);
    phase=angle(hilbert(squeeze(nanmean(tmpshowall,2))));
    phase_all{k}=phase(showt==0,:);
    subplot(1,4,1); hold on;
    disp(num2str(size(tmpshowall,3)));
    plot(showt,nanmean(nanmean(tmpshowall,3),2));
end
subplot(1,4,2);
circ_plot(phase_all{1}','hist',[],20,true,true,'linewidth',2,'color','r')
subplot(1,4,3);
circ_plot(phase_all{2}','hist',[],20,true,true,'linewidth',2,'color','r')
subplot(1,4,4);
circ_plot(phase_all{3}','hist',[],20,true,true,'linewidth',2,'color','r')
[p_phase1,z_phase1]=circ_rtest(phase_all{1});
[p_phase2,z_phase2]=circ_rtest(phase_all{2});
[p_phase3,z_phase3]=circ_rtest(phase_all{3});
%circ_wwtest(cat(2,phase_all{1},phase_all{2},phase_all{3}),[ones(1,length(phase_all{1})),2*ones(1,length(phase_all{2})),3*ones(1,length(phase_all{3}))]);
if j==1
    phase_all_int=phase_all;
    amp_all_int=amp_all;
    amp_lat_int=amp_lat;
    p_phase_int=[p_phase1,p_phase2,p_phase3];
    z_phase_int=[z_phase1,z_phase2,z_phase3];
else
    phase_all_pyr=phase_all;
    amp_all_pyr=amp_all;
    amp_lat_pyr=amp_lat;
    p_phase_pyr=[p_phase1,p_phase2,p_phase3];
    z_phase_pyr=[z_phase1,z_phase2,z_phase3];
end
end
anovadata=cat(2,amp_all_pyr{1},amp_all_pyr{2},amp_all_pyr{3},amp_all_int{1},amp_all_int{2},amp_all_int{3});
anovadata2=cat(2,amp_lat_pyr{1},amp_lat_pyr{2},amp_lat_pyr{3},amp_lat_int{1},amp_lat_int{2},amp_lat_int{3});
group=cat(1,ones(length(amp_all_pyr{1}),1),2*ones(length(amp_all_pyr{2}),1),3*ones(length(amp_all_pyr{3}),1),4*ones(length(amp_all_int{1}),1),5*ones(length(amp_all_int{2}),1),6*ones(length(amp_all_int{3}),1));
groupneuron=cat(1,ones(length(amp_all_pyr{1}),1),ones(length(amp_all_pyr{2}),1),ones(length(amp_all_pyr{3}),1),2*ones(length(amp_all_int{1}),1),2*ones(length(amp_all_int{2}),1),2*ones(length(amp_all_int{3}),1));
groupfunc=cat(1,ones(length(amp_all_pyr{1}),1),2*ones(length(amp_all_pyr{2}),1),3*ones(length(amp_all_pyr{3}),1),1*ones(length(amp_all_int{1}),1),2*ones(length(amp_all_int{2}),1),3*ones(length(amp_all_int{3}),1));
out=SRH_test([anovadata',groupneuron,groupfunc],'neuron','func'); % SRH test for phase
k=1;amp_stat=cat(2,amp_all_pyr,amp_all_int); statname={'pyr_inten','pyr_side','pyr_both','int_inten','int_side','int_both'};
[mu,sigma]=cellfun(@(x) normfit(x),amp_stat,'UniformOutput',0);
pnorm=cellfun(@(x,y,z) normcdf(z,x,y),mu,sigma,amp_stat,'UniformOutput',0);
gaussianp=cellfun(@(x,y) kstest(x(:),[x(:),y(:)]),amp_stat,pnorm,'UniformOutput',1);
for i=1:5
    for j=i+1:6
        if gaussianp(i)==0&&gaussianp(j)==0
            if vartest2(amp_stat{i},amp_stat{j})>0
            [~,pstat(k)]=ttest2(amp_stat{i},amp_stat{j},'Vartype','unequal');
            else
                [~,pstat(k)]=ttest2(amp_stat{i},amp_stat{j},'Vartype','equal');
            end
        else
        pstat(k)=ranksum(amp_stat{i},amp_stat{j});
        end
        statpair{k,1}=[statname{i},statname{j}];
        k=k+1;
    end
end
pstat=pstat([1,2,6,13,14,15,3,8,12])';
statpair=statpair([1,2,6,13,14,15,3,8,12]);
pstat_SD=[fwer_sidak(pstat(1:3),0.05);fwer_sidak(pstat(4:6),0.05);fwer_sidak(pstat(7:9),0.05)];
pphase=[p_phase_pyr,p_phase_int]';
pphase_SD=[fwer_sidak(pphase,0.05)];
[p_observe]=anovan(anovadata,{groupneuron,groupfunc});
[p_observe2]=anovan(anovadata2,group);
[~,p_observe3]=ttest2(amp_all_pyr{1},amp_all_int{1});
%%
function [phase,mrl,rayleighz,rayleighp,amplitude]=getPhaseAmpfromSTA(triggeredLFP,showt,display)
% triggeredLFP is a matrix, timelfp * channel *neuron 
  invalidindex=squeeze(nanmean(nanmean(triggeredLFP,2),1))==0;
  phase=nan(length(invalidindex),1);
  amplitude=nan(length(invalidindex),1);
  triggeredLFP(:,:,invalidindex)=[];
   amp=squeeze(nanmean(triggeredLFP,2));
    [amplitude(~invalidindex),latency]=min(amp(showt>0,:));
%     invalidindex=amplitude>mean(amplitude)+3*std(amplitude)|amplitude<mean(amplitude)-3*std(amplitude);
    showttmp=showt(showt>0&showt<0.025);
%     amplitude(invalidindex)=[];    
%     latency(invalidindex)=[];  
    tmpphase=angle(hilbert(squeeze(nanmean(triggeredLFP,2))));
    phase(~invalidindex)=tmpphase(showt==0,:);
    [rayleighp,rayleighz]=circ_rtest(phase(~invalidindex)');
    mrl=circ_r(phase(~invalidindex));
    if display==1
        figure; subplot(1,2,1)
        plot(showt,nanmean(nanmean(triggeredLFP,3),2));
        subplot(1,2,2);
        circ_plot(phase(~invalidindex),'hist',[],20,true,true,'linewidth',2,'color','r');
    end
end
function SPKsummary=normalized(SPKsummary,spkt,condition,option)
for i=1:length(condition)
    tmpdata=eval(['SPKsummary.',condition{i}(1:end-4)]);
    tmpdata=cellfun(@(x) basecorrect(mean(x,2),spkt,-2,0,option),tmpdata,'UniformOutput',0);
    eval(['SPKsummary.',condition{i}(1:end-4),'=tmpdata;']);
end
end
function adddepthline(timerange)
  separate=[1.5,5.5,8.5,10.5,13.5];
    hold on;line(gca,[timerange(1),timerange(2)],[separate(1),separate(1)],'LineWidth',1); % separate layer I and II/III
    line(gca,[timerange(1),timerange(2)],[separate(2),separate(2)],'LineWidth',1); % separate layer II/III and layer IV
      line(gca,[timerange(1),timerange(2)],[separate(3),separate(3)],'LineWidth',1); % separate layer IV and layer Va
       line(gca,[timerange(1),timerange(2)],[separate(4),separate(4)],'LineWidth',1,'LineStyle','--'); % separate layer Va and layer Vb
       line(gca,[timerange(1),timerange(2)],[separate(5),separate(5)],'LineWidth',1); % separate layer Vb and layer VI
       set(gca,'YTick',[separate(1)/2,separate(1)+(separate(2)-separate(1))/2, separate(2)+(separate(3)-separate(2))/2,separate(3)+(separate(4)-separate(3))/2,separate(4)+(separate(5)-separate(4))/2,separate(5)+(16-separate(5))/2]);
       set(gca,'YTickLabel',{'Layer I','Layer II/III','Layer IV','Layer Va','Layer Vb','Layer VI'});
       xlabel('troughToPeak');
end
function index=conditionfilter(data,varargin)
p=inputParser;
addParameter(p,'putativeCellType',[]);
addParameter(p,'Layer',[]);
addParameter(p,'firingRate',[]);
addParameter(p,'invalid',[]);
parse(p,varargin{:}{:});
indexType=true(1,length(data.Layer));
indexLayer=true(1,length(data.Layer));
indexfiring=true(1,length(data.Layer));
indexinvalid=true(1,length(data.Layer));
binplotoutput=[];
if ~isempty(p.Results.putativeCellType)
    indexType=ismember(data.putativeCellType,p.Results.putativeCellType);
end
if ~isempty(p.Results.Layer)
    for j=1:length(p.Results.Layer)
        indexLayer(j,:)=ismember(data.Layer,p.Results.Layer{j});
    end
end
if ~isempty(p.Results.firingRate)
    indexfiring=data.firingRate>p.Results.firingRate(1)&data.firingRate<p.Results.firingRate(2);
end
if ~isempty(p.Results.invalid)
    indexinvalid=~p.Results.invalid;
end
if ~isempty(p.Results.Layer)
for j=1:length(p.Results.Layer)
    index{j}=indexType&indexfiring&indexinvalid&indexLayer(j,:);
end
else
    index{1}=indexType&indexfiring&indexinvalid;
end
end
function spikepresent_show(SPKsummary, spkt, condition, index,timelimit)
    for i=1:1:length(condition)
        for c=1:length(index)
        tmpspk=eval(['SPKsummary.',condition{i}(1:end-4),'{index(c)};']);  
%         if i==3; tmpspk(:,6)=[]; end;if i==1; tmpspk(:,3)=[];end
        if ~issparse(tmpspk) 
        tmpspk=basecorrect(mean(tmpspk,2),spkt,-0.1,0,'Subtract');
        tmpspk_all(:,c)=smoothdata(tmpspk,'gaussian',5);
        elseif c==1
            subplot(2,2,i)
            [~,xPoints,yPoints]=plotSpikeRaster(logical(tmpspk)','PlotType','vertline2','TimePerBin',1/40000);
            plot(gca,xPoints/40000-2,yPoints);  axis  tight; xlim(timelimit);
            title(condition{i});
        end
        end
        hold on; plot(spkt,mean(tmpspk_all,2)); xlim(timelimit);
    end
end
