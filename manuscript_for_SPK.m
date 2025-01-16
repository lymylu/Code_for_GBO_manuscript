%% main function of the S1 silicon SPK data
datamat=matfile('workspaceofSPK_for_manuscript.mat');
depth=linspace(100,1600,16);
condition_laser={'Left_3_5_S1_3.mat','Right_3_5_S1_3.mat','Left_3_S1_3.mat','Right_3_S1_3.mat'};
condition_elec={'Left_1_S1_3.mat','Right_1_S1_3.mat','Left_75_S1_3.mat','Right_75_S1_3.mat'};
% spkt=linspace(-2,4,6/0.01+1); % for the bin size of 10 ms
%% using the SPKsummary in workspaceofSPK_rawraster.mat to plot the raster show;
%using the SPKsummary in workspaceofSPK.mat to plot the gaussian firing
%rate show 
SPKsummary=datamat.SPKsummary;
index=ismember(SPKsummary.Laserresponse,'Side');
%index=ismember(SPKsummary.Electricresponse,'Both');
index=find(index==1);
firingRate=SPKsummary.firingRate(index);
SPKsummary_raster=datamat.SPKsummary_raster;
SPKsummary_bin=datamat.SPKsummary_bin;
%% Figure S3A and Figure S5 A
SPKsummary_raster.spkt=linspace(-2,4,240001);
spikepresent_show(SPKsummary_raster,SPKsummary_raster.spkt',condition_laser,index(4),[-0.1,0.3]);% index(n) for n, Laser using 26,4,16, electric using  22,1,15 (for Intensity, Side, Both)
% SPKsummary_raster.spkt=linspace(-2,4,240001);
% spikepresent_show(SPKsummary_raster,SPKsummary_raster.spkt',condition_elec,index(15),[-0.02,0.05]);% index(n) for n, Laser using 26,4,16, electric using  23,1,15 (for Intensity, Side, Both)
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
%% for Figure S3 F-H S4 F-H Figure 4 J-L Figure 5 J-L
SpikeTriggeredsummary_elec=datamat.SpikeTriggeredsummary_elec_Gamma;
SpikeTriggeredsummary_laser=datamat.SpikeTriggeredsummary_laser_Gamma;
% SpikeTriggeredsummary_elec=datamat.SpikeTriggeredsummary_elec_ERP;
% SpikeTriggeredsummary_laser=datamat.SpikeTriggeredsummary_laser_ERP;
showt=linspace(-0.05,0.05,101); % for ERP showt=linspace(-0.1,0.1,201);
type={'Intensity','Side','Both','none'};
celltype={'Interneuron','Pyramidal Cell'};
for j=1:2
    figure; hold on;
for k=1:3
      tmpdata=eval(['SpikeTriggeredsummary_laser;']);
    tmpcondition=tmpdata.spikelfp(ismember(SPKsummary.Laserresponse,type{k})&ismember(SPKsummary.putativeCellType,celltype{j}));

tmpshow=tmpcondition;
 tmpshowall=[];
for c=1:length(tmpshow)
    tmpshowall=cat(3,tmpshowall,tmpshow{c});
end
invalidindex=find(squeeze(nanmean(nanmean(tmpshowall,2),1))==0); % remove the neurons with no firing in the intervals (LFP all zeros) 
tmpshowall(:,:,invalidindex)=[];
amp=squeeze(nanmean(tmpshowall,2));
  [~,invalidindex2]=deleteoutliers(max(amp)-min(amp),0.01); % for several neurons, due to the low spikes in the interval, the amplitude of STA will be extreme large and should be exclude.
tmpshowall(:,:,invalidindex2)=[];
 amp=squeeze(nanmean(tmpshowall,2));
 amp_all{k}=[];amp_lat{k}=[];
[amp_all{k},amp_lat{k}]=min(amp(showt>0&showt<0.1,:));
    [~,invalidindex3]=deleteoutliers(amp_all{k},0.01);
    amp_all{k}(invalidindex3)=[];
        showttmp=showt(showt>0&showt<0.025);
     phase=angle(hilbert(squeeze(nanmean(tmpshowall,2))));
    phase_all{k}=phase(showt==0,:);
subplot(1,4,1); hold on;
    disp(num2str(size(tmpshowall,3)));
    plot(showt,nanmean(nanmean(tmpshowall,3),2)); 
    if length(showt)==101
        xlim([-0.025,0.025]);
    else
        xlim([-0.1,0.1]);
    end
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
%%
function spikepresent_show(SPKsummary, spkt, condition, index,timelimit)
    for c=1:length(index)
        figure;
        for i=1:1:length(condition)
        
        tmpspk=eval(['SPKsummary.',condition{i}(1:end-4),'{index(c)};']);  
%         if i==3; tmpspk(:,6)=[]; end;if i==1; tmpspk(:,3)=[];end
        if ~issparse(tmpspk) 
        tmpspk=basecorrect(mean(tmpspk,2),spkt,-0.2,0,'zscore');
        tmpspk_all(:,c)=smoothdata(tmpspk,'gaussian',5);
        hold on; plot(spkt,mean(tmpspk_all,2)); xlim(timelimit);
        else
            subplot(2,2,i)
            [~,xPoints,yPoints]=plotSpikeRaster(logical(tmpspk)','PlotType','vertline2','TimePerBin',1/40000);
            plot(gca,xPoints/40000-2,yPoints);  axis  tight; xlim(timelimit);
            title(condition{i});
        end
        end
         
       
    end
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
       xlabel('troughToPeak');
end