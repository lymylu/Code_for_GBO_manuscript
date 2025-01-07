% close all
% clc
datamataverage=matfile('EEGdata_for_manuscript.mat','Writable',true);
%% pool all subject
% left and right pool all trials into one matrix
 tmpscore=catsubject(datamataverage,'score_all',1);
 tmplevel=catsubject(datamataverage,'level_all',1);
 [tmpGBO_Cz,subject_Cz,trial_Cz]=catsubject(datamataverage,'GBO_Cz',3);
 datamataverage.tmpGBO_Cz=tmpGBO_Cz;
 datamataverage.subject_Cz=subject_Cz;
 datamataverage.trial_Cz=trial_Cz;
%% get normalized GBO and score from level.
[tmpGBO_Cz_partial,subject_Cz,trial_Cz,tmpscore_partial]=catsubject_partial(datamataverage,'GBO_Cz','level_all','score_all');
datamataverage.tmpGBO_Cz_partial=tmpGBO_Cz_partial;
datamataverage.tmpscore_partial=tmpscore_partial;
datamataverage.subject_Cz=subject_Cz;
datamataverage.trial_Cz=trial_Cz;
%% LMM analysis
%% partial LMM 
tmpscore_partial=datamataverage.tmpscore_partial; 
tmpGBO_Cz_partial=datamataverage.tmpGBO_Cz_partial;
scoreall=cat(1,tmpscore{1},tmpscore{2});validindex=true(size(scoreall));
[T,P]=LMM.LMM_analysis_bilateral(subject_Cz,trial_Cz,tmpscore_partial,cellfun(@(x) permute(x(lfpt>-0.2&lfpt<0.5,:,:),[3,1,2]),tmpGBO_Cz_partial,'UniformOutput',0),validindex);
datamataverage.T_partial_LMM_bilateral_y_score=T;
datamataverage.P_partial_LMM_bilateral_y_score=P;
%% score LMM
tmpscore=datamataverage.tmpscore; 
tmpGBO_Cz=datamataverage.tmpGBO_Cz;
scoreall=cat(1,tmpscore{1},tmpscore{2});validindex=true(size(scoreall));
[T,P]=LMM.LMM_analysis_bilateral(subject_Cz,trial_Cz,tmpscore,cellfun(@(x) permute(x(lfpt>-0.2&lfpt<0.5,:,:),[3,1,2]),tmpGBO_Cz,'UniformOutput',0),validindex);
datamataverage.T_LMM_bilateral_y_score=T;
datamataverage.P_LMM_bilateral_y_score=P;
%% level LMM
tmpGBO_Cz=datamataverage.tmpGBO_Cz;
scoreall=cat(1,tmpscore{1},tmpscore{2});validindex=true(size(scoreall));
[T,P]=LMM.LMM_analysis_bilateral(subject_Cz,trial_Cz,tmplevel,cellfun(@(x) permute(x(lfpt>-0.2&lfpt<0.5,:,:),[3,1,2]),tmpGBO_Cz,'UniformOutput',0),validindex);
datamataverage.T_level_LMM_y_level=T;
datamataverage.P_level_LMM_y_level=P;
%% ROI select from T and P value (combined right and left trials)
variablesuffix={'_LMM_bilateral_y_score','_partial_LMM_bilateral_y_score','_level_LMM_y_level'};
roisuffix={'_score','_partial','_level'};
lfpt=linspace(-1,2,3001);
lfptroi=lfpt(lfpt>-0.2&lfpt<0.5);
for m=1:3
    T=eval(['datamataverage.T',variablesuffix{m}]);
    P=eval(['datamataverage.P',variablesuffix{m}]);
    tfroi{m}=roiselect(T,P,lfptroi,[0,0.5],[50,100]);
    if m==2
        tmpGBO=cellfun(@(x) x(lfpt>-0.2&lfpt<0.5,:,:),tmpGBO_Cz_partial,'UniformOutput',0);
    else
        tmpGBO=cellfun(@(x) x(lfpt>-0.2&lfpt<0.5,:,:),tmpGBO_Cz,'UniformOutput',0); 
    end
    tmpGBO_roi=cellfun(@(x) squeeze(nanmean(nanmean(x.*repmat(tfroi{m},[1,1,size(x,3)]),1),2)),tmpGBO,'UniformOutput',0);
    %eval(['datamataverage.GBOroi',roisuffix{m},'=tmpGBO_roi;']);
end
%% 
figure;
for m=1:3
    T=eval(['datamataverage.T',variablesuffix{m}]);
    P=eval(['datamataverage.P',variablesuffix{m}]);
    [a]=fdr_BH(P,0.05);
    Ttmp=T(:); Ttmp=Ttmp(:);
    sigT(m)=max(Ttmp(a>0.05-0.0001&a<0.05+0.0001));
    Ttmp2=squeeze(mean(T(:,50:100),2));
    plot(lfptroi,Ttmp2);
    hold on;
    plot(lfptroi,repmat(sigT(m),[1,length(lfptroi)]));
    sigtime(m)=min(lfptroi(find(Ttmp2>sigT(m))));
end
savefig(gcf,'averageLMMT.fig'); % figure 3B
%% figure 3C, shuffle the data to generate the distribution of sigtime.
% *_shuffle_LMM_1.mat contains 140 shuffles (20 times * 7 [for
% 90%,80%,70%,60%,50%,40% and 30% of the origin data]), in this manuscript, using 80% of the data.
for m=1:3
    shufflemat=matfile(fullfile(basepath,[roisuffix{m}(2:end),'_shuffle_LMM_1.mat']));
    sigtime_shuffle=[];
    for s=1:140
        T=eval(['shufflemat.Shuffle_T_',num2str(s),';']);
        P=eval(['shufflemat.Shuffle_P_',num2str(s),';']);
        [a]=LMM.fdr_correct(P,0.05); a=a(:);
        Ttmp=T(:);
%         a(a>0.05)=1;a(a<0.05)=0;
%         sigtime_shuffle(m,s)=min(lfptroi(find(mean(a(:,50:100),2)<0.8)));
        sigT=min(Ttmp(a>0.05&a<0.05+0.0001));
        Ttmp2=squeeze(mean(T(:,50:100),2));
        try
        sigtime_shuffle(s)=min(lfptroi(find(Ttmp2>sigT&Ttmp2>0)));
        catch
            sigtime_shuffle(s)=nan;  
        end
    end
    sigtime_shuffletmp=reshape(sigtime_shuffle,7,20);
    sigtime_shuffle_all(m,:)=sigtime_shuffletmp(2,:);
end
%% Figure 3D LMM topoplot (reanalysis the LMM from the roi of GBO among channels to get the topographic distribution)
titlelist={'score','score|level','level'};
figure;
    GBOroi=catsubject(datamataverage,['GBOroi',roisuffix{1}],2);
    [T,P]=LMM.LMM_analysis_bilateral(subject_Cz,trial_Cz,tmpscore,cellfun(@(x) x',GBOroi,'UniformOutput',0),validindex);
    subplot(3,2,1); topoplot(T,EEG.chanlocs);caxis([0,12]); title('power~score');
    pcorrect=fdr_correct(P,0.05);
    subplot(3,2,2); topoplot(pcorrect,EEG.chanlocs);caxis([0,0.01]);
    GBOroi=catsubject(datamataverage,['GBOroi',roisuffix{2}],2);
     [GBOroi,subject_Cz,trial_Cz,tmpscore_partial]=catsubject_partial(datamataverage,['GBOroi',roisuffix{3}],'level_all','score_all');
    [T,P]=LMM.LMM_analysis_bilateral(subject_Cz,trial_Cz,tmpscore_partial,cellfun(@(x) x',GBOroi,'UniformOutput',0),validindex);
    subplot(3,2,3); topoplot(T,EEG.chanlocs);caxis([0,10]); title('(power~level)~(score~level)');
     pcorrect=fdr_correct(P,0.05);
    subplot(3,2,4); topoplot(pcorrect,EEG.chanlocs);caxis([0,0.01]);
     GBOroi=catsubject(datamataverage,['GBOroi',roisuffix{3}],2);
    [T,P]=LMM.LMM_analysis_bilateral(subject_Cz,trial_Cz,tmplevel,cellfun(@(x) x',GBOroi,'UniformOutput',0),validindex);
    subplot(3,2,5); topoplot(T,EEG.chanlocs);caxis([0,10]);title('power~level');
    pcorrect=fdr_correct(P,0.05);
    subplot(3,2,6); topoplot(pcorrect,EEG.chanlocs);caxis([0,0.01]);
%% Figure 3E
close all;
for m=1:3
    invalidscore=false(1,length(tmpscore{1})+length(tmpscore{2}));
    tmpGBO_roi=catsubject(datamataverage,['GBOroi',roisuffix{m}],2);
    if m==3
    [scoreindex,scorelegend]=divideindex(cat(1,tmplevel{1},tmplevel{2}),1:4,false);
    elseif m==2
    [tmpGBO_roi,subject_Cz,trial_Cz,tmpscore_partial]=catsubject_partial(datamataverage,['GBOroi',roisuffix{m}],'level_all','score_all'); 
    tmpscore_all=cat(1,tmpscore{1},tmpscore{2});
    tmpscore_partial_all=cat(1,tmpscore_partial{1},tmpscore_partial{2});
    invalidscore=invalidscore(tmpscore_all>=4);
    [scoreindex,scorelegend]=divideindex(tmpscore_all(~invalidscore),[4:9;5:10],false);
    else    
        [scoreindex,scorelegend]=divideindex(cat(1,tmpscore{1},tmpscore{2}),[4:9;5:10],false);
    end
    tmpGBO_roi_all=cat(2,tmpGBO_roi{1}(13,:),tmpGBO_roi{2}(13,:)); % select Cz channel
    tmpGBO_roi_all=tmpGBO_roi_all(~invalidscore)';
    figure;
    boxGBO=nan(length(scoreindex{1}),length(scoreindex));
    amounttrial=[];
    for j=1:length(scoreindex)
        boxGBO(scoreindex{j},j)=tmpGBO_roi_all(scoreindex{j});
        amounttrial(j)=sum(~isnan(boxGBO(:,j)));
    end
    subplot(1,2,1);
    x=1:length(scoreindex);
    errorbar(x,nanmean(boxGBO,1),std(boxGBO,'omitnan')./sqrt(amounttrial)); title(['linearfit',roisuffix{m}]); %ylim([-0.02,0.05]); 
    xlim([0.5,j+0.5]);
    [f1,gof1]=fit(x',nanmean(boxGBO,1)','poly1');
    hold on; plot(f1,x,nanmean(boxGBO,1),'b','predfunc');
    
    if m==2
        text(1,0.15,struct2str(gof1));
    else
        text(1,0.02,struct2str(gof1));
    end
    subplot(1,2,2);
     errorbar(x,nanmean(boxGBO,1),std(boxGBO,'omitnan')./sqrt(amounttrial)); title(['exponentialfit',roisuffix{m}]); %ylim([-0.02,0.05]); 
     xlim([0.5,j+0.5]);
    if m==2
    ft = fittype('a*exp(b*x)+c');
        [f2,gof2]=fit(x',nanmean(boxGBO,1)',ft,'Startpoint',[1,2,-1]);
        hold on; plot(f2,x,nanmean(boxGBO,1),'b','predfunc');
        text(1,0.15,struct2str(gof2));
    else
        [f2,gof2]=fit((1:length(scoreindex))',nanmean(boxGBO,1)','exp1','Startpoint',[1,2]);
        hold on; plot(f2,1:length(scoreindex),nanmean(boxGBO,1),'b','predfunc');
        text(1,0.02,struct2str(gof2));
    end 
end
%% Figure 2
tmptmp_gamma=datamataverage.tmptmp_gamma;
C3spec=datamataverage.C3gamma;
C4spec=datamataverage.C4gamma;
tmpmat=matfile('EEGdata_for_manuscript_Fz.mat'); % for ERP, the data were rereferenced to Fz channel.
tmptmp=tmpmat.tmptmp;
EEGFz=tmpmat.EEG;
EEG=datamataverage.EEG;
location={'Right','Left'};
% plot TFD of C3 and C4 channel
figure;c=1;
for i=1:length(location)
        subplot(4,2,c)
        imagesc(lfpt,1:50,squeeze(mean(C3spec(i,:,:,(1:50)),2))');c=c+1;xlim([-0.2,0.5]); caxis([-1,2]); axis xy; title(['C3 low ',location{i}]);
        subplot(4,2,c)
        imagesc(lfpt,50:100,squeeze(mean(C3spec(i,:,:,(50:100)),2))');c=c+1;xlim([-0.2,0.5]); caxis([-0.01,0.02]);axis xy;title(['C3 high ',location{i}]);
        subplot(4,2,c)
        imagesc(lfpt,1:50,squeeze(mean(C4spec(i,:,:,(1:50)),2))');c=c+1;xlim([-0.2,0.5]); caxis([-1,2]); axis xy;title(['C4 low ',location{i}]);
        subplot(4,2,c)
        imagesc(lfpt,50:100,squeeze(mean(C4spec(i,:,:,(50:100)),2))');c=c+1;xlim([-0.2,0.5]); caxis([-0.01,0.02]);axis xy;title(['C4 high ',location{i}]);
end
savefig(gcf,'AverGBO.fig');close(gcf); % figure 2C
figure; 
c=1;
for i=1:length(location)
            tinterval=lfpt>=0.15&lfpt<0.2;
            tmpplot=squeeze(mean(tmptmp_gamma(i,:,:,tinterval,:),[2,4,5])); 
            subplot(1,2,c); topoplot(tmpplot,EEG.chanlocs);  title([location{i}]); caxis([-0.05,0.05]);
            c=c+1;
end
savefig(gcf,'TopoGBO.fig');close(gcf); % figure 2C
channellabel=struct2table(EEGFz.chanlocs);
channellabel=channellabel.labels;
for i=1:length(location)
        tmpplot=squeeze(mean(tmptmp(i,:,ismember(channellabel,'C3'),:,:),2));
        subplot(1,2,1);hold on;plot(lfpt,mean(tmpplot,2)); title('C3 channel'); xlim([-0.5,1]);ylim([-25,25]); set(gca,'YDir','reverse');
        tmpplot=squeeze(mean(tmptmp(i,:,ismember(channellabel,'C4'),:,:),2));
        subplot(1,2,2);hold on;plot(lfpt,mean(tmpplot,2)); title('C4 channel'); xlim([-0.5,1]);ylim([-25,25]);set(gca,'YDir','reverse');
end
savefig(gcf,'AverERP.fig');close(gcf); % figure 2B
figure;  tindex=lfpt>0.15&lfpt<0.16; c=1;
for i=1:length(location)
        tmpplot=squeeze(mean(tmptmp(i,:,:,tindex,:),2));
        subplot(1,2,c); topoplot(squeeze(mean(mean(tmpplot,2),3)),EEGFz.chanlocs); caxis([-15,15]); title([location{i}]);
        c=c+1;
end
savefig(gcf,'TopoERP.fig');close(gcf); % figure 2B

% % %
function [tmpdata,subject,trial,partialscore]=catsubject_partial(datamat,datavarname,levelvarname,scorename)
    data=evalvar(datamat,datavarname);% spectrogram data 2*95 cell time*frequency*trial
    level=evalvar(datamat,levelvarname); % level data 2*95 cell, trial*1 range 1/2/3/4
    score=evalvar(datamat,scorename); % score data 2*95 cell, trial*1 range 0-10
    tmpdata=cell(1,size(data,1));subject=cell(1,size(data,1));
    trial=cell(1:size(data,1));
    partialscore=cell(1:size(data,1)); % normalized score
    for j=1:size(data,1) % 2 right and left
        for i=1:size(data,2) % 95 subject
            tmplevel=level{j,i};
            for c=1:4
                index=tmplevel==c; % get each rank of level index.
                if size(data{j,i},3)>1
                data{j,i}(:,:,index)=(data{j,i}(:,:,index)-repmat(squeeze(mean(data{j,i}(:,:,index),3)),[1,1,size(data{j,i}(:,:,index),3)]))...
                    ./repmat(squeeze(std(data{j,i}(:,:,index),0,3)),[1,1,size(data{j,i}(:,:,index),3)]);
                else
                       data{j,i}(:,index)=(data{j,i}(:,index)-repmat(squeeze(mean(data{j,i}(:,index),2)),[1,size(data{j,i}(:,index),2)]))...
                    ./repmat(squeeze(std(data{j,i}(:,index),0,2)),[1,size(data{j,i}(:,index),2)]);
                end
                score{j,i}(index)=(score{j,i}(index)-mean(score{j,i}(index)))/std(score{j,i}(index));
            end
            if size(data{j,i},3)>1
            tmpdata{j}=cat(3,tmpdata{j},data{j,i}); 
            subject{j}=cat(1,subject{j},ones(size(data{j,i},3),1)*i); % subject label
            trial{j}=cat(2,trial{j},1:size(data{j,i},3)); % trial label
            partialscore{j}=cat(1,partialscore{j},score{j,i}); % cat the normalized score.
            else
            tmpdata{j}=cat(2,tmpdata{j},data{j,i});
            subject{j}=cat(1,subject{j},ones(size(data{j,i},2),1)*i); % subject label
            trial{j}=cat(2,trial{j},1:size(data{j,i},2)); % trial label
            partialscore{j}=cat(1,partialscore{j},score{j,i}); % cat the normalized score.
            end
        end
    end
end
function [tmpdata,subject,trial]=catsubject(datamat,varname,catdimension)
    data=evalvar(datamat,varname);
    tmpdata=cell(1,size(data,1));subject=cell(1,size(data,1));
    trial=cell(1:size(data,1));
    for j=1:size(data,1)
    for i=1:size(data,2)
        tmpdata{j}=cat(catdimension,tmpdata{j},data{j,i});
        subject{j}=cat(1,subject{j},ones(size(data{j,i},catdimension),1)*i);
        trial{j}=cat(2,trial{j},1:size(data{j,i},catdimension));
    end
    end
end
function [index,legend]=divideindex(data,window,prc)
if size(window,1)==2
    timebegin=window(1,:);
    timeend=window(2,:);
    for i=1:length(window)
        if i==length(window)
            if prc
                index{i}=data>=prctile(data,timebegin(i))&data<=prctile(data,timeend(i));
            else
                index{i}=data>=timebegin(i);
            end
        legend{i}=strcat('>=',num2str(timebegin(i)),'<=',num2str(timeend(i)));
        else
            if prc
                index{i}=data>=prctile(data,timebegin(i))&data<prctile(data,timeend(i));
            else
                index{i}=data>=timebegin(i)&data<timeend(i);
            end
             legend{i}=strcat('>=',num2str(timebegin(i)),'<',num2str(timeend(i)));
        end
    end
else
    for i=1:length(window)
        index{i}=data==window(i);
          legend{i}=strcat('=',num2str(window(i)));
    end
end
end
function data=evalvar(datamat,varname)
    data=eval(['datamat.',varname]);
end
function pcorrect=fdr_correct(p,alpha)
    %[pcorrect]=fdr_TN(p,alpha);
    dim=size(p);
    [p]=fdr_BH(p(:),alpha);
    pcorrect=reshape(p,dim);
end
function tfroi=roiselect(T,P,lfptroi,timerange,freqrange)
    % select the roi according to the overlap of significant P and given
    % time-frequnecy range.
        figure;
        pcorrect=fdr_correct(P,0.05);
        imagesc(lfptroi,1:100,T');axis xy;caxis([-2,8]);
        hold on; contour(lfptroi,1:100,pcorrect',[0.01,0.01]);colormap jet;
        v=nan(size(pcorrect));
        v(pcorrect<0.05)=1;
        roi=nan(size(v));
        roi(lfptroi>timerange(1)&lfptroi<timerange(2),freqrange(1):freqrange(2))=1;
        tfroi=roi.*v;
end
function output=struct2str(gof)
    varname=fieldnames(gof);
    for i=1:length(varname)
        output{i}=[varname{i},num2str(eval(['gof.',varname{i}]))];
    end
end
