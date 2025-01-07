classdef LMM
    % function package for LMM and its stability
    properties
        
    end
    methods(Static)
      function [T,P]=LMM_analysis_bilateral(subject_all,trial_all,score_all,data_all)
    T=nan(size(data_all,2),size(data_all,3));
    P=nan(size(data_all,2),size(data_all,3));
    T=T(:);
    P=P(:);
    data_all_calculate=reshape(data_all,size(data_all,1),[]);
    parfor i=1:size(data_all_calculate,2)
        tbl=table(subject_all,trial_all',score_all,data_all_calculate(:,i),'VariableNames',{'subject','trial','score','power'});
        lme=fitlme(tbl,'score~1+power+(1+power|subject)','FitMethod','REML','CovariancePattern','Diagonal');
        result=dataset2table(lme.Coefficients);
        T(i)=result.tStat(2);
        P(i)=result.pValue(2);
    end
    T=reshape(T,size(data_all,2),size(data_all,3));
    P=reshape(P,size(data_all,2),size(data_all,3));
end
function LMM_analysis_stability(subject,trial,score,data,validindex,savemat)
     subject_all=cat(1,subject{1},subject{2});
     trial_all=cat(2,trial{1},trial{2});
     score_all=cat(1,score{1},score{2});
     data_all=cat(1,data{1},data{2});
     varnames=fieldnames(savemat);
     multiWaitbar('Process the suffle LMM',0);
     for i=1:length(validindex)
         subject=subject_all(validindex{i});
         score=score_all(validindex{i});
         trial=trial_all(validindex{i});
         data=data_all(validindex{i},:,:);
          if ~ismember(varnames,['Shuffle_T_',num2str(i)])
              [T,P]=LMM.LMM_analysis_bilateral(subject,trial,score,data);
              eval(['savemat.Shuffle_T_',num2str(i),'=T;']);
              eval(['savemat.Shuffle_P_',num2str(i),'=P;']);
          end
          multiWaitbar('Process the shuffle LMM',i/length(validindex));
     end
end
function tfroi=roiselect(P,lfptroi,timerange,freqrange,color)
    % select the roi according to the overlap of significant P and given
    % time-frequnecy range.
        %figure;
        pcorrect=LMM.fdr_correct(P,0.05);
        %imagesc(lfptroi,1:100,T');axis xy;caxis([-2,8]);
        if nargin<5
        %hold on; contour(lfptroi,1:100,pcorrect',[0.01,0.01]);
        else
            %hold on; contour(lfptroi,1:100,pcorrect',[0.01,0.01],'linewidth',4,'color',color);
        end
        v=nan(size(pcorrect));
        v(pcorrect<0.05)=1;
        roi=nan(size(v));
        roi(lfptroi>timerange(1)&lfptroi<timerange(2),freqrange(1):freqrange(2))=1;
        tfroi=roi.*v;
end
function pcorrect=fdr_correct(p,alpha)
    %[pcorrect]=fdr_TN(p,alpha);
    dim=size(p);
    [p]=fdr_BH(p(:),alpha);
    pcorrect=reshape(p,dim);
end
function [tmpdata,subject,trial]=catsubject(datamat,varname,catdimension)
    data=LMM.evalvar(datamat,varname);
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
function data=evalvar(datamat,varname)
    data=eval(['datamat.',varname]);
end
function r=cal_ICC(T_origin,T_shuffle)
        r=[];
        shufflenum=size(T_shuffle,1);
        T_shuffle=T_shuffle(:);
        T_shuffle(isnan(T_shuffle))=0;
        T_shuffle=reshape(T_shuffle,shufflenum,[]);
        T_origin=T_origin(:);
        for i=1:shufflenum
            r(i)=ICC([squeeze(T_shuffle(i,:))',T_origin],'1-1');
        end
        r=reshape(r',shufflenum,[]);
end
function catdata(file,fileadded)
    filemat=matfile(file);
    for i=1:50
        T_tmp(i,:,:)=eval(['filemat.Shuffle_T_',num2str(i),';']);
        P_tmp(i,:,:)=eval(['filemat.Shuffle_P_',num2str(i),';']);
    end
    T_reshape=reshape(T_tmp,5,10,700,100);
    P_reshape=reshape(P_tmp,5,10,700,100);
    filemat2=matfile(fileadded);
    for i=1:20
        T_add(i,:,:)=eval(['filemat2.Shuffle_T_',num2str(i),';']);
        P_add(i,:,:)=eval(['filemat2.Shuffle_P_',num2str(i),';']);
    end
    T_reshape_add=reshape(T_add,2,10,700,100);
    P_reshape_add=reshape(P_add,2,10,700,100);
    T_all=cat(1,T_reshape,T_reshape_add);
    P_all=cat(1,P_reshape,P_reshape_add);
    T_all=reshape(T_all,[],700,100);
    P_all=reshape(P_all,[],700,100);
    filematall=matfile([file(1:end-4),'all.mat'],'Writable',true);
    for i=1:size(T_all,1)
        eval(['filematall.Shuffle_T_',num2str(i),'=squeeze(T_all(i,:,:));']);
        eval(['filematall.Shuffle_P_',num2str(i),'=squeeze(P_all(i,:,:));']);
    end
end
    end
end