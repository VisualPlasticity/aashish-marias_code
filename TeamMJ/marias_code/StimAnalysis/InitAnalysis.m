function [data] = InitAnalysis(TraceType,bl,block,plane2,sigthresh,plottraces,wantzscore)

if nargin < 7
    wantzscore = 0;
end

if nargin < 6
    plottraces = 0;
end

if nargin < 5
    sigthresh = 1;
end

data.sigthresh = sigthresh;


% Load electrode positions
data.zdiff = 20; % um
data.epos = load('Z:\Teams-MJ-data\epos.mat');


% begin analysis

nblock = length(block);

data.frames = cell(nblock,1);
data.AmplVals = cell(nblock,1);
data.ncell = NaN(nblock,1);
data.cells = cell(nblock,1);
data.Neur = cell(nblock,1);
data.zcells = cell(nblock,1);
data.zNeur = cell(nblock,1);
data.locs = cell(nblock,1);
data.rellocs = cell(nblock,1); % rel locs = location relative to electrode tip
data.meanimage = cell(nblock,1);



for k = 1:nblock
    % get stim type info
    tmpblock = load(block{k});
    frm = tmpblock.info.frame;
    df = diff(frm);
    framestmp = [frm(df>100)+1];
    data.frames{k} = round(framestmp/3);
    
    % get stim Amplitude info
    tmpbl = load(bl{k});
    Ampl = tmpbl.params.Ampl;
    da = diff(Ampl);
    tmpAmplVals = Ampl(find(da>1)+1);
    data.AmplVals{k}=tmpAmplVals(1:end-1);
    
    tmpplane2 = load(plane2{k});
    iscell = logical(cat(1,tmpplane2.stat.iscell));
    tmplocs = cat(1,tmpplane2.stat.med);
    
    tmplocplanes = cat(1,tmpplane2.stat.iplane);

% converted from pixels to microns

    data.locs{k} = (double(tmplocs(iscell,:)) .* repmat(data.epos.micronppixel/data.epos.mag(k),size(tmplocs(iscell,:),1),1));
    data.rellocs{k} = (data.locs{k} - repmat(data.epos.epos(k,1:2),size(data.locs{k},1),1)).* repmat(data.epos.micronppixel/data.epos.mag(k),size(data.locs{k},1),1);
    data.locs{k}(:,3) = (tmplocplanes(iscell))*data.epos.micronpplane;
    data.rellocs{k}(:,3)=(tmplocplanes(iscell) - data.epos.epos(k,3))*data.epos.micronpplane;   

    if strcmp(TraceType,'Neu')
        data.cells{k} = tmpplane2.FcellNeu{1}(iscell,:);
    elseif strcmp(TraceType,'Cell')
        data.cells{k} = tmpplane2.Fcell{1}(iscell,:);
    elseif strcmp(TraceType,'pt7Neu')
        data.cells{k} = tmpplane2.Fcell{1}(iscell,:)-0.7*tmpplane2.FcellNeu{1}(iscell,:);
    elseif strcmp(TraceType,'1Neu')
        data.cells{k} = tmpplane2.Fcell{1}(iscell,:)-tmpplane2.FcellNeu{1}(iscell,:);
    end
    data.Neur{k} = tmpplane2.FcellNeu{1}(iscell,:);
    data.ncell(k) = size(data.cells{k},1);
    
    if wantzscore
        data.zcells{k} = zscore(data.cells{k}')';
        data.zNeur{k} = zscore(data.Neur{k}')';
    else
        data.zcells{k} = data.cells{k};
        data.zNeur{k} = data.Neur{k};
    end
try
    for iplane=1:3
    data.meanimage{k}(:,:,iplane) = tmpplane2.ops{iplane}.mimg1(tmpplane2.ops{iplane}.yrange,tmpplane2.ops{iplane}.xrange);
    end
catch
    for iplane=1:3
    data.meanimage{k}(:,:,iplane) = tmpplane2.ops{iplane}.meanImg(tmpplane2.ops{iplane}.yrange,tmpplane2.ops{iplane}.xrange);
    end
end

data.ncell(k) = size(data.cells{k},1);
end

if plottraces == 1
    for m = 1:nblock
        figure(111+m), clf,hold on
        for k = 1:data.ncell(m)
            plot(data.cells{m}(k,:)+1000*k)
        end
%         for g = 1:length(data.frames{m})
%             plot([data.frames{m}(g) data.frames{m}(g)],[0 1000*data.ncell(m)],'g-')
%             text(data.frames{m}(g),500,num2str(data.AmplVals{m}(g)))
%         end
        title(['Expt. ',m])
        set(gcf,'Position',[100 700 1000 1000])
        saveas(gcf,['Traces_',TraceType,num2str(m),'.pdf'])
    end
end


%% Get cell resps
data.resps = cell(nblock,1);
data.respsorg = cell(nblock,1);
data.maxresp = cell(nblock,1);
data.maxrespidx = cell(nblock,1);
data.maxrespNeur = cell(nblock,1);
data.priorresp = cell(nblock,1);
data.sigact = cell(nblock,1);
data.sigactp = cell(nblock,1);
data.nsamp = 20;
data.amplid = unique(data.AmplVals{1});
data.nampl = length(data.amplid);
%colorlist = colorbrewerRGB(data.nampl, 'qualitative');
c = turbo(data.nampl);
data.AvgResp = cell(nblock,1);
data.VarResp = cell(nblock,1);
data.NeurResp = cell(nblock,1);
data.NeurVar = cell(nblock,1);
data.RespEdges = [-10:1:10];
cc = turbo(length(data.RespEdges));

for m = 1:nblock
    ntrials = length(data.frames{m});
    ncell = size(data.zcells{m},1);
    data.resps{m} = NaN(ntrials,ncell*data.nsamp);
    data.respsorg{m} = NaN(ntrials,ncell,data.nsamp);
    data.maxresp{m} = NaN(ntrials,ncell);
    data.maxrespidx{m} = NaN(ntrials,ncell);
    data.maxrespNeur{m} = NaN(ntrials,ncell);
    data.priorresp{m} = NaN(ntrials,ncell);
    data.sigact{m} = zeros(ntrials,ncell);
    data.sigactp{m} = NaN(data.nampl,ncell);
    
    
    for k = 1:ntrials
        for g = 1:ncell
            data.resps{m}(k,((g-1)*data.nsamp+1):g*data.nsamp) = data.zcells{m}(g,data.frames{m}(k):(data.frames{m}(k)+data.nsamp-1));
            data.respsorg{m}(k,g,:) = data.zcells{m}(g,data.frames{m}(k):(data.frames{m}(k)+data.nsamp-1));
            data.maxresp{m}(k,g) = data.respsorg{m}(k,g,2);
            data.maxrespidx{m}(k,g) = 2;
            data.priorresp{m}(k,g) = mean(data.zcells{m}(g,data.frames{m}(k)-6:data.frames{m}(k)-1));
            tmpNeur = data.zNeur{m}(g,data.frames{m}(k)+1);
            priorNeur = mean(data.zNeur{m}(g,data.frames{m}(k)-6:data.frames{m}(k)-1));
            data.maxrespNeur{m}(k,g) = (tmpNeur-priorNeur)/priorNeur;
            % dF/F is a normalization, tells you how many baselines away
            % the cell changes to
            data.sigact{m}(k,g) = (data.maxresp{m}(k,g) - data.priorresp{m}(k,g))/data.priorresp{m}(k,g);
%             data.sigact{m}(k,g) = (data.maxresp{m}(k,g) - data.priorresp{m}(k,g));
        end
    end
    
    if 0
        % compute p across ALL trials using signrank
        for g = 1:nampl
            for q = 1:ncell
                chtr = find(data.AmplVals{m}==amplid(g));
                priorcell = data.priorresp{m}(chtr,q);
                maxcell = data.maxresp{m}(chtr,q);
                data.sigactp{m}(g,q) = signrank(priorcell,maxcell);
            end % for q = 1:ncell
        end % for g = 1:nampl
    end
    
    
    if 0
        figure(m), imagesc(resps{m})
        xlabel('Time')
        ylabel('Cell #')
        
        %         figure(2000+m)
        %         imagesc(sigactp{m}<0.001)
    end
    
    
    % average over trial types
    data.AvgResp{m} = NaN(data.nampl,ncell);
    data.VarResp{m} = NaN(data.nampl,ncell);
    data.NeurResp{m} = NaN(data.nampl,ncell);
    data.NeurVar{m} = NaN(data.nampl,ncell);
    for k = 1:data.nampl
        data.AvgResp{m}(k,:) = nanmean(data.sigact{m}(data.AmplVals{m}==data.amplid(k),:));
        data.VarResp{m}(k,:) = nanvar(data.sigact{m}(data.AmplVals{m}==data.amplid(k),:));
        data.NeurResp{m}(k,:) = nanmean(data.maxrespNeur{m}(data.AmplVals{m}==data.amplid(k),:));
        data.NeurResp{m}(k,:) = nanvar(data.maxrespNeur{m}(data.AmplVals{m}==data.amplid(k),:));
    end
    %     figure(m)
    %     subplot(1,3,1)
    %     imagesc(AvgResp{m}), colorbar
    %     title(['Average Responses Mouse ',num2str(m)])
    %     xlabel('Cell #')
    %     ylabel('Ampl Val')
    %     subplot(1,3,2)
    %     imagesc(VarResp{m}), colorbar
    %     title(['Variance Responses Mouse ',num2str(m)])
    %     xlabel('Cell #')
    %     ylabel('Ampl Val')
    %     subplot(1,3,3), hold on
    %     for k = 1:data.nampl
    %         plot(AvgResp{m}(k,:),VarResp{m}(k,:),'.','color',colorlist(k,:))
    %     end
    %     legend(num2str(data.amplid));
    %     xlabel('Avg Resp')
    %     ylabel('Var Resp')
    %     title(['Single cells in M',num2str(m)])
    %     set(gcf,'Position',[10 10 1000 300])
    
    %     figure(1000)
    %     subplot(2,3,m), hold on
    %     for k = 1:data.nampl
    %         plot(AvgResp{m}(k,:),VarResp{m}(k,:),'.','color',colorlist(k,:))
    %     end
    %     legend(num2str(data.amplid),'location','west');
    %     xlabel('Avg Resp')
    %     ylabel('Var Resp')
    %     title(['Single cells in M',num2str(m)])
    %     set(gcf,'Position',[10 100 1000 500])
    
    %     figure(1001)
    %     subplot(2,3,m), hold on
    %     for k = 1:data.nampl
    %         plot(AvgResp{m}(k,:),sqrt(VarResp{m}(k,:))./AvgResp{m}(k,:),'.','color',colorlist(k,:))
    %     end
    %     legend(num2str(data.amplid),'location','west');
    %     xlabel('Avg Resp')
    %     ylabel('CV Resp')
    %     title(['Single cells in M',num2str(m)])
    %     set(gcf,'Position',[10 100 1000 500])
    
    
    
    if 0
        figure(10+m)
        for g = 1:data.nampl
            tmpavgresp = data.AvgResp{m}(g,:);
            tmpavgresp(tmpavgresp>10) = 10;
            tmpavgresp(tmpavgresp<-10) = -10;
            BINS = discretize(tmpavgresp,data.RespEdges);
            subplot(2,5,g), hold on
            imagesc(data.meanimage{m})
            for k = 1:ncell
                if data.AvgResp{m}(g,k)>= data.sigthresh
                    plot(data.locs{m}(k,2),data.locs{m}(k,1),'o',...
                        'markersize',tmpavgresp(k)+1,'color',cc(BINS(k),:),...
                        'markerfacecolor',cc(BINS(k),:))
                end
            end
            title(['M',num2str(m),', ',num2str(data.amplid(g)),' uA'])
            axis([61 724 0 506])
            plot(data.epos.epos(m,1),data.epos.epos(m,2),'rd','linewidth',2)
        end
        set(gcf,'Position',[100 500 1200 400])
        saveas(gcf,['Distr_of_Resps_',TraceType,'_Mouse_',num2str(m),'_sigthresh_',num2str(data.sigthresh),'.pdf'])
    end
    
end

% function end


