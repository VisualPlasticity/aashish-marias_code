% ANALYZE STIM-EVOKED DATA
% changed distance to 3D in microns - 4th August 2021 - JA
clear all
close all
clc

% TRACE TYPE
TraceType = 'Neu';
% TraceType = 'Cell';
% TraceType = 'pt7Neu';
% TraceType = '1Neu';
plottraces = 0;

% Normalization
wantzscore = 0;

% get STIM start
block{1} = '.\GAD1\190214\002\mgad_000_002.mat';
block{2} = '.\GADA\190128\002\mGADA_000_002.mat';
block{3} = '.\GADB\190128\004\mGADb_000_004.mat';
block{4} = '.\GADD\190128\002\mGADd_000_002.mat';
block{5} = '.\GADC\190130\002\mGADC_000_002.mat';


% get STIM ID
bl{1} = '.\GAD1\190214\002\Block-198.mat';
bl{2} = '.\GADA\190128\002\Block-200.mat';
bl{3} = '.\GADB\190128\004\Block-206.mat';
bl{4} = '.\GADD\190128\002\Block-202.mat';
bl{5} = '.\GADC\190130\002\Block-208.mat';
% bl = load('Block-208.mat');

% % Get cell data (MD)
% plane2{1} = '.\GAD1\190214\002\F_GAD1_190124_plane2.mat';
% plane2{2} = '.\GADA\190128\002\F_GADA_190128_plane2.mat';
% plane2{3} = '.\GADB\190128\004\F_GADB_190128_plane2.mat';
% plane2{4} = '.\GADD\190128\002\F_GADD_190128_plane2.mat';
% plane2{5} = '.\GADC\190130\002\F_GADC_190130_plane2.mat';

% Get cell data allpairs JA
plane2{1} = '.\GAD1\190214\002\suite2p\F_stim_GAD1_allplanes_JA.mat';
plane2{2} = '.\GADA\190128\002\tiffs\suite2p\F_stim_GADA_allplanes_JA.mat';
plane2{3} = '.\GADB\190128\004\tiffs\suite2p\F_stim_GADB_allplanes_JA.mat';
plane2{4} = '.\GADC\190130\002\tiffs\suite2p\F_stim_GADC_allplanes_JA.mat';
plane2{5} = '.\GADD\190128\002\tiffs\suite2p\F_stim_GADD_allplanes_JA.mat';


% Load electrode positions
zdiff = 20; % um
epos = load('.\epos.mat');

% begin analysis
nblock = length(block);

frames = cell(nblock,1);
AmplVals = cell(nblock,1);
ncell = NaN(nblock,1);
cells = cell(nblock,1);
Neur = cell(nblock,1);
zcells = cell(nblock,1);
zNeur = cell(nblock,1);
locs = cell(nblock,1);
rellocs = cell(nblock,1); % rel locs = location relative to electrode tip
meanimage = cell(nblock,1);

%%
for k = 3:5
    % get stim type info
    tmpblock = load(block{k});
    frm = tmpblock.info.frame;
    df = diff(frm);
    framestmp = [frm(df>100)+1];
    frames{k} = round(framestmp/3);
    
    % get stim Amplitude info
    tmpbl = load(bl{k});
    Ampl = tmpbl.params.Ampl;
    da = diff(Ampl);
    tmpAmplVals = Ampl(find(da>1)+1);
    AmplVals{k}=tmpAmplVals(1:end-1);
    
    tmpplane2 = load(plane2{k});
    iscell = logical(cat(1,tmpplane2.stat.iscell));
    tmplocs = double(cat(1,tmpplane2.stat.med));
    try 
        tmplocplanes = cat(1,tmpplane2.stat.iplane);
    catch
        tmpplanenum = str2num(tmpplane2.ops.fast_disk(end)) + 1;
        tmplocplanes = repmat(tmpplanenum,size(tmpplane2.stat,2),1);
    end
    locs{k} = tmplocs(iscell,:).* repmat(epos.micronppixel/epos.mag(k),size(tmplocs(iscell,:),1),1); 
    rellocs{k} = (locs{k} - epos.epos(k,1:2)) .* repmat(epos.micronppixel/epos.mag(k),size(locs{k},1),1); 
    %converted from pixels to microns
    locs{k}(:,3) = (tmplocplanes(iscell))*epos.micronpplane;
    rellocs{k}(:,3)=(tmplocplanes(iscell) - epos.epos(k,3))*epos.micronpplane;
    % 20micron apart between each plane
    if strcmp(TraceType,'Neu')
        cells{k} = tmpplane2.FcellNeu{1}(iscell,:);
    elseif strcmp(TraceType,'Cell')
        cells{k} = tmpplane2.Fcell{1}(iscell,:);
    elseif strcmp(TraceType,'pt7Neu')
        cells{k} = tmpplane2.Fcell{1}(iscell,:)-0.7*tmpplane2.FcellNeu{1}(iscell,:);
    elseif strcmp(TraceType,'1Neu')
        cells{k} = tmpplane2.Fcell{1}(iscell,:)-tmpplane2.FcellNeu{1}(iscell,:);
    end
    Neur{k} = tmpplane2.FcellNeu{1}(iscell,:);
    ncell(k) = size(cells{k},1);

    if wantzscore
        zcells{k} = zscore(cells{k}')';
        zNeur{k} = zscore(Neur{k}')';
    else
        zcells{k} = cells{k};
        zNeur{k} = Neur{k};
    end
try
    for iplane = 1:3
        meanimage{k}(:,:,iplane) = tmpplane2.ops{iplane}.mimg1(tmpplane2.ops{iplane}.yrange,tmpplane2.ops{iplane}.xrange);
    end
catch
    for iplane = 1:3
        tmpplane2.ops{iplane}.mimg1 = tmpplane2.ops{iplane}.meanImg;
        meanimage{k}(:,:,iplane) = tmpplane2.ops{iplane}.mimg1(tmpplane2.ops{iplane}.yrange,tmpplane2.ops{iplane}.xrange);
    end
end
    ncell(k) = size(cells{k},1);
end
%%
if plottraces == 1
    for m = 1:nblock
        figure(111+m), clf,hold on
        for k = 1:ncell(m)
            plot(cells{m}(k,:)+1000*k)
        end
        for g = 1:length(frames{m})
            plot([frames{m}(g) frames{m}(g)],[0 1000*ncell(m)],'g-')
            text(frames{m}(g),500,num2str(AmplVals{m}(g)))
        end
        title(['Expt. ',plane2{m}])
        set(gcf,'Position',[100 700 1000 1000])
        saveas(gcf,['Traces_',TraceType,plane2{m}(3:6),'.pdf'])
    end
end

%% Get cell resps
sigthresh = 1;

resps = cell(nblock,1);
respsorg = cell(nblock,1);
maxresp = cell(nblock,1);
maxrespidx = cell(nblock,1);
maxrespNeur = cell(nblock,1);
priorresp = cell(nblock,1);
sigact = cell(nblock,1);
sigactp = cell(nblock,1);
nsamp = 20;
amplid = unique(AmplVals{1});
nampl = length(amplid);
colorlist = colorbrewerRGB(nampl, 'qualitative');
c = turbo(nampl);

AvgResp = cell(nblock,1);
VarResp = cell(nblock,1);
NeurResp = cell(nblock,1);
RespEdges = [-10:1:10];
cc = turbo(length(RespEdges));

for m = 1:nblock
    ntrials = length(frames{m});
    ncell = size(zcells{m},1);
    resps{m} = NaN(ntrials,ncell*nsamp);
    respsorg{m} = NaN(ntrials,ncell,nsamp);
    maxresp{m} = NaN(ntrials,ncell);
    maxrespidx{m} = NaN(ntrials,ncell);
    maxrespNeur{m} = NaN(ntrials,ncell);
    priorresp{m} = NaN(ntrials,ncell);
    sigact{m} = zeros(ntrials,ncell);
    sigactp{m} = NaN(nampl,ncell);
    
    
    for k = 1:ntrials
        for g = 1:ncell
            resps{m}(k,((g-1)*nsamp+1):g*nsamp) = zcells{m}(g,frames{m}(k):(frames{m}(k)+nsamp-1));
            respsorg{m}(k,g,:) = zcells{m}(g,frames{m}(k):(frames{m}(k)+nsamp-1));
            maxresp{m}(k,g) = respsorg{m}(k,g,2);
            maxrespidx{m}(k,g) = 2;
            priorresp{m}(k,g) = mean(zcells{m}(g,frames{m}(k)-6:frames{m}(k)-1));
            tmpNeur = zNeur{m}(g,frames{m}(k)+1);
            priorNeur = mean(zNeur{m}(g,frames{m}(k)-6:frames{m}(k)-1));
            maxrespNeur{m}(k,g) = (tmpNeur-priorNeur)/priorNeur;
            % dF/F is a normalization, tells you how many baselines away
            % the cell changes to
            sigact{m}(k,g) = (maxresp{m}(k,g) - priorresp{m}(k,g))/priorresp{m}(k,g);
        end
    end
    
    if 0
        % compute p across ALL trials using signrank
        for g = 1:nampl
            for q = 1:ncell
                chtr = find(AmplVals{m}==amplid(g));
                priorcell = priorresp{m}(chtr,q);
                maxcell = maxresp{m}(chtr,q);
                sigactp{m}(g,q) = signrank(priorcell,maxcell);
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
    AvgResp{m} = NaN(nampl,ncell);
    VarResp{m} = NaN(nampl,ncell);
    NeurResp{m} = NaN(nampl,ncell);
    for k = 1:nampl
        AvgResp{m}(k,:) = nanmean(sigact{m}(AmplVals{m}==amplid(k),:));
        VarResp{m}(k,:) = nanvar(sigact{m}(AmplVals{m}==amplid(k),:));
        NeurResp{m}(k,:) = nanmean(maxrespNeur{m}(AmplVals{m}==amplid(k),:));
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
%     for k = 1:nampl
%         plot(AvgResp{m}(k,:),VarResp{m}(k,:),'.','color',colorlist(k,:))
%     end
%     legend(num2str(amplid));
%     xlabel('Avg Resp')
%     ylabel('Var Resp')
%     title(['Single cells in M',num2str(m)])
%     set(gcf,'Position',[10 10 1000 300])
    
%     figure(1000)
%     subplot(2,3,m), hold on
%     for k = 1:nampl
%         plot(AvgResp{m}(k,:),VarResp{m}(k,:),'.','color',colorlist(k,:))
%     end
%     legend(num2str(amplid),'location','west');
%     xlabel('Avg Resp')
%     ylabel('Var Resp')
%     title(['Single cells in M',num2str(m)])
%     set(gcf,'Position',[10 100 1000 500])
    
%     figure(1001)
%     subplot(2,3,m), hold on
%     for k = 1:nampl
%         plot(AvgResp{m}(k,:),sqrt(VarResp{m}(k,:))./AvgResp{m}(k,:),'.','color',colorlist(k,:))
%     end
%     legend(num2str(amplid),'location','west');
%     xlabel('Avg Resp')
%     ylabel('CV Resp')
%     title(['Single cells in M',num2str(m)])
%     set(gcf,'Position',[10 100 1000 500])
    

    
    if 0
        figure(10+m)
        for g = 1:nampl
            tmpavgresp = AvgResp{m}(g,:);
            tmpavgresp(tmpavgresp>10) = 10;
            tmpavgresp(tmpavgresp<-10) = -10;
            BINS = discretize(tmpavgresp,RespEdges);
            subplot(2,5,g), hold on
            imagesc(meanimage{m})
            for k = 1:ncell
                if AvgResp{m}(g,k)>= sigthresh
                    plot(locs{m}(k,2),locs{m}(k,1),'o',...
                        'markersize',tmpavgresp(k)+1,'color',cc(BINS(k),:),...
                        'markerfacecolor',cc(BINS(k),:))
                end
            end
            title(['M',num2str(m),', ',num2str(amplid(g)),' uA'])
            axis([61 724 0 506])
            plot(epos.epos(m,1),epos.epos(m,2),'rd','linewidth',2)
        end
        set(gcf,'Position',[100 500 1200 400])
        saveas(gcf,['Distr_of_Resps_',TraceType,'_Mouse_',num2str(m),'_sigthresh_',num2str(sigthresh),'.pdf'])
    end
    
end

% figure(1000), saveas(gcf,['Var_of_Resps_',TraceType,'.pdf'])
% figure(1001), saveas(gcf,['CV_of_Resps_',TraceType,'.pdf'])

%% Compute significant activation
amplid = unique(AmplVals{1});
nampl = length(amplid);
resact = cell(nblock,1);
sigthresh = 1;

% figure(222)
for m = 1:nblock % per mouse
    resact{m} = NaN(nampl,size(sigact{m},2));
    for k = 1:nampl % per amplitude
        % calculate fraction of sig. activations per cell
        ampltrials = find(AmplVals{m}==amplid(k));
        resact{m}(k,:) = mean(sigact{m}(ampltrials,:)>=sigthresh);
    end
    % sort cells by max response to 50 uA
    [sortval,sortidx] = sort(maxresp{m}(nampl,:));
%     subplot(2,3,m)
%     imagesc(resact{m}(:,sortidx),[0 1]),colorbar
%     title(num2str(m))
%     set(gca,'YTick',[1:nampl],'YTickLabel',amplid)
%     xlabel('Cell #')
%     ylabel('Stim Ampl (uA)')
end
% set(gcf,'Position',[10 300 1000 500])

sigsig = cell(nampl,1);
for m = 1:nblock
    s = sum(sigact{m}>=sigthresh,2)/size(sigact{m},2);
    for k = 1:nampl
        chtr = find(AmplVals{m}==amplid(k));
        sigsig{k} = [sigsig{k}(:); s(chtr)];
    end
end
% figure(223), hold on
meansig = NaN(nampl,1);
stdsig = NaN(nampl,1);
maxsig = NaN(nampl,1);
pctiles = [.05, 0.10, .25, .75, .90, 0.95];
pctvals = NaN(length(pctiles),nampl);
for k = 1:nampl
%     plot(amplid(k)*ones(length(sigsig{k}))+randn(length(sigsig{k}),1),sigsig{k},'.')
    meansig(k) = mean(sigsig{k});
    stdsig(k) = std(sigsig{k});
    y = quantile(sigsig{k},[pctiles]);
    pctvals(:,k) = y(:);
end
% plot(amplid,meansig,'ko-','linewidth',2)
for g = 1:length(pctiles)
%     plot(amplid,pctvals(g,:),'k-')
end
% title(['Significant Activation by Stim Amplitude: ',TraceType,', Thresh =',num2str(sigthresh)])
% xlabel('Fraction of cells Sig. Acticated')
% ylabel('Stim Ampl (uA)')
% saveas(gcf,['SigActivation_',TraceType,'_Thresh_',num2str(sigthresh),'.pdf'])
    
%% Plot sample responses for cells with large max responses
amplid = unique(AmplVals{1});
nampl = length(amplid);
ncelluse = 20;
for m = 1:nblock
    figure(m)
    for g = 1:ncelluse       
        for k = 1:nampl
            subplot(nampl,ncelluse,(k-1)*ncelluse + g)
            usetrials = find(AmplVals{m}==amplid(k));
            traces = squeeze(respsorg{m}(usetrials,g,:));
            plot(traces','k')
            title(['C',num2str(g),', ',num2str(amplid(k)),' uA'])
        end
        set(gcf,'Position',[10 10 1500 800])
    end
    saveas(gcf,['Sample_Evoked_Traces_',TraceType,'M',num2str(m),'.pdf'])
end
%% Plot cell resps as fxn of distance from electrode
amplid = unique(AmplVals{1});
nampl = length(amplid);
colorlist = colorbrewerRGB(nampl, 'qualitative');
markerlist = ['.','s','o','d','v','^'];
db = 50;
% distlist = [0:db:500, 600];
distlist = [0:db:1000];
distcent = distlist(2)/2 + distlist(1:end-1);
allCOM = NaN(nblock,3);
relcomvals = cell(nblock,1);
relcombins = cell(nblock,1);
cellbins = cell(nblock,1);
for m = 1:nblock
    dist = sqrt(rellocs{m}(:,1).^2 + rellocs{m}(:,2).^2+rellocs{m}(:,3).^2);%3D
    [N,EDGES,BIN] = histcounts(dist,distlist);
    display(max(dist))
    display(min(dist))
    cellbins{m} = BIN;
    binlist = unique(BIN);
    nbin = length(binlist);
    % Compute center of mass of max neuropil response during 50 uA stim
    COM = nansum([locs{m}(:,1).*NeurResp{m}(nampl,:)',locs{m}(:,2).*NeurResp{m}(nampl,:)',locs{m}(:,3).*NeurResp{m}(nampl,:)']/nansum(NeurResp{m}(nampl,:)));
%     COM = nansum([locs{m}(:,1).*AvgResp{m}(nampl,:)',locs{m}(:,2).*AvgResp{m}(nampl,:)']/nansum(AvgResp{m}(nampl,:)));
    allCOM(m,:) = COM;
    relcom = sqrt((locs{m}(:,1)-COM(1)).^2 + (locs{m}(:,2)-COM(2)).^2+(locs{m}(:,3)-COM(3)).^2);
    relcomvals{m} = relcom;
    [Nc,EDGESc,BINc] = histcounts(relcom,distlist);
    relcombins{m} = BINc;
end

% PLOT DATA
AmplResps = cat(2,AvgResp{:});
allbin = cell2mat(cellbins);
allbincom = cell2mat(relcombins);
binidxs = cell(nbin,1);
bincidx = cell(nbin,1);
for g = 1:nbin
    binidxs{g} = find(allbin==binlist(g));
    bincidx{g} = find(allbincom==binlist(g));
end
for k = 1:nampl
    figure(101), hold on
    subplot(5,2,k), hold on
    plot(distlist(allbin)+db/2+5*randn(size(distlist(allbin))),AmplResps(k,:),...
        '.','color',cc(k,:))
    
    %         mn = NaN(nbin,1);
    %         stdv = NaN(nbin,1);
    %         for g = 1:nbin
    %             mn(g) = nanmean(AvgResp{m}(k,BIN==binlist(g)));
    %             stdv(g) = nanstd(AvgResp{m}(k,BIN==binlist(g)));
    %         end
    %         errorbar(distlist(1:end-1)+db/2,mn,stdv,'ko-','linewidth',2)
    axis([0 600 0 10])
    if k == 1
        title([TraceType,': dF/F Max Resp. as FXN distance from electrode tip'])
    end
    
    figure(102), hold on
    subplot(5,2,k), hold on
    plot(distlist(allbincom)+db/2+5*randn(size(distlist(allbincom))),AmplResps(k,:),...
        '.','color',cc(k,:))

    axis([0 600 0 10])
    if k == 1
        title([TraceType,': dF/F Max Resp. as FXN distance from COM, 3 uA'])
    else
%         title([num2str(amplid(k)),' uA'])
    end
end % for k = 1:nampl

meancom = NaN(nampl,nbin);
stdcom = NaN(nampl,nbin);
meand = NaN(nampl,nbin);
stdd = NaN(nampl,nbin);
for k = 1:nampl
    for g = 1:nbin
%         errorbar(distlist(g)+db/2+k,nanmean(AmplResps(k,binidxs{g})),...
%             nanstd(AmplResps(k,binidxs{g})),'o-','linewidth',2,'color',colorlist(k,:))
        meand(k,g) = nanmean(AmplResps(k,binidxs{g}));
        stdd(k,g) = nanstd(AmplResps(k,binidxs{g}));
    end
    figure(101)
    subplot(5,2,10), hold on
%     errorbar(distlist(1:end-1) + db/2, meand(k,:), stdd(k,:),'o-','linewidth',2,'color',colorlist(k,:))
  
    plot(distcent(binlist), meand(k,:),'o-','linewidth',2,'color',colorlist(k,:))
    for g = 1:nbin
%         errorbar(distlist(g)+db/2+k*2,nanmean(AmplResps(k,bincidx{g})),...
%             nanstd(AmplResps(k,bincidx{g})),'o-','linewidth',2,'color',colorlist(k,:))
        meancom(k,g) = nanmean(AmplResps(k,bincidx{g}));
        stdcom(k,g) = nanstd(AmplResps(k,bincidx{g}));
    end
    figure(102)
    subplot(5,2,10), hold on
%     errorbar(distlist(1:end-1) + db/2, meancom(k,:), stdcom(k,:),'o-','linewidth',2,'color',colorlist(k,:))
    plot(distcent(binlist), meancom(k,:),'o-','linewidth',2,'color',colorlist(k,:))
end
figure(101)
xlabel('pixels (um)')
ylabel('Avg. change in dF/F with Stim')
set(gcf,'Position',[10 10 1000 700])
saveas(gcf,['Distr_of_Avg_Resps_From_Electrode_Tip_',TraceType,'.pdf'])
figure(102)
xlabel('pixels (um)')
ylabel('Avg. change in dF/F with Stim')
set(gcf,'Position',[10 10 1000 700])
saveas(gcf,['Distr_of_Avg_Resps_From_COM_',TraceType,'.pdf'])
%%
figure(103)
colorlist = colorbrewerRGB(nampl, 'qualitative');
for k = 1:nampl
%     subplot(2,2,1), hold on
%      errorbar(distlist(1:end-1) + db/2, meand(k,:), stdd(k,:),'o-','linewidth',2,'color',colorlist(k,:))
    subplot(2,1,1), hold on
    errorbar(distcent(binlist), meancom(k,:), stdcom(k,:),'o-','linewidth',2,'color',colorlist(k,:))
end
xlabel('Cell Distance')
ylabel('Average Evoked Response')
title('Cell activation as a function of Stimulation Amplitude and Cell Location COM')
legend(num2str(amplid))
cc = turbo(nbin);
for g = 1:nbin
%     subplot(2,2,3), hold on
%     errorbar(amplid, meand(:,g), stdd(:,g),'o-','linewidth',2,'color',cc(g,:))
    subplot(2,1,2), hold on
    errorbar(amplid, meancom(:,g), stdcom(:,g),'o-','linewidth',2,'color',cc(g,:))
end
xlabel('Stimulation Amplitude')
ylabel('Average Evoked Response')
legend(num2str(distlist(:)+db/2))
title('Cell activation as a function of Stimulation Amplitude and Cell Location COM')
set(gcf,'Position',[10 10 500 700])
saveas(gcf,['Distr_of_Resps_Ampl_Dist_COM_',TraceType,'.pdf'])

%% Plot correlation between cells and Neuropil
% PER MOUSE

for m = 1:nblock
    c = maxresp{m};
    n = maxrespNeur{m};
    figure(201)
    subplot(1,nblock,m)
    CR = corr(c,n);
    imagesc(CR), colorbar
    title('Corr cells and Neuropil')
    figure(202)
    subplot(1,nblock,m)
    CR = corr(c,c);
    imagesc(CR), colorbar
    title('Corr cells')
    figure(203)
    subplot(1,nblock,m)
    CR = corr(n,n);
    imagesc(CR), colorbar
    title('Corr Neuropili')
end

%% Look at variabililty in evoked response as a fxn of distance from the electrode tip

figure(102), hold on
for m = 1:nblock
    dist = sqrt(rellocs{m}(:,1).^2 + rellocs{m}(:,2).^2+ rellocs{m}(:,3).^2);
    for k = 1:nampl
        plot(dist,VarResp{m}(k,:),'.','color',colorlist(k,:),'marker',markerlist(m))
    end
end
legend(num2str(amplid))
title('Var Resp. as FXN distance from electrode tip')

figure(103), hold on
for m = 1:nblock
    dist = sqrt(rellocs{m}(:,1).^2 + rellocs{m}(:,2).^2+rellocs{m}(:,3).^2);
    for k = 1:nampl
        plot(dist,sqrt(VarResp{m}(k,:))./AvgResp{m}(k,:),'.','color',colorlist(k,:),'marker',markerlist(m))
    end
end
legend(num2str(amplid))
title('CV Resp. as FXN distance from electrode tip')


%%
% sort responses by stimulation amplitude
% for m = 1:nblock
%     [B,I] = sort(AmplVals{m});
%     if m == 1
%         figure(30), hold on
%         for k = 24
%                 subplot(5,5,k), hold on
%             subplot(1,2,1),hold on
%             for g = 1:ntrials
%                 plot(squeeze(respsorg{m}(I(g),k,:))+g*1000)
%                 hold on
%                 plot(maxrespidx{m}(I(g),k),maxresp{m}(I(g),k)+g*1000,'k*')
%             end
%             title(['sort ',num2str(k)])
%             subplot(1,2,2)
%             plot(AmplVals{k},...
%                 max(squeeze(respsorg{m}(:,k,:)),[],2)-squeeze(respsorg{m}(:,k,1)),...
%                 'ko')
%         end
%     end
% end
%

% DO PCA
rng(0)
colorlist = colorbrewerRGB(nblock, 'qualitative');
markerid = ['o','s','^','d','.'];
ncomponents = 10;
amplid = unique(AmplVals{1});
bvals = cell(nblock,1);

for m = 1:nblock
    [coeff, score, latent] = pca(maxresp{m},'NumComponents',ncomponents);
%     score = tsne(maxresp{m});
    ntrials = size(score,1);
    X = [ones(ntrials,1),score];
    
    
    % PLOT responses in PCA space
    amplid = unique(AmplVals{m});
    nampl = length(amplid);
    c = turbo(nampl);
    
    figure(1000), hold on
    subplot(2,3,m), hold on
    for k = 1:nampl
        plot(score(AmplVals{m}==amplid(k),1),score(AmplVals{m}==amplid(k),2),...
            '.','marker',markerid(m),'color',c(k,:),'markerfacecolor',c(k,:))
    end
    %     legend(num2str(amplid))
    xlabel('PC1')
    ylabel('PC2')
    title(['Mouse #',num2str(m)])
    
    subplot(2,3,nblock+1)
    for k = 1:nampl
        plot(score(AmplVals{m}==amplid(k),1),score(AmplVals{m}==amplid(k),2),'.','marker',markerid(m),...
            'color',c(k,:),'markerfacecolor',c(k,:))
    end
    xlabel('PC1')
    ylabel('PC2')
    title('All Superimposed')
    
    
    [sortvals,sortvalsid] = sort(AmplVals{m});
    ntrain = 60;
    Trainsub = randperm(length(AmplVals{m}));
    alltrials = 1:length(AmplVals{m});
    train = ismember(alltrials,Trainsub(1:ntrain));
    test = ~train;
    
    X2 = [ones(length(AmplVals{m}),1),score];
    y = AmplVals{m}(train==1);
    b = regress(y,X2(train==1,:));
    bvals{m} = b;
    yhat = X2(test(:)==1,:)*b;
    yt = AmplVals{m}(test==1);
    SStot = sum((yt - mean(yt)).^2);
    SSres = sum((yt - yhat).^2);
    Rsq = 1 - SSres/SStot;
    figure(2000), plot(AmplVals{m}(test==1)+0.3*m,yhat,'.','color',colorlist(m,:),...
        'marker',markerid(m),'markerfacecolor',colorlist(m,:))
    hold on, plot(AmplVals{m},AmplVals{m},'k-')
    
    %     hold on, plot(sortvals,yhat(sortvalsid),'-','color',colorlist(m,:),...
    %         'markerfacecolor',colorlist(m,:))
    xlabel('Stimulation Amplitude (uA)')
    ylabel('Evoked fluorescence (PCA coord.)')
    %     title('Regression on held-out data')
    title(sprintf('Regression on held-out data; Ntrain = %d, Ntest = %d',ntrain,length(AmplVals{m})-ntrain))
    text(10,30-2*m,['Rsq',num2str(m),'=',num2str(Rsq)])
    saveas(gcf,'LinearEvoked_correctd1NEU.pdf')
    
    figure(3000), hold on
    plot(score(:,1),AmplVals{m},'.','color',colorlist(m,:),...
        'marker',markerid(m),'markerfacecolor',colorlist(m,:))
    xlabel('PC 1')
    ylabel('Ampl (uA)')
%     plot3(score(:,1),score(:,2),X2*b,'-','color',colorlist(m,:))
%     text(10,30-2*m,['Rsq',num2str(m),'=',num2str(Rsq)])
    saveas(gcf,'Ampl_vs_PC1.pdf')
    
    
end
figure(1000),saveas(gcf,'PCA_traces.pdf')
figure(2000),saveas(gcf,'LinearEvoked_traces.pdf')


[bvals{1},bvals{2},bvals{3},bvals{4}]

%%
% DO PCA
rng(0)
colorlist = colorbrewerRGB(nblock, 'qualitative');
markerid = ['o','s','^','d','.'];
ncomponents = 10;
bigPCA = NaN(length(amplid)*9,ncomponents*nblock);
allAmplVals = NaN(length(amplid)*9*nblock,1);
for m = 1:nblock
    [coeff, score, latent] = pca(maxresp{m},'NumComponents',ncomponents);
    
    % PLOT responses in PCA space
    amplid = unique(AmplVals{m});
    nampl = length(amplid);
    c = turbo(nampl);
    
    % make sure you grab 9 of each trialtype
    usetrials = [];
    for k = 1:nampl
        tmp = find(AmplVals{m}==amplid(k));
        usetrials = [usetrials;tmp(1:9)];
    end
    usetrials = sort(usetrials);
    
    subPCA = score(usetrials,:);
    subAmplVals = AmplVals{m}(usetrials,:);
    
    % now sort Ampl Vals so they are in same order for each mouse
    [sortvals,sortvalsid] = sort(AmplVals{m});
    
    ntrials = length(subAmplVals);
    X = [ones(ntrials,1),subPCA];
    bigPCA(:,(m-1)*ncomponents+1:ncomponents*m) = subPCA;
    allAmplVals((m-1)*ntrials+1:ntrials*m)=subAmplVals;
    
    figure(1000), hold on
    subplot(2,3,m), hold on
    for k = 1:nampl
        for g = 1:sum(subAmplVals==amplid(k))
            tmpord = find(subAmplVals==amplid(k));
            plot(subPCA(tmpord(g),1),subPCA(tmpord(g),2),'o',...
                'color',c(k,:),'markerfacecolor',c(k,:),'markersize',1+g)
        end
    end
    %     legend(num2str(amplid))
    xlabel('PC1')
    ylabel('PC2')
    title(['Mouse #',num2str(m)])
    
    
    
    
    
    ntrain = 60;
    Trainsub = randperm(length(AmplVals{m}));
    alltrials = 1:length(AmplVals{m});
    train = ismember(alltrials,Trainsub(1:ntrain));
    test = ~train;
    
    X2 = [ones(length(AmplVals{m}),1),score];
    y = AmplVals{m}(train==1);
    b = regress(y,X2(train==1,:));
    yhat = X2(test(:)==1,:)*b;
    yt = AmplVals{m}(test==1);
    SStot = sum((yt - mean(yt)).^2);
    SSres = sum((yt - yhat).^2);
    Rsq = 1 - SSres/SStot;
    figure(2000), plot(AmplVals{m}(test==1)+0.3*m,yhat,'.','color',colorlist(m,:),...
        'marker',markerid(m),'markerfacecolor',colorlist(m,:))
    hold on, plot(AmplVals{m},AmplVals{m},'k-')
    
    %     hold on, plot(sortvals,yhat(sortvalsid),'-','color',colorlist(m,:),...
    %         'markerfacecolor',colorlist(m,:))
    xlabel('Stimulation Amplitude (uA)')
    ylabel('Evoked fluorescence (PCA coord.)')
    %     title('Regression on held-out data')
    title(sprintf('Regression on held-out data; Ntrain = %d, Ntest = %d',ntrain,length(AmplVals{m})-ntrain))
    text(10,30-2*m,['Rsq',num2str(m),'=',num2str(Rsq)])
    %     saveas(gcf,'LinearEvoked.pdf')
    
end

[coeff, finalscore, latent] = pca(bigPCA,'NumComponents',3);
figure(1000)
subplot(2,3,nblock+1), hold on
for k = 1:nampl
    tmpord = find(subAmplVals==amplid(k));
    for g = 1:sum(subAmplVals==amplid(k))
        plot(finalscore(tmpord(g),1),finalscore(tmpord(g),2),'o',...
            'color',c(k,:),'markerfacecolor',c(k,:),'markersize',1+g)
    end
end
xlabel('PC1')
ylabel('PC2')
title('All Superimposed')
legend(num2str(amplid))
saveas(gcf,'LinearEvoked.pdf')

% Plot how well you can do regression with ALL PCA
AllAmplVals = cell2mat(AmplVals);
[sortvals,sortvalsid] = sort(AllAmplVals);
ntrain = round(0.8*length(AllAmplVals));
Trainsub = randperm(length(AllAmplVals));
alltrials = 1:length(AllAmplVals);
train = ismember(alltrials,Trainsub(1:ntrain));
test = ~train;

X2 = [ones(length(AmplVals{m}),1),score];
y = AmplVals{m}(train==1);
b = regress(y,X2(train==1,:));
yhat = X2(test(:)==1,:)*b;
yt = AmplVals{m}(test==1);
SStot = sum((yt - mean(yt)).^2);
SSres = sum((yt - yhat).^2);
Rsq = 1 - SSres/SStot;
figure(1001)

%% plot max evoked response as function of time
amp = unique(AmplVals{1});
figure(222), clf, hold on
for g = 1:length(AmplVals)
    for k = 1:length(amp)
        tmp = find(AmplVals{g}==amp(k));
        x = amp(k)*ones(length(tmp),1);
        x = x + linspace(0,5,length(tmp))';
        allcell = maxresp{g}(tmp,:);
        
        plot(x,mean(allcell,2),'-')
        
        % do stats on first vs. last
        y1 = allcell(1,:);
        y2 = allcell(length(tmp),:);
        d = mean(y2) - mean(y1);
        %     p = signrank(y1,y2);
        %     [h,p] = ttest(y1,y2)
        %     text(double(x(1)),mean(allcell(1,:)),num2str(d))
        %     if p < 1E-3
        %         text(double(x(1)),mean(allcell(1,:))+0.2,'*')
        %     end
    end
end
title('Averaged across all cells')
xlabel('Stimulation Amplitude (uA)')
ylabel('Ampl. of evoked response')
set(gcf,'Position',[100 200 600 200])
saveas(gcf,'Evoked_as_fxn_time.pdf')


