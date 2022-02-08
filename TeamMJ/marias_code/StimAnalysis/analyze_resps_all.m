% ANALYZE STIM-EVOKED DATA
clear all
%close all
clc

% TRACE TYPE
% TraceType = 'Neu';
% TraceType = 'Cell';
% TraceType = 'pt7Neu';
TraceType = '1Neu';
plottraces = 0;


% get STIM start
block{1} = './GAD1/190214/002/mgad_000_002.mat';
block{2} = './GADA/190128/002/mGADA_000_002.mat';
block{3} = './GADB/190128/004/mGADb_000_004.mat';
block{4} = './GADD/190128/002/mGADd_000_002.mat';
% block = load('mGADC_000_002.mat');

% get STIM ID
bl{1} = './GAD1/190214/002/Block-198.mat';
bl{2} = './GADA/190128/002/Block-200.mat';
bl{3} = './GADB/190128/004/Block-206.mat';
bl{4} = './GADD/190128/002/Block-202.mat';
% bl = load('Block-208.mat');

% Get cell data
plane2{1} = './GAD1/190214/002/F_GAD1_190124_plane2.mat';
plane2{2} = './GADA/190128/002/F_GADA_190128_plane2.mat';
plane2{3} = './GADB/190128/004/F_GADb_190128_plane2.mat';
plane2{4} = './GADD/190128/002/F_GADD_190128_plane2.mat';
% plane2 = load('F_GADC_190130_plane1.mat');


% begin analysis
nblock = length(block);

frames = cell(nblock,1);
AmplVals = cell(nblock,1);
ncell = NaN(nblock,1);
cells = cell(nblock,1);
zcells = cell(nblock,1);
for k = 1:nblock
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
    iscell = cat(1,tmpplane2.stat.iscell);
    if strcmp(TraceType,'Neu')
        cells{k} = tmpplane2.FcellNeu{1}(iscell,:);
    elseif strcmp(TraceType,'Cell')
        cells{k} = tmpplane2.Fcell{1}(iscell,:);
    elseif strcmp(TraceType,'pt7Neu')
        cells{k} = tmpplane2.Fcell{1}(iscell,:)-0.7*tmpplane2.FcellNeu{1}(iscell,:);
    elseif strcmp(TraceType,'1Neu')
        cells{k} = tmpplane2.Fcell{1}(iscell,:)-tmpplane2.FcellNeu{1}(iscell,:);
    end
    
    ncell(k) = size(cells{k},1);
    zcells{k} = zscore(cells{k}')';
end

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
resps = cell(nblock,1);
respsorg = cell(nblock,1);
maxresp = cell(nblock,1);
maxrespidx = cell(nblock,1);
priorresp = cell(nblock,1);
sigact = cell(nblock,1);
sigactp = cell(nblock,1);
nsamp = 20;
for m = 1:nblock
    ntrials = length(frames{m});
    ncell = size(cells{m},1);
    resps{m} = NaN(ntrials,ncell*nsamp);
    respsorg{m} = NaN(ntrials,ncell,nsamp);
    maxresp{m} = NaN(ntrials,ncell);
    maxrespidx{m} = NaN(ntrials,ncell);
    priorresp{m} = NaN(ntrials,ncell);
    sigact{m} = zeros(ntrials,ncell);
    sigactp{m} = NaN(ntrials,ncell);
    for k = 1:ntrials
        for g = 1:ncell
            resps{m}(k,((g-1)*nsamp+1):g*nsamp) = zcells{m}(g,frames{m}(k):(frames{m}(k)+nsamp-1));
            respsorg{m}(k,g,:) = zcells{m}(g,frames{m}(k):(frames{m}(k)+nsamp-1));
            [maxresp{m}(k,g),maxrespidx{m}(k,g)] = max(squeeze(respsorg{m}(k,g,:)));
            priorresp{m}(k,g) = mean(zcells{m}(g,frames{m}(k)-6:frames{m}(k)-1));
            sigact{m}(k,g) = (maxresp{m}(k,g) - priorresp{m}(k,g))/priorresp{m}(k,g);
        end
        % compute p 
%         for g = 1:ncell
%             var
%         end
    end
    if 0
        figure(m), imagesc(resps{m})
        xlabel('Time')
        ylabel('Cell #')
    end
end

%% Compute significant activation
amplid = unique(AmplVals{m});
nampl = length(amplid);
resact = cell(nblock,1);
sigthresh = 3;

figure(222)
for m = 1:nblock % per mouse
    resact{m} = NaN(nampl,size(sigact{m},2));
    for k = 1:nampl % per amplitude
        % calculate fraction of sig. activations per cell
        ampltrials = find(AmplVals{m}==amplid(k));
        resact{m}(k,:) = mean(sigact{m}(ampltrials,:)>=sigthresh);
    end
    % sort cells by max response to 50 uA
    [sortval,sortidx] = sort(maxresp{m}(nampl,:));
    subplot(2,3,m)
    imagesc(resact{m}(:,sortidx),[0 1]),colorbar
    title(num2str(m))
    set(gca,'YTick',[1:nampl],'YTickLabel',amplid)
    xlabel('Cell #')
    ylabel('Stim Ampl (uA)')
end
set(gcf,'Position',[10 300 1000 500])

sigsig = cell(nampl,1);
for m = 1:nblock
    s = sum(sigact{m}>=sigthresh,2)/size(sigact{m},2);
    for k = 1:nampl
        chtr = find(AmplVals{m}==amplid(k));
        sigsig{k} = [sigsig{k}(:); s(chtr)];
    end
end
figure(223), hold on
meansig = NaN(nampl,1);
stdsig = NaN(nampl,1);
maxsig = NaN(nampl,1);
pctiles = [.05, 0.10, .25, .75, .90, 0.95];
pctvals = NaN(length(pctiles),nampl);
for k = 1:nampl
    plot(amplid(k)*ones(length(sigsig{k}))+randn(length(sigsig{k}),1),sigsig{k},'.')
    meansig(k) = mean(sigsig{k});
    stdsig(k) = std(sigsig{k});
    y = quantile(sigsig{k},[pctiles]);
    pctvals(:,k) = y(:);
end
plot(amplid,meansig,'ko-','linewidth',2)
for g = 1:length(pctiles)
    plot(amplid,pctvals(g,:),'k-')
end

%% Plot sample responses for cells with large max responses
amplid = unique(AmplVals{1});
nampl = length(amplid);
ncelluse = 20;
for m = 1:nblock
    figure(1)
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
end
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