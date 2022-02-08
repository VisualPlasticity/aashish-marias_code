function [] = ComputeSigAct(data,TraceType,savedatachoice)
%function [] = ComputeSigAct(data,TraceType,savedatachoice)
if nargin < 3
    savedatachoice = 0;
end

% Significant Activation
amplid = data.amplid;
nampl = length(amplid);
nblock = length(data.AmplVals);
resact = cell(nblock,1);
sigthresh = data.sigthresh;


% figure(222)
for m = 1:nblock % per mouse
    resact{m} = NaN(nampl,size(data.sigact{m},2));
    for k = 1:nampl % per amplitude
        % calculate fraction of sig. activations per cell
        ampltrials = find(data.AmplVals{m}==amplid(k));
        resact{m}(k,:) = mean(data.sigact{m}(ampltrials,:)>=sigthresh);
    end
    % sort cells by max response to 50 uA
    [sortval,sortidx] = sort(data.sigact{m}(nampl,:));
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
    s = sum(data.sigact{m}>=sigthresh,2)/size(data.sigact{m},2);
    for k = 1:nampl
        chtr = find(data.AmplVals{m}==amplid(k));
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
    meansig(k) = median(sigsig{k});
    stdsig(k) = std(sigsig{k});
    y = quantile(sigsig{k},[pctiles]);
    pctvals(:,k) = y(:);
end
plot(amplid,meansig,'ko-','linewidth',2)
for g = 1:length(pctiles)
    plot(amplid,pctvals(g,:),'k-')
end
title(['Median Significant Activation by Stim Amplitude Per Trial: ',TraceType,', Thresh =',num2str(sigthresh)])
ylabel('Fraction of cells Sig. Acticated')
xlabel('Stim Ampl (uA)')
if savedatachoice
    saveas(gcf,['SigActivation_',TraceType,'_Thresh_',num2str(sigthresh),'.pdf'])
end