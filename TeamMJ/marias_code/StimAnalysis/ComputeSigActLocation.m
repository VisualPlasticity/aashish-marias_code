function [] = ComputeSigActLocation(data,TraceType,savedatachoice)
% function [] = ComputeSigAct(data,TraceType,savedatachoice)

if nargin < 3
    savedatachoice = 0;
end
% Significant Activation

amplid = data.amplid;
nampl = length(amplid);
nblock = length(data.AmplVals);
sigthresh = data.sigthresh;


% compute COM locs
db = 50;
distlist = [0:db:500, 600];
relcomvals = cell(nblock,1);
relcombins = cell(nblock,1);
cellbins = cell(nblock,1);
for m = 1:nblock
    dist = sqrt(data.rellocs{m}(:,1).^2 + data.rellocs{m}(:,2).^2);
    [N,EDGES,BIN] = histcounts(dist,distlist);
    cellbins{m} = BIN;
    binlist = unique(BIN);
    nbin = length(binlist);
    % Compute center of mass of max neuropil response during 50 uA stim
    COM = nansum([data.locs{m}(:,1).*data.NeurResp{m}(nampl,:)',...
        data.locs{m}(:,2).*data.NeurResp{m}(nampl,:)']/nansum(data.NeurResp{m}(nampl,:)));
    relcom = sqrt((data.locs{m}(:,1)-COM(1)).^2 + (data.locs{m}(:,2)-COM(2)).^2);
    relcomvals{m} = relcom;
    [Nc,EDGESc,BINc] = histcounts(relcom,distlist);
    relcombins{m} = BINc;
end

% compute significant activations in each bin

sigsig = cell(nampl,nbin);
for m = 1:nblock
    for g = 1:nbin
        subcell = find(relcombins{m}==binlist(g));
        s = sum(data.sigact{m}(:,subcell)>=sigthresh,2)/length(subcell);
        for k = 1:nampl
            sigsig{k,g} = [sigsig{k,g}(:); s(data.AmplVals{m}==amplid(k))];
        end
    end
end

sigsigavg = NaN(nampl,nbin);
for g = 1:nbin
    for k = 1:nampl
        sigsigavg(k,g) = mean(sigsig{k,g});
    end
end
figure(), imagesc(sigsigavg), colorbar
set(gca,'XTick',[1:nbin],'XTickLabel',[distlist(1:end-1)+db/2])
set(gca,'YTick',[1:nampl],'YTickLabel',amplid)
xlabel('Distance from COM of Evoked Neuropil Activity')
ylabel('Stimulation Amplitude (uA)')
title(['Relationship between Stim Amplitude and Significant Evoked Activity: ',...
    TraceType,', Sig Thresh =',num2str(sigthresh)])
if savedatachoice
    saveas(gcf,['SigActivation_byLoc_',TraceType,'_Thresh_',num2str(sigthresh),'.pdf'])
end