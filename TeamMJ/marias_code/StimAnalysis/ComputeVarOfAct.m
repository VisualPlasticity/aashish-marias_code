function [] = ComputeVarOfAct(data,TraceType,savedatachoice)
%function [] = ComputeVarOfAct(data,TraceType,savedatachoice)
if nargin < 3
    savedatachoice = 0;
end

% Variance of Responses
amplid = data.amplid;
nampl = length(amplid);
nblock = length(data.AmplVals);


% compute locations relative to COM of neuropil activation at 50 uA
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


allAvgResp = cell(nampl,nbin);
allVarResp = cell(nampl,nbin);
for m = 1:nblock
    for g = 1:nbin
        subcell = find(relcombins{m}==binlist(g));
        for k = 1:nampl
            allVarResp{k,g} = [allVarResp{k,g},data.VarResp{m}(k,subcell)];
            allAvgResp{k,g} = [allAvgResp{k,g},data.AvgResp{m}(k,subcell)];
        end
    end
end
%colorlist = colorbrewerRGB(nampl, 'qualitative');
colorlist = turbo(data.nampl);
figure(1234), hold on
set(gcf,'Position',[100 100 300 800])
for g = 1:nbin
    for k = 1:nampl
        subplot(4,3,k),hold on
        plot(allAvgResp{k,g},allVarResp{k,g},'.','color',colorlist(k,:))
    end
end

