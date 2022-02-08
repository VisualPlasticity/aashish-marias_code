% ANALYZE STIM-EVOKED DATA

% get STIM start
block = load('mgad_000_002.mat');
frm = block.info.frame;
df = diff(frm);
frames = [frm(1);frm(df>100)+1];
frames = round(frames/3);

% get STIM ID
bl = load('Block-198.mat');
Ampl = bl.params.Ampl;
da = diff(Ampl);
AmplVals = Ampl(find(da>1)+1);

% Get cell data
plane2 = load('F_GAD1_190124_plane2.mat');
iscell = cat(1,plane2.stat.iscell);
cells = plane2.Fcell{1}(iscell,:);

zcells = zscore(cells')';

% Get cell resps
nsamp = 20;
ntrials = length(frames);
ncell = size(cells,1);
resps = NaN(ntrials,ncell*nsamp);
respsorg = NaN(ntrials,ncell,nsamp);
% num trials, cell number, num samples
for k = 1:ntrials
    for g = 1:ncell
        resps(k,((g-1)*nsamp+1):g*nsamp) = cells(g,frames(k):(frames(k)+nsamp-1));
        %resps(trial number, 20 frames in order (0:20, 21:40...)) =         %cell number 
        respsorg(k,g,:) = cells(g,frames(k):(frames(k)+nsamp-1));
    end
end

figure(1), imagesc(resps)

figure(2), 
for k = 1:64 %here k is a cell number from 1 to 64
    subplot(8,8,k)
    imagesc(squeeze(respsorg(:,k,:)))

    %displays an image of the squeezed 2D array of each trial (y) vs 20
    %frames (x) from that trial for the cell number k
    title(['Cell ' num2str(k)])    
    xlabel('frame number after stimulation')
    ylabel('trial number')
end

% sort responses by stimulation amplitude
[B,I] = sort(AmplVals);
figure(3), hold on
for k = 24 %cell number 24
%     subplot(5,5,k), hold on
subplot(1,2,1),hold on
    for g = 1:ntrials
        plot(squeeze(respsorg(I(g),k,:))+g*500)
    end
    title(['Response for cell ',num2str(k), ' for 20 frames after stimuation across 90 trials'])
    xlabel('frame number after stimulation')
    ylabel('trial number')
subplot(1,2,2)
    plot(AmplVals,...
        max(squeeze(respsorg(:,k,:)),[],2)-squeeze(respsorg(:,k,1)),...
        'ko')
    %returns max response value from the 2D array of ntrials by nsamples formed from removing dimensions of
    %length 1, from the 2nd dimension MINUS the response from all trials at
    %frame 1 after stimulation
    title(['Stimulation amplitude against response amplitude for cell ',num2str(k)])
    xlabel('stimulation amplitude')
    ylabel('response amplitude (max response - baseline)')
   %squeezes a 3D array into a 2D array of 90 trials by 20 samples (frames) 
end


% DO PCA
[coeff, score, latent] = pca(squeeze(respsorg(:,1,:)));
X = [ones(ntrials,1),score(:,1:10)];

% PLOT responses in PCA space
amplid = unique(AmplVals);
nampl = length(amplid);
c = turbo(nampl);

figure(), hold on
for k = 1:nampl
    title('PCA');
    xlabel('PC1');
    ylabel('PC2');
    plot(score(AmplVals==amplid(k),1),score(AmplVals==amplid(k),2),...
        'x','color', c(k,:))
    grid on;
    legend('Stimulation Amplitude 1','Stimulation Amplitude 2','Stimulation Amplitude 3','Stimulation Amplitude 4','Stimulation Amplitude 5','Stimulation Amplitude 6','Stimulation Amplitude 7','Stimulation Amplitude 8','Stimulation Amplitude 9')
end

