function [block,bl,plane2,plottraces,wantzscore,sigthresh] = db(TraceType)

if nargin < 1
    % TRACE TYPE
    TraceType = 'Neu';
    % TraceType = 'Cell';
    % TraceType = 'pt7Neu';
    % TraceType = '1Neu';
end

sigthresh = 1;

plottraces = 0;

% Normalization
wantzscore = 0;

% get STIM start
block{1} = './GAD1/190214/002/mgad_000_002.mat';
block{2} = './GADA/190128/002/mGADA_000_002.mat';
block{3} = './GADB/190128/004/mGADb_000_004.mat';
block{4} = './GADD/190128/002/mGADd_000_002.mat';
block{5} = './GADC/190130/002/mGADC_000_002.mat';

% get STIM ID
bl{1} = './GAD1/190214/002/Block-198.mat';
bl{2} = './GADA/190128/002/Block-200.mat';
bl{3} = './GADB/190128/004/Block-206.mat';
bl{4} = './GADD/190128/002/Block-202.mat';
bl{5} = './GADC/190130/002/Block-208.mat';

% Get cell data
plane2{1} = '.\GAD1\190214\002\suite2p\F_stim_GAD1_allplanes_JA.mat';
plane2{2} = '.\GADA\190128\002\tiffs\suite2p\F_stim_GADA_allplanes_JA.mat'
plane2{3} = '.\GADB\190128\004\tiffs\suite2p\F_stim_GADB_allplanes_JA.mat'
plane2{4} = '.\GADD\190128\002\tiffs\suite2p\F_stim_GADD_allplanes_JA.mat'
plane2{5} = './GADC/190130/002/tiffs/suite2p/F_stim_GADC_allplanes_JA.mat';
