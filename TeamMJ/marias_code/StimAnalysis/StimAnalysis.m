% ANALYZE STIM-EVOKED DATA
clear all
%close all
clc

% Trace Type
% TraceType = 'Neu';
% TraceType = 'Cell';
TraceType = 'pt7Neu';
% TraceType = '1Neu';
savedatachoice = 0;

% initialize params
[block,bl,plane2,plottraces,wantzscore,sigthresh] = db(TraceType);

%% Run initial analysis

data = InitAnalysis(TraceType,bl,block,plane2,sigthresh,plottraces,wantzscore);


%% Plot Average Response as a function of Stim Amplitude, include single-cell data (NO LOCATION)
analyze_dist_resps
%% Plot Significant Activation as a function of Stim Amplitude, include single-trial data (no location)
% in this plot, each point signifies a single trial
% data is complied for all 5 mice together
ComputeSigAct(data,TraceType,savedatachoice);

%% Plot Significant Activation as a function of Stim Amplitude AND Location
% in this plot, each point signifies a single trial
% data is complied for all 5 mice together
ComputeSigActLocation(data,TraceType,1);

%PlotSampleTracesSigActLocation(data,TraceType,1);
%% Compute Variance of Responses, look as function of Location
ComputeVarOfAct(data,TraceType,savedatachoice);

