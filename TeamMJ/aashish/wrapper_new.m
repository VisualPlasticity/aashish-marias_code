epos = load('epos.mat');
se{1} = load('F_GAD1_allplanes_JA.mat');
se{2} = load('F_GADA_allplanes_JA.mat');
se{3} = load('F_GADB_allplanes_JA.mat');
se{4} = load('F_GADC_allplanes_JA.mat');
se{5} = load('F_GADD_allplanes_JA.mat');
clc
%% Create pre post
[pre, post, neur_coeff, redcell, npair, nframe] = create_pre_post(se)
%% Set movie frames
[movie_frames_pre, movie_frames_post, movie_sect_length] = set_movie_frames(pre, post, nframe)
%% Plot PCA per mouse
[score] = plot_pca(pre, post, movie_frames_pre, movie_frames_post, movie_sect_length)
%% Plot PCA corr pre vs post
plot_pca_corr(pre, post, movie_sect_length, score)
%% Plot skewness 
plot_skewness(pre, post)
%% Plot skew corr pre vs post
plot_skewness_popn_corr(pre, post, movie_frames_pre, movie_frames_post, redcell, npair)
%% Plot skew rc 
plot_skewness_rc(pre, post, redcell)
%% Plot skewness change vs distance 
[rellocs] = plot_skewness_dist_prepost(pre, post, npair, epos, redcell)