% Name: pre_post_analysis
% Author: Aashish Khimasia
% Date: 08/08/21
%% Load data 
[epos, se] = create_se()
%% Create structures 'pre', 'post'
[pre, post, neur_coeff, redcell, npair, nframe] = create_pre_post(se);
%% Set movie frames
[movie_frames_pre, movie_frames_post, movie_sect_length] = set_movie_frames(pre, post, nframe);
%% Plot PCA per mouse
[score] = plot_pca(pre, post, movie_frames_pre, movie_frames_post, movie_sect_length);
%% Plot PCA corr pre vs post
plot_pca_corr(pre, post, movie_sect_length, score);
%% Plot skewness 
plot_skewness(pre, post);
%% Plot skew corr pre vs post
plot_skewness_popn_corr(pre, post, movie_frames_pre, movie_frames_post, redcell, npair);
%% Plot skew rc 
plot_skewness_rc(pre, post, redcell);
%% Plot skewness change vs distance 
[rellocs] = plot_skewness_dist_prepost(pre, post, npair, epos, redcell);