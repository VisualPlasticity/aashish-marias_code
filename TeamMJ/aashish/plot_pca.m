%% pca analysis
function [score] = plot_pca(pre, post, movie_frames_pre, movie_frames_post, movie_sect_length)

for m = 1:5
    %normalising the data
    normalised_pre = zscore(pre{m}.F(:,movie_frames_pre{m}), 0, 2);
    normalised_post = zscore(post{m}.F(:,movie_frames_post{m}), 0, 2);
    normalised = cat(2,normalised_pre,normalised_post);
    
    %run pca analysis
    [coeff{m}, score{m}, latent{m}, ~, explained{m}] = pca(normalised');

    %eigenvectors, scores(remapping of data), eigenvalues, ~, percentage of
    %variance explained by the PC

    c = turbo(movie_sect_length{m});

    %BOTH SESSIONS
    figure;hold on;
    title('PCA scores');
    %scores
    subplot(1,2,1); hold on;
    for i = 1:movie_sect_length{m}
        plot3(score{m}(i,1), score{m}(i,2), score{m}(i,3), '.', 'color', c(i,:))

    end
    figaxes(1) = gca;
    figaxes(1).XLabel.String = 'PC1';
    figaxes(1).YLabel.String = 'PC2';
    figaxes(1).ZLabel.String = 'PC3';
    figaxes(1).Title.String = [pre{m}.ops.mouse_name ' pre-stim'] % thr_status];
    view(135, 20);
    grid on;
    subplot(1,2,2); hold on;
    for i = (movie_sect_length{m} + 1):(movie_sect_length{m} + movie_sect_length{m})
        plot3(score{m}(i,1), score{m}(i,2), score{m}(i,3), '.', 'color', c(i - movie_sect_length{m},:))

    end
    figaxes(2) = gca;
    figaxes(2).XLabel.String = 'PC1';
    figaxes(2).YLabel.String = 'PC2';
    figaxes(2).ZLabel.String = 'PC3';
    figaxes(2).Title.String = [post{m}.ops.mouse_name ' post-stim'] % thr_status];
    view(135, 20);
    axes_array_1(1) = figaxes(1).ZLim(1);
    axes_array_1(2) = figaxes(2).ZLim(1);
    axes_array_2(1) = figaxes(1).ZLim(2);
    axes_array_2(2) = figaxes(2).ZLim(2);
    axes_array_sorted_1 = sort(axes_array_1);
    axes_array_sorted_2 = sort(axes_array_2);
    Z_Lim(1) = axes_array_sorted_1(1);
    Z_Lim(2) = axes_array_sorted_2(2);
    figaxes(1).ZLim = Z_Lim;
    figaxes(2).ZLim = Z_Lim;
    grid on;

    linkaxes(figaxes);
    color_bar = colorbar;
    color_bar.Ticks = [0 1];
    color_bar.TickLabels = {'1', num2str(movie_sect_length{m})};
    color_bar.Label.String = 'Movie Frame Num';
    colormap turbo;

end