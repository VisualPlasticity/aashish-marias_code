%% calculate pca correlation between pre and post
function plot_pca_corr(pre, post, movie_sect_length, score)

for m = 1:5
    for i = 1:size(pre{m}.F,1)
        R{m} = corrcoef(score{m}(1:movie_sect_length{m},i), score{m}((movie_sect_length{m} + 1):movie_sect_length{m} * 2,i));
        pca_corr_coeff{m}(i) = R{m}(1,2);
    end

    figure;hold on;
    for i = 1:size(pre{m}.F,1)
        plot(i, pca_corr_coeff{m}(i), 'x', 'color', 'r')
    end
    title([pre{m}.ops.mouse_name ': corr pre vs post for each PC']);
    xlabel('PC');
    ylabel('correlation_coefficient');
end