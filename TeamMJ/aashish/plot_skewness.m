%% plot skewness
function plot_skewness(pre, post)
%% 
% for m = 1:5
%     for i = 1:numel(pre{m}.stat);
%         pre_skew{m}(i) = pre{m}.stat{i}.skew;
%     end
%     for j = 1:numel(post{m}.stat);
%         post_skew{m}(j) = post{m}.stat{j}.skew;
%     end
%     figure; hold on;
%     subplot(2,3,m)
%     boxplot([pre_skew{m}', post_skew{m}']);
%     plot([1;2], [pre_skew{m}', post_skew{m}'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5);
%     title('Skewness, pre and post, all mice')
%     xlim([.5 2.5]);
%     xticks(1:2);
%     xticklabels({'pre-stimulation-skew','post-stimulation-skew'})
% end
%
    for m = 1:5
        for i = 1:numel(pre{m}.stat);
            pre_skew{m}(i) = pre{m}.stat{i}.skew;
        end
        for j = 1:numel(post{m}.stat);
            post_skew{m}(j) = post{m}.stat{j}.skew;
        end
    end
    all_pre_skew = [pre_skew{1},pre_skew{2},pre_skew{3},pre_skew{4},pre_skew{5}];
    all_post_skew = [post_skew{1},post_skew{2},post_skew{3},post_skew{4},post_skew{5}];
    
    figure; hold on;
    boxplot([all_pre_skew', all_post_skew']);
    plot([1;2], [all_pre_skew', all_post_skew'],'x-', 'Color',[.5 .5 .5 .3], 'LineWidth',.5);
    title('Skewness, pre and post, all mice')
    xlim([.5 2.5]);
    xticks(1:2);
    xticklabels({'pre-stimulation-skew','post-stimulation-skew'})