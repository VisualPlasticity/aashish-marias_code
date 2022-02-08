%% create pre and post
function [pre, post, neur_coeff, redcell, npair, nframe] = create_pre_post(se)

neur_coeff = 0.7;
for m = 1:5 {m}
    pre{m}.F =se{m}.F_pre - neur_coeff*se{m}.Fneu_pre;
    post{m}.F =se{m}.F_post - neur_coeff*se{m}.Fneu_post;
    pre{m}.stat = se{m}.stat_pre; 
    post{m}.stat = se{m}.stat_post; 

    redcell{m} = se{m}.redcell;
    pre{m}.velocity = se{m}.velocity_pre;
    post{m}.velocity = se{m}.velocity_post;

    % pre.ops and post.ops created only using ops from plane1
    pre{m}.ops = se{m}.ops_pre{1};
    post{m}.ops = se{m}.ops_post{1};

    npair{m} = size(pre{m}.F,1);
    nframe{m} = size(pre{m}.F,2); 
end