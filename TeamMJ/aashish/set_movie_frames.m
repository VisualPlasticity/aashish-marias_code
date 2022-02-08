%% set movie frames
function [movie_frames_pre, movie_frames_post, movie_sect_length] = set_movie_frames(pre, post, nframe) 

for m = 1:5
    movie_frames_pre{m} = [];
    movie_frames_post{m} = [];

    all_start_frame_pre{m} = ceil((pre{m}.ops.start_frame-1)/3+1);
    all_start_frame_post{m} = ceil((post{m}.ops.start_frame-1)/3+1);

    if all_start_frame_pre{m} >= all_start_frame_post{m}
        start_frame_later{m} = all_start_frame_pre{m};
    else
        start_frame_later{m} = all_start_frame_post{m};
    end

    movie_sect_length{m} = nframe{m} - start_frame_later{m} + 1;

    start_frame_pre{m} = all_start_frame_pre{m};
    start_frame_post{m} = all_start_frame_post{m};
    end_frame_pre{m} = all_start_frame_pre{m} + movie_sect_length{m} - 1;
    end_frame_post{m} = all_start_frame_post{m} + movie_sect_length{m} - 1;

    movie_frames_pre{m} = start_frame_pre{m}:end_frame_pre{m};
    movie_frames_post{m} = start_frame_post{m}:end_frame_post{m};
end
