flag_terminate = false;
num_successful_runs = 0;
total_time = 0;

while ~flag_terminate
    [Y,U,t_total,t_update,Xobs] = forwardIntegrate();
    info = getTrajectoryInfo(Y,U,Xobs,t_update,TestTrack);
    %disp(info.percent_of_track_completed)
    total_time = total_time + info.t_score;
    if info.percent_of_track_completed == 1
        num_successful_runs = num_successful_runs + 1
        avg_time = total_time / num_successful_runs
        %disp(info.t_score)
    else
        flag_terminate = true;
        plotPart2
    end
end
