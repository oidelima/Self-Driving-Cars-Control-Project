function info = getTrajectoryInfo(Y,U,Xobs,t_update,TestTrack)
% info = getTrajectoryInfo(Y,U,Xobs)
%
% Given a trajectory, input, and a cell array of obstacles, return
% information about whether or not the trajectory left the track or crashed
% into any obstacles, and if the input limits were exceeded.
%
% NOTE: the trajectory is only for the vehicle's center of mass, so we
% are not checking if the "corners" of the vehicle leave the track or hit
% any obstacles.
%
% INPUTS:
%   Y           an N-by-2 trajectory in the x and y coordinates of the
%               vehicle's states, where the first column is x and the
%                second column is y OR an N-by-6 trajectory in the full
%               state, where the first column is x and the third column is
%               y
%
%   U           an N-by-2 vector of inputs, where the first column is the
%               steering angle and the second column is the rear wheel
%               driving force
%   
%   Xobs        a 1-by-Nobs cell array where each cell contains a 4-by-2
%               obstacle definition, as would be generated by the
%               generateRandomObstacles function (this is an optional
%               argument, so leave it out if your trajectory doesn't avoid
%               any obstacles)
%
%   TestTrack   a TestTrack object for which TestTrack.cline is the
%               centerline (this is an optional argument; the function will
%               by default try to load the provided TestTrack.mat file)
%
%   t_update    a M-by-1 vector of time that records the time consumption
%               when the control input generation function is called
%
% OUTPUTS:
%   info        a struct containing the following catergories 
%               
%               info.Y : The trajectory given as an input argument to the 
%               function.
%
%               info.U : The inputs given as an input argument to the
%               function.
%
%               info.t_finished : Time in seconds when the finish line is 
%               crossed. An empty vector is returned if finish line is not
%               crossed.
%
%               info.t_end : Time in seconds at end of trajectory.
%
%               info.left_track_position : 2-by-1 vector with x and y 
%               coordinates of location where vehicle first leaves the 
%               track. An empty vector is returned if vehicle never leaves 
%               the track.
%
%               info.left_track_time : Time in seconds when vehicle first 
%               leaves the track. An empty vector is returned if vehicle 
%               never leaves the track.
%
%               info.left_percent_of_track_completed : Percentage of the
%               track completed before vehicle first leaves the track. An 
%               empty vector is returned if vehicle never leaves the track.
%
%               info.crash_position : 2-by-1 vector with x and y
%               coordinates of location where vehicle first hits an
%               obstacle. An empty vector is returned if vehicle never hits
%               an obstacle.
%
%               info.crash_time :  Time in seconds when vehicle first hits
%               an obstacle. An empty vector is returned if vehicle never
%               hits an obstacle.
%
%               info.crash_percent_of_track_completed : Percentage of the
%               track completed before vehicle first hits an obstacle. An
%               empty vector is returned if vehicle never hits an obstacle.
%
%               info.input_exceeded : 1-by-2 boolean with the first entry
%               being true if constraints on input 1 are violated and entry
%               2 being true if constraints on input 2 are violated.
%
%               info.percent_of_track_completed : Percentage of track
%               complete prior to leaving the track, hitting an obstacle,
%               or the control inputs end.
%
%               info.t_score : Score of computational time as 
%                      info.t_finished + M * max(num_exceed_time_limit,0),
%               where M(=10) is some penalty value, num_exceed_time_limit
%               is the number of elements in t_update that are longer than
%               0.5 second. An empty vector is returned if finish line is 
%               not crossed.
%
% Written by: Shreyas Kousik and Matthew Porter
% Created: 14 Dec 2017
% Modified: 24 Oct 2018 
% Modified: 22 Nov 2019 (JL)
