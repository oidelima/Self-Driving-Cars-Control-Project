clear
load TestTrack.mat

places = 3;

ineqGrid = ones(

for cline_index = 1:size(TestTrack.cline, 2)
    current_center = TestTrack.cline(:,cline_index);
    current_close_enough = round(current_center, places);
end