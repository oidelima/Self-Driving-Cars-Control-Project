function new_cline = subdivideTrack(centerline)
    new_cline_preformat = (centerline(:,1:end-1) + centerline(:,2:end))/2;
    new_cline = reshape([centerline(:,1:end-1);new_cline_preformat], 2, []);
    new_cline = [new_cline,centerline(:,end)];
end