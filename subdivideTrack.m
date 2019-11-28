function new_cline = subdivideTrack(centerline)
    num_rows = size(centerline, 1);
    new_cline_preformat = (centerline(:,1:end-1) + centerline(:,2:end)) ./ 2;
    new_cline = reshape([centerline(:,1:end-1);new_cline_preformat], num_rows, []);
    new_cline = [new_cline,centerline(:,end)];
end