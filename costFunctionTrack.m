function j = costFunctionTrack(pos, centerLine, ~)
    Idc = knnsearch(centerLine', pos');
    %j1 = sum((pos - centerLine(:,Idc)) .^ 2, 1);
    j = sum(100*(1-(Idc / size(centerLine, 2))).^2);
end