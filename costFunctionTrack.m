function j = costFunctionTrack(pos, centerLine, ~)
    Idc = knnsearch(centerLine', pos');
    j1 = (sum((pos - centerLine(:,Idc)) .^ 2, 1));
    j = 10*j1 + (10000*(1-(Idc' ./ size(centerLine, 2))));
end