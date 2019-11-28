function g = inequalityConstraintTrackPrecomp(leftBound, rightBound)
    xlim1 = 200;
    xlim2 = 1500;
    ylim1 = -200;
    ylim2 = 825;
    
    [rx, ry] = meshgrid(xlim1:.1:xlim2,ylim1:.1:ylim2);
    rx2 = reshape(rx,[1,size(rx,1)*size(rx,2)]);
    ry2 = reshape(ry,[1,size(ry,1)*size(ry,2)]);
    g = inequalityConstraintTrack([rx2;ry2], leftBound, rightBound);
end