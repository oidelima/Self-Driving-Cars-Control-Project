function dg = DinequalityConstraintTrack(pos, leftBound, rightBound)
    dg = gradient(inequalityConstraintTrack(pos, leftBound, rightBound), 0.05);
end