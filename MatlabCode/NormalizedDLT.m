%% The normalized DLT algorithm
function H = NormalizedDLT(pts1, pts2)
    [normedPts1, normedPts2, T1, T2] = NormalizePoints(pts1, pts2);
    H = DLT(normedPts1, normedPts2);
    H = inv(T2) * H * T1;
end