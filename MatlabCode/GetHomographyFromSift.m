function H = GetHomographyFromSift(F, e2, alpha, beta, scale, pt1, pt2, A)

    As = GetA(alpha, beta, scale, pt1, pt2, F, 1);
    
    if As == 0
        H = 0;
        return;
    end
    
    H = zeros(3,3,0);
    for ri = 1 : size(As, 1) / 2            
        Ai = As((ri - 1) * 2 + 1 : ri * 2, :);       
        H(:,:,ri) = ComputeHAF(reshape(Ai', 1, 4), e2, F, pt1, pt2);
    end     
end
