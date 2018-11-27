function A = GetA(alpha, beta, scale, pt1, pt2, F, mode)

    ca = cos(beta);
    sa = sin(beta);
    cb = cos(alpha);
    sb = sin(alpha);
    u1 = pt1(1);
    v1 = pt1(2);
    u2 = pt2(1);
    v2 = pt2(2);

    f1 = F(1,1);
    f2 = F(1,2);
    f3 = F(1,3);
    f4 = F(2,1);
    f5 = F(2,2);
    f6 = F(2,3);
    f7 = F(3,1);
    f8 = F(3,2);
    f9 = F(3,3);

    A1 = (ca * cb * u1 * f1 + ca * cb * v1 * f2 + ca * cb * f3 + sa * cb * u1 * f4 + sa * cb * v1 * f5 + sa * cb * f6);
    B1 = (-sa * sb * u1 * f1 - sa * sb * v1 * f2 - sa * sb * f3 + ca * sb * f4 * u1 + ca * sb * v1 * f5 + ca * sb * f6);
    C1 = (ca * sb * u1 * f1 + ca * sb * v1 * f2 + ca * sb * f3 + sa * sb * f4 * u1 + sa * sb * v1 * f5 + sa * sb * f6);
    D1 = u2 * f1 + v2 * f4 + f7;

    A2 = (- ca * sb * u1 * f1 - ca * sb * v1 * f2 - ca * sb * f3 - sa * sb * u1 * f4 - sa * sb * f6 - sa * sb * v1 * f5);
    B2 = (- sa * cb * u1 * f1 - sa * cb * v1 * f2 - sa * cb * f3 + ca * cb * u1 * f4 + ca * cb * f6 + ca * cb * v1 * f5);
    C2 = (ca * cb * u1 * f1 + ca * cb * v1 * f2 + ca * cb * f3 + sa * cb * u1 * f4 + sa * cb * f6 + sa * cb * v1 * f5);
    D2 = u2 * f2 + v2 * f5 + f8;

    a = (B2 - C2 * B1 / C1);
    b = (D2 - C2 * D1 / C1);
    c = scale * (A2 - C2 * A1 / C1);

    if isinf(a) || isnan(a) || isinf(b) || isnan(b) || isinf(c) || isnan(c) 
        A = 0;
        return;
    end
    
    r = roots([a b c]);

    A = [];
    best_A = 0;
    best_err = 1e10;
    for ri = 1 : length(r)
        qvi = r(ri);
        if abs(a * qvi^2 + b * qvi + c) > 1e-4
            continue;
        end

        qui = scale / qvi;
        wi = -A1 / C1 * qui - B1 / C1 * qvi - D1 / C1 ;

        Ai = [cos(beta), -sin(beta); sin(beta) cos(beta)] * ...
            [qui, wi; 0, qvi] * ...
            [cos(alpha), -sin(alpha); sin(alpha) cos(alpha)];

        A = [A; Ai];
    end
end