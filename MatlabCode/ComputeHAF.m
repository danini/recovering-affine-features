%-----------------------------------------------------------------%
% Copyright 2014-2016, Daniel Barath  barath.daniel@sztaki.mta.hu %
%-----------------------------------------------------------------%
% Affines
% 	Type: Matrix
%	Size: N * 4
%	Structure: Each row contains the parameteters of an affine transformation.
%	a11, a12, a21, a22
% e2
%	Type: Vector
%	Size: 3 * 1
%	The epipoles on the second images
% F
%	Type: Matrix
%	Size: 3 * 3
%	The fundamental matrix
% pts1, pts2
%	Type: Vector
%	Size: N * 2
%	The points on the first and second images
function H = ComputeHAF(Affines, e2, F, pts1, pts2)

A	= zeros(size(Affines,1)*6, 3);
b   = zeros(size(Affines,1)*6, 1);

for i = 1 : 6 : size(Affines,1)*6
	currRow		= (i-1) / 6 + 1;
	pt1			= pts1(currRow,:);
	pt2			= pts2(currRow,:);
	
	A(i, 1)		= Affines(currRow,1) * pt1(1) + pt2(1) - e2(1);
	A(i, 2)		= Affines(currRow,1) * pt1(2);
	A(i, 3)		= Affines(currRow,1);
	
	A(i+1, 1)	= Affines(currRow,2) * pt1(1);
	A(i+1, 2)	= Affines(currRow,2) * pt1(2) + pt2(1) - e2(1);
	A(i+1, 3)	= Affines(currRow,2);
	
	A(i+2, 1)	= Affines(currRow,3) * pt1(1) + pt2(2) - e2(2);
	A(i+2, 2)	= Affines(currRow,3) * pt1(2);
	A(i+2, 3)	= Affines(currRow,3);
	
	A(i+3, 1)	= Affines(currRow,4) * pt1(1);
	A(i+3, 2)	= Affines(currRow,4) * pt1(2) + pt2(2) - e2(2);
	A(i+3, 3)	= Affines(currRow,4);
	
	A(i+4, 1)	= (pt1(1) * e2(1) - pt1(1) * pt2(1));
	A(i+4, 2)	= (pt1(2) * e2(1) - pt1(2) * pt2(1));
	A(i+4, 3)	= (e2(1) - pt2(1));
	
	A(i+5, 1)	= (pt1(1) * e2(2) - pt1(1) * pt2(2));
	A(i+5, 2)	= (pt1(2) * e2(2) - pt1(2) * pt2(2));
	A(i+5, 3)	= (e2(2) - pt2(2));
    
    b(i)        = F(2,1);
    b(i+1)      = F(2,2);
    b(i+2)      = -F(1,1);
    b(i+3)      = -F(1,2);
    b(i+4)      = -(pt1(1) * F(2,1) + pt1(2) * F(2,2) + F(2,3));
    b(i+5)      = (pt1(1) * F(1,1) + pt1(2) * F(1,2) + F(1,3));
end;
		
V = pinv(A) * b;

h31		= V(1);
h32		= V(2);
h33		= V(3);

h21		= e2(2) * h31 - 1 * F(1,1);
h22		= e2(2) * h32 - 1 * F(1,2);
h23		= e2(2) * h33 - 1 * F(1,3);
h11		= e2(1) * h31 + 1 * F(2,1);
h12		= e2(1) * h32 + 1 * F(2,2);
h13		= e2(1) * h33 + 1 * F(2,3);

H		= [h11, h12, h13; h21, h22, h23; h31, h32, h33];
end

