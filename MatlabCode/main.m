addpath 'data'
rng(0)

test                = 'adam';
threshold           = 3.0;
truncated_threshold = threshold * 3 / 2;

%% Load data
img1 = imread(strcat(test, 'A.png'));
img2 = imread(strcat(test, 'B.png'));
M = load(strcat(test, '.pts'));
N = size(M, 1);

all_pts1 = M(:,1:3)';
all_pts2 = M(:,4:6)';

%% Estimate fundamental matrix
[F,inliersIndex,status] = estimateFundamentalMatrix(M(:,1:2), M(:,4:5), 'Method','MSAC',...
    'NumTrials', 5000, 'DistanceThreshold', 0.35);

if status ~= 0
    disp('Error while estimating the fundamental matrix.');
    return;
end

% Estimate the epipole on the second image
e2 = null(F');
e2 = e2 / e2(3);

%% 1-point RANSAC homography fitting 
best_inliers = [];
best_H = [];
best_score = 0;
avg_time = 0; 
evaluation_time = 0;

disp("1-point RANSAC homography estimation started...")
for iter = 1 : N    
    q1 = M(iter, 7);
    q2 = M(iter, 8);
    alpha = M(iter, 10);
    beta = -M(iter, 9);
    scale = q2 / q1;
        
    tic;
    Hs = GetHomographyFromSift(F, e2, alpha, beta, scale, M(iter, 1:2), M(iter, 4:5));
    time_one_estimation = toc;
    evaluation_time = evaluation_time + time_one_estimation;
    avg_time = avg_time + time_one_estimation;
    
    tic;
    for i = 1 : size(Hs, 3)
        Hi = Hs(:,:,i);
        
        % Count the inliers
        pts2_t = Hi * all_pts1;
        pts2_t = rdivide(pts2_t, pts2_t(3,:));
        residuals = vecnorm(all_pts2 - pts2_t);
        
        inliers = find(residuals < truncated_threshold);        
        score = sum(1 - residuals(inliers) / truncated_threshold);
                
        % Update the so-far-the-best model if needed
        if score > best_score
            best_score = score;
            best_inliers = inliers;
            best_H = Hi;
        end
    end
    evaluation_time = evaluation_time + toc;
end
disp("...finished")

avg_time = avg_time / N;

% LSQ fitting
H = NormalizedDLT(M(best_inliers, 1:3), M(best_inliers, 4:6));
pts2_t = H * all_pts1;
pts2_t = rdivide(pts2_t, pts2_t(3,:));
residuals = vecnorm(all_pts2 - pts2_t);
best_inliers = find(residuals < truncated_threshold);    

disp(strcat('Average time of one estimation = ', num2str(avg_time), ' secs'))
disp(strcat('Time of the full procedure= ', num2str(evaluation_time), ' secs'))
disp(strcat('Average re-reprojection error = ', num2str(mean(residuals(best_inliers))), ' px'))


% Draw points
close all;
image = [img1 img2];

imshow(image)
hold on;

for i = 1 : length(best_inliers)
    color = rand(1,3);
    
    plot([M(best_inliers(i), 1), size(img1, 2) + M(best_inliers(i), 4)], [M(best_inliers(i), 2) M(best_inliers(i), 5)], 'Color', color)
    
    scatter(M(best_inliers(i), 1), M(best_inliers(i), 2), 40, color, 'filled');
    scatter(size(img1, 2) + M(best_inliers(i), 4), M(best_inliers(i), 5), 40, color,'filled');
end
hold off;