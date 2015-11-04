function [cluster] = KmeansCluster(tot_mat, num_centroids, dist_mode)

% Count the number of rows and columns of the input matrix 
[num_rows, num_feats] = size(tot_mat);

% Find initial centroids
% centroids = randperm(num_rows, num_centroids);
centroids = randperm(num_rows);
centroids = centroids(1:num_centroids);

% Calculate distance between each point and all the centroids
dist_mat = zeros(num_centroids, num_rows);
cluster = zeros(num_rows, 1);
for i=1:num_centroids
    dist_mat(i,:) = calcDist(tot_mat(centroids(i),:), tot_mat, 2);
    cluster(centroids(i)) = find(centroids == centroids(i));
end

% Assign all the data points into the closest cluster
for i=1:num_rows    
    min_dist = min(dist_mat(:,i));
    min_idx = min(find(dist_mat(:,i) == min_dist));
    cluster(i) = min_idx;
end

num_identical = 0;
prev_num_identical = -1;
converge_cnt = 0;
while true 
    if (num_identical == num_rows) || (converge_cnt > 4)
        break;
    end
        
    prev_cluster = cluster;
    prev_num_identical = num_identical;
    
    % Calculate a mean point for each cluster
    mean_points = zeros(num_centroids, num_feats);
    for i=1:num_centroids
        curr_cluster_idcs = find(cluster == i);
        mean_calc_template = tot_mat(curr_cluster_idcs, :);
        if isempty(mean_calc_template)
            mean_points(i,:) = tot_mat(cluster(i), :);
        else
            mean_points(i,:) = mean(mean_calc_template);
        end
    end
    
    % Calculate distance between each point and all the centroids
    for i=1:num_centroids
        dist_mat(i,:) = calcDist(mean_points(i,:), tot_mat, 2);
    end
    
    % Assign all the data points into the closest cluster
    for i=1:num_rows    
        min_dist = min(dist_mat(:,i));
        min_idx = min(find(dist_mat(:,i) == min_dist));
        cluster(i) = min_idx;
    end
    
    num_identical = length(find(prev_cluster == cluster));
    if abs(prev_num_identical - num_identical) < 10
        converge_cnt = converge_cnt + 1;
    end
end

end