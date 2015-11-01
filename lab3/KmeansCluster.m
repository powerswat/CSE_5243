function [cluster] = KmeansCluster(tot_mat, num_centroids, dist_mode)

% Count the number of rows and columns of the input matrix 
[num_rows, num_feats] = size(tot_mat);

% Find initial centroids
centroids = randperm(num_rows, num_centroids);

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

prev_cluster = 0;
while length(find(prev_cluster == cluster)) < num_rows
    prev_cluster = cluster;
    
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
        if length(find(isnan(mean_points(i,:)))) > 0
            a = 1;
        end
    end
    
    % Calculate distance between each point and all the centroids
    for i=1:num_centroids
        dist_mat(i,:) = calcDist(mean_points(i,:), tot_mat, 2);
        if sum(dist_mat(i,:)) == 0
            a = 1;
        end
    end
    
    % Assign all the data points into the closest cluster
    for i=1:num_rows    
        min_dist = min(dist_mat(:,i));
        min_idx = min(find(dist_mat(:,i) == min_dist));
        cluster(i) = min_idx;
    end
end

end