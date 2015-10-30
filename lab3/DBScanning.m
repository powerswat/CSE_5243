function [cluster] = DBScanning(tot_mat, min_pts, eps, dist_opt)

tot_row_size = size(tot_mat,1);

dist_mat = zeros(tot_row_size, tot_row_size);

cluster = zeros(tot_row_size,1);
cluster_id = 1;
cluster_pts = zeros(tot_row_size,1);
push_idx = 1;
pop_idx = 1;

for i=1:tot_row_size
    num_min_pts = 0;
    if cluster(i)== 0
        cluster(i) = cluster_id;
    else
        continue;
    end
    
    j = 1;
    k = i;
    while true
        if k==j
            j = j + 1;
            continue;
        end
        if j>=tot_row_size
            rev_pop_idx = min(find(cluster_pts>0));
            if isempty(rev_pop_idx)
                break;
            else
                pop_idx = rev_pop_idx;
            end
            k = cluster_pts(pop_idx);
            cluster_pts(pop_idx) = 0;
            pop_idx = mod(pop_idx + 1, tot_row_size);
            if k == 0
                break;
            end
            num_min_pts = 0;
            j = 1;
            continue;
        end
        if cluster_pts(2)<0
            cluster_pts = zeros(min_pts,1);
        end 
        
        if strcmp(dist_opt, '')
            dist_mat(k,j) = pdist([tot_mat(k,:); tot_mat(j,:)], 'minkowski', 2);
        else
            dist_mat(k,j) = pdist([tot_mat(k,:); tot_mat(j,:)], 'minkowski', 1);
        end
        dist_mat(j,k) = dist_mat(k,j);
        
        if(dist_mat(k,j) <= eps) && cluster(j) == 0
            cluster(j) = cluster_id;
            num_min_pts = num_min_pts + 1;
            cluster_pts(push_idx) = j;
            push_idx = mod(push_idx + 1, tot_row_size);
        end
        
        if num_min_pts >= min_pts
            k = cluster_pts(pop_idx);
            cluster_pts(pop_idx) = 0;
            pop_idx = mod(pop_idx + 1, tot_row_size);
            if k == 0
                break;
            end
            num_min_pts = 0;
            j = 1;
            
            continue;
        end
        
        j = j + 1;
    end
    
    cluster_id = cluster_id + 1;
end

end