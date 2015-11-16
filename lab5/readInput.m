function [tot_mat, tot_vec_lbl, rest_tot_mat, rest_vec_lbl] = readInputMat(baseDir)

% Read the necessary mat files
disp('Load existing Body features');
load([baseDir, 'bdyFeat_fin.mat']);
disp('Load existing topic features');
load([baseDir, 'tpcFeat_fin.mat']);

% Transform the read datasets so it is able to be analyzed
disp('Transform the read datasets so it is able to be analyzed');
[bdyFeatMat, ~, bdyVectLabel] = transformTpcBdyFeat(bdyVectLabel, tpcVectLabel, ...
                                bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js);  

% Keep the rows that contain at least one feature in it.
keep_row_idcs = find(sum(bdyFeatMat,2));

% tot_mat = bdyFeatMat(keep_row_idcs,:);
tot_mat = bdyFeatMat(keep_row_idcs,:);
tot_vec_lbl = bdyVectLabel(:,1);

% Remove heuristically found unnecessary data.
rmv_col = find(strcmp(tot_vec_lbl, 'nil'));

tot_vec_lbl(rmv_col,:) = [];
tot_mat(:,rmv_col) = [];

% Add some excluded features so the feature matrix covers all the documents
[tot_mat, tot_vec_lbl] = genSmallerDataset(tot_mat, tot_vec_lbl, bdyVectLabel);

end