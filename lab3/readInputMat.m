function [tot_mat, tot_vec_lbl, bdyVectLabel] = readInputMat(baseDir)

% Read the necessary mat files
disp('Load existing Body features');
load([baseDir, 'bdyFeat_fin.mat']);
disp('Load existing topic features');
load([baseDir, 'tpcFeat_fin.mat']);

% Transform the read datasets so it is able to be analyzed
disp('Transform the read datasets so it is able to be analyzed');
[bdyFeatMat, tpcFeatMat, bdyVectLabel] = transformTpcBdyFeat(bdyVectLabel, tpcVectLabel, ...
                                bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js);  

% Integrate the two matrix so it is able to be clustered (Only contains the
% documents that has both of topics and bodies.
tpc_idx_map = zeros(size(tpcFeatMat,1),2);
bdy_idx_map = zeros(size(bdyFeatMat,1),2);
keep_tpc_idcs = find(sum(tpcFeatMat,2));
keep_bdy_idcs = find(sum(bdyFeatMat,2));
keep_idcs = intersect(keep_tpc_idcs, keep_bdy_idcs);

% tot_mat = [tpcFeatMat(keep_idcs,:)];
% tot_vec_lbl = [tpcVectLabel];
tot_mat = [tpcFeatMat(keep_idcs,:), bdyFeatMat(keep_idcs,:)];
tot_vec_lbl = [tpcVectLabel; bdyVectLabel(:,1)];
% tot_mat = bdyFeatMat(keep_idcs,:);
% tot_vec_lbl = bdyVectLabel(:,1);

% Remove heuristically found unnecessary data.
rmv_col = find(strcmp(tot_vec_lbl, 'nil'));
% zero_vec_row = find(sum(tot_mat, 2)==0);

tot_vec_lbl(rmv_col,:) = [];
tot_mat(:,rmv_col) = [];
% tot_mat(zero_vec_row,:) = [];

% Add some excluded features so the feature matrix covers all the documents
[tot_mat, tot_vec_lbl, bdyVectLabel] = genHalfDataset(tot_mat, tot_vec_lbl, bdyVectLabel);

end