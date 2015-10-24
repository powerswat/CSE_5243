function [bdyFeatMat bdyVectLabel, tpcFeatMat, tpcVectLabel] = readInputMat(baseDir)

% Read the necessary mat files
disp('Load existing Body features');
load([baseDir, 'bdyFeat_fin.mat']);
disp('Load existing topic features');
load([baseDir, 'tpcFeat_fin.mat']);

% Transform the read datasets so it is able to be analyzed
disp('Transform the read datasets so it is able to be analyzed');
[bdyFeatMat, tpcFeatMat, bdyVectLabel] = transformTpcBdyFeat(bdyVectLabel, tpcVectLabel, ...
                                bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js);  

% Integrate the two matrix so it is able to be clustered
totMat = [tpcFeatMat,bdyFeatMat];
                            
a = 1;

end