function lab2Driver

% % Driver program for the lab 2. (Please use this to run the program)

ratio_set = [0.7 0.8 0.9];
on_eff = zeros(length(ratio_set),4);
off_eff = zeros(length(ratio_set),4);
accu = zeros(length(ratio_set),4);
    
%% Read all the necessary input data
baseDir = 'C:\Temp\CSE_5243\';
disp('Read all the necessary input data');
if exist([baseDir, 'TFIDF_fin.mat'], 'file')
    disp('Load existing TFIDF');
    load([baseDir, 'TFIDF_fin.mat']);
end

[TFnewIdx, bdy_is, bdy_js, bdy_vals, bdyVectLabel, tpc_is, tpc_js, ...
    tpcVectLabel, cntVec] = readInputTxt(baseDir);

%% Transform the read datasets so it is able to be analyzed
disp('Transform the read datasets so it is able to be analyzed');
[bdyFeatMat, tpcFeatMat, bdyVectLabel] = transformDataSet(TFnewIdx, bdyVectLabel, tpcVectLabel, ...
                                bdy_is, bdy_js, bdy_vals, tpc_is, tpc_js);  
    
%% Run the classifiers with different training/testing ratios
disp('Run the classifiers with different training/testing ratios');
for i=1:length(ratio_set)
    disp(['Start experiment with (train/test): ', num2str(ratio_set(i)), '/', ...
            num2str(1-ratio_set(i))])
    [off_eff(i,:), on_eff(i,:), accu(i,:)] = textClassify(ratio_set(i), TFnewIdx, ...
                bdyVectLabel, tpcVectLabel, bdy_is, bdy_js, bdy_vals, tpc_is, ...
                tpc_js, TFIDF, tpcFeatMat, bdyFeatMat);
end

%% Save all the results
disp('Save all the results')
saveResults(baseDir, on_eff, off_eff, accu, ratio_set);

end