function [TFnewIdx, bdy_is, bdy_js, bdy_vals, bdyVectLabel, tpc_is, tpc_js, ...
            tpcVectLabel, cntVec] = readInputTxt(baseDir)

[files,paths] = findFiles(baseDir, '.txt', 1);

% Remove TFIDF files from the list
non_TFIDF_idcs = find(cellfun(@isempty, strfind(files, 'TFIDF\TFIDF_')));
files = files(non_TFIDF_idcs);
paths = paths(non_TFIDF_idcs);

% Read TFnewIdx.txt
TFnew_idx = find(~cellfun(@isempty, strfind(files, 'TFnewIdx.txt')));
file_id = fopen(paths{TFnew_idx}, 'r');
tmp_TFnewIdx = textscan(file_id, '%f', 'Delimiter', '\t');
TFnewIdx = tmp_TFnewIdx{1,1};
fclose(file_id);

% Read tpcFeatMat.txt
tpcFeatMat_idx = find(~cellfun(@isempty, strfind(files, 'tpcFeatMat.txt')));
file_id = fopen(paths{tpcFeatMat_idx}, 'r');
tmp_tpcFeatMat = textscan(file_id, '%f%f', 'Delimiter', '\t');
tpc_is = tmp_tpcFeatMat{1,1};
tpc_js = tmp_tpcFeatMat{1,2};
fclose(file_id);

% Read bdyFeatMat.txt
bdyFeatMat_idx = find(~cellfun(@isempty, strfind(files, 'bdyFeatMat.txt')));
file_id = fopen(paths{bdyFeatMat_idx}, 'r');
tmp_bdyFeatMat = textscan(file_id, '%f%f%f', 'Delimiter', '\t');
bdy_is = tmp_bdyFeatMat{1,1};
bdy_js = tmp_bdyFeatMat{1,2};
bdy_vals = tmp_bdyFeatMat{1,3};
fclose(file_id);

% Read tpcVectLabel.txt
tpcVectLabel_idx = find(~cellfun(@isempty, strfind(files, 'tpcVectLabel.txt')));
file_id = fopen(paths{tpcVectLabel_idx});
tmp_tpcVectLabel = textscan(file_id, '%s');
tpcVectLabel = tmp_tpcVectLabel{1,1};
fclose(file_id);

% Read bdyVectLabel.txt
bdyVectLabel_idx = find(~cellfun(@isempty, strfind(files, 'bdyVectLabel.txt')));
file_id = fopen(paths{bdyVectLabel_idx});
tmp_bdyVectLabel = textscan(file_id, '%s%f', 'Delimiter', '\t');
bdyVectLabel(:,1) = tmp_bdyVectLabel{1,1};
bdyVectLabel(:,2) = num2cell(tmp_bdyVectLabel{1,2});
fclose(file_id);

% Read cntVec.txt
len_cntVec_col = length(tpcVectLabel);
cntVec_format = '';
for i=1:len_cntVec_col
    cntVec_format = [cntVec_format '%f'];
end

cntVec_idx = find(~cellfun(@isempty, strfind(files, 'cntVec.txt')));
file_id = fopen(paths{cntVec_idx});
tmp_cntVec = textscan(file_id, cntVec_format, 'Delimiter', '\t');
cntVec = zeros(length(tmp_cntVec{1}),len_cntVec_col);
for i=1:len_cntVec_col
    cntVec(:,i) = tmp_cntVec{1,i};
end
fclose(file_id);

end