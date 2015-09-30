function convertToTxt 

% Convert mat file type feature information to txt file type feature
% information.

baseDir = 'C:\Temp\CSE_5243\';

%% Write the body feature vector information
% if exist([baseDir, 'bdyFeat_fin.mat'], 'file')
%     disp('Load existing body features');
%     load([baseDir, 'bdyFeat_fin.mat']);
% end
% 
% file_id = fopen([baseDir, 'bdyVectLabel.txt'], 'w');
% for i=1:length(bdyVectLabel)
%     fprintf(file_id, '%s\t%d\r\n', bdyVectLabel{i,1}, bdyVectLabel{i,2});
% end
% fclose(file_id);
% 
% file_id = fopen([baseDir, 'bdyFeatMat.txt'], 'w');
% for i=1:length(bdy_is)
%     fprintf(file_id, '%d\t%d\t%d\r\n', bdy_is(i), bdy_js(i), bdy_vals(i));
% end
% fclose(file_id);
% 
% %% Write the topic feature vector information
% if exist([baseDir, 'tpcFeat_fin.mat'], 'file')
%     disp('Load existing topic features');
%     load([baseDir, 'tpcFeat_fin.mat']);
% end
% 
% file_id = fopen([baseDir, 'tpcVectLabel.txt'], 'w');
% for i=1:length(tpcVectLabel)
%     fprintf(file_id, '%s\r\n', tpcVectLabel{i});
% end
% fclose(file_id);
% 
% file_id = fopen([baseDir, 'tpcFeatMat.txt'], 'w');
% for i=1:length(tpc_is)
%     fprintf(file_id, '%d\t%d\r\n', tpc_is(i), tpc_js(i));
% end
% fclose(file_id);
% 
% file_id = fopen([baseDir, 'cntVec.txt'], 'w');
% for i=1:size(cntVec,1)
%     for j=1:size(cntVec,2)-1
%         fprintf(file_id, '%d\t', cntVec(i,j));
%     end
%     fprintf(file_id, '%d\r\n', cntVec(i,j+1));
% end
% fclose(file_id);

%% Write the TFIDF information
if exist([baseDir, 'TFIDF_fin.mat'], 'file')
    disp('Load existing TFIDF');
    load([baseDir, 'TFIDF_fin.mat']);
end

file_id = fopen([baseDir, 'TFnewIdx.txt'], 'w');
for i=1:length(TFnewIdx)
    fprintf(file_id, '%d\r\n', TFnewIdx(i));
end
fclose(file_id);

for i=1:length(TFIDF)
    file_id = fopen([baseDir, 'TFIDF\TFIDF_', num2str(i), '.txt'], 'w');
    for j=1:size(TFIDF{i},1)
        fprintf(file_id, '%s\t%d\t%f\t%f\t%f\r\n', char(TFIDF{i}(j,1)), ...
                    cell2mat(TFIDF{i}(j,2)), cell2mat(TFIDF{i}(j,3)), ...
                    cell2mat(TFIDF{i}(j,4)), cell2mat(TFIDF{i}(j,5)));
    end
    fclose(file_id);
end

end
