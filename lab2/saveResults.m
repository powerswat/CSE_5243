function saveResults(baseDir, on_eff, off_eff, accu, ratio_set)

% Save all the results

file_id = fopen([baseDir, 'online_efficiency.txt'], 'w');
fprintf(file_id, '%s\t%s\t%s\t%s\t\t%s\r\n', 'Train/Test', 'KNN_full', ...
            'KNN_half', 'NB_full', 'NB_half');
for i=1:size(on_eff,1)
    fprintf(file_id, '%s\t\t', [num2str(ratio_set(i)) '/' num2str((1-ratio_set(i)))]);
    fprintf(file_id, '%f\t%f\t%f\t%f\r\n', on_eff(i,:));
end
fclose(file_id);

file_id = fopen([baseDir, 'offline_efficiency.txt'], 'w');
fprintf(file_id, '%s\t%s\t%s\t%s\t\t%s\r\n', 'Train/Test', 'KNN_full', ...
            'KNN_half', 'NB_full', 'NB_half');
for i=1:size(off_eff,1)
    fprintf(file_id, '%s\t\t', [num2str(ratio_set(i)) '/' num2str((1-ratio_set(i)))]);
    fprintf(file_id, '%f\t%f\t%f\t%f\r\n', off_eff(i,:));
end
fclose(file_id);

file_id = fopen([baseDir, 'accuracy.txt'], 'w');
fprintf(file_id, '%s\t%s\t%s\t%s\t\t%s\r\n', 'Train/Test', 'KNN_full', ...
            'KNN_half', 'NB_full', 'NB_half');
for i=1:size(accu,1)
    fprintf(file_id, '%s\t\t', [num2str(ratio_set(i)) '/' num2str((1-ratio_set(i)))]);
    fprintf(file_id, '%f\t%f\t%f\t%f\r\n', accu(i,:));
end
fclose(file_id);

a= 1;

end