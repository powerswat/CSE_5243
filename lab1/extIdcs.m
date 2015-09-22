function [idcSets, nums] = extIdcs(rawTexts, stTags, edTags)

idcSets = cell(length(stTags),1);
tmpIdcs = cell(length(rawTexts),2);
nums = zeros(length(idcSets),1);

for i=1:length(idcSets)
    for j=1:length(rawTexts)
        nums(i) = nums(i) + length(strfind(rawTexts{j}, stTags{i}));

        tmpIdcs{j,1} = strfind(rawTexts{j}, stTags{i})';
        if ~strcmp(stTags{i}, 'newid')
            tmpIdcs{j,2} = strfind(rawTexts{j}, edTags{i})';
        end
        tmpIdcs{j,3} = cell(length(tmpIdcs{j,1}),1);
    end
    idcSets{i} = tmpIdcs;
end

for i=1:length(rawTexts)
    for j=1:length(idcSets{1}{i})-1
        newid_st = idcSets{1}{i}(j);
        newid_ed = idcSets{1}{i}(j+1);
        
        idcSets{2}{i,3}{j} = intersect(find(idcSets{2}{i} < newid_ed), ...
                                    find(idcSets{2}{i} > newid_st));
        idcSets{3}{i,3}{j} = intersect(find(idcSets{3}{i} < newid_ed), ...
                                    find(idcSets{3}{i} > newid_st));
        idcSets{4}{i,3}{j} = intersect(find(idcSets{4}{i} < newid_ed), ...
                                    find(idcSets{4}{i} > newid_st));
        idcSets{5}{i,3}{j} = intersect(find(idcSets{5}{i} < newid_ed), ...
                                    find(idcSets{5}{i} > newid_st));
    end
    newid_st = idcSets{1}{i}(j+1);
    idcSets{2}{i,3}{j+1} = find(idcSets{2}{i} > newid_st);
    idcSets{3}{i,3}{j+1} = find(idcSets{3}{i} > newid_st);
    idcSets{4}{i,3}{j+1} = find(idcSets{4}{i} > newid_st);
    idcSets{5}{i,3}{j+1} = find(idcSets{5}{i} > newid_st);
end

a = 1;