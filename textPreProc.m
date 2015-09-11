function textPreProc

inputDir = 'C:\Users\Young Suk Cho\Dropbox\CSE 5243 Datamining\CSE5243-master\reuters\';

% Read all the sgm file names  (Import findFiles)
[files,paths] = findFiles(inputDir, '.sgm', 1);
[files, s_idx] = sort(files);
paths = paths(s_idx);

% Check whether there is already a dictionary or not. If there is, then
% read it.
baseDir = 'C:\Temp\CSE_5243\';
isDictComplete = 0;
if exist([baseDir, 'stemDict_fin.mat'], 'file')
    disp('Load Dictionary');
    load([baseDir, 'stemDict_fin.mat']);
    isDictComplete = 1;
end
isBodyComplete = 0;
if exist([baseDir, 'bodyTxt_fin.mat'], 'file')
    disp('Load Body text');
    load([baseDir, 'bodyTxt_fin.mat']);
    isBodyComplete = 1;
end

if (exist([baseDir, 'stemDictionary.mat'], 'file') && isDictComplete == 0)
    load([baseDir, 'stemDictionary.mat']);
    dictPart1 = stemDict;
else
    dictPart1 = '';
end

% Data sets
rawTexts = cell(length(paths),1);
for i=1:length(paths)
    fileID = fopen(paths{i},'r');
    rawTexts{i} = lower(fileread(paths{i}));
    extIdc = regexp(rawTexts{i}, '[\w\s\<\>\/]');
    rawTexts{i} = rawTexts{i}(extIdc);
    disp(['Read txt file: ', num2str(i) , '/', num2str(length(paths))]);
end

% Locate every <TOPIC>, <PLACES>, and <BODY> index
numBodies = 0;
bdy_idcs = cell(length(paths),2);
tpc_idcs = cell(length(paths),2);
plc_idcs = cell(length(paths),2);

for i=1:length(paths)
    numBodies = numBodies + length(strfind(rawTexts{i}, '<body>'));
    
    bdy_idcs{i,1} = strfind(rawTexts{i}, '<body>');
    tpc_idcs{i,1} = strfind(rawTexts{i}, '<topics>');
    plc_idcs{i,1} = strfind(rawTexts{i}, '<places>');
    
    bdy_idcs{i,2} = strfind(rawTexts{i}, '</body>');
    tpc_idcs{i,2} = strfind(rawTexts{i}, '</topics>');
    plc_idcs{i,2} = strfind(rawTexts{i}, '</places>');
    
end

% Collect <TOPIC>, <PLACES>, and <BODY> information together and link them
% into a single data entry
if isBodyComplete == 0
    bodyTxt = cell(numBodies,1);
    term_idx = 0;
    for i=1:length(paths)
        for j=1:length(bdy_idcs{i})

            term_idx = term_idx + 1;

            bodyTxt{term_idx} = strread(rawTexts{i}(bdy_idcs{i,1}(j)+length('<body>') ...
            :bdy_idcs{i,2}(j)-1),'%s','delimiter',' ');
            txtOnlyIdcs = find(~cellfun(@isempty, regexp(bodyTxt{i}, '[A-z,a-z]')));
            bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
            txtOnlyIdcs = find(cellfun(@isempty, strfind(bodyTxt{i}, '<')));
            bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
            txtOnlyIdcs = find(cellfun(@isempty, strfind(bodyTxt{i}, '>')));
            bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
            txtOnlyIdcs = find(cellfun(@isempty, strfind(bodyTxt{i}, '/')));
            bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);

            % Stem the word if necessary
            if isDictComplete == 0
                [bodyTxt{term_idx}, stemDict] = stemmer('', bodyTxt{term_idx}, term_idx, dictPart1);
                bodyTxt{term_idx}(find(cellfun(@isempty, bodyTxt{term_idx}))) = [];
                dictPart1 = stemDict;            
            end

    %         [uqWds,~,wdIdx]=unique(bodyTxt{(i-1)*j+j});
    %         nUqWds = length(uqWds);
    %         counts = hist(wdIdx,1:nUqWds);
    %         plot(counts)
    %         text([1:nUqWds], counts, uqWds)
    %         drawnow;

            if (mod(term_idx, 100)==0)
                disp(['Current Iter: ', num2str(term_idx), '/', num2str(numBodies)]);
            end
        end
    end
end

% Optimize the stemdict by removing the rows of which an original word and
% its corresponding queried word stays the same
if isDictComplete == 0
    fin_row = find(~strcmp(stemDict(:,1), stemDict(:,2)));
    stemDict = stemDict(fin_row, 1:2);
    fin_row = find(~cellfun(@isempty, stemDict(:,2)));
    stemDict = stemDict(fin_row, 1:2);
    isLetter = ~cellfun(@isempty, regexp(stemDict, '[A-Z,a-z]'));
    fin_row = find(and(isLetter(:,1), isLetter(:,2)));
    stemDict = stemDict(fin_row, 1:2);

    save([baseDir, 'stemDict_fin.mat'], 'stemDict');
    
elseif isBodyComplete == 0
    % If the terms are not stemmed, stem them 
    for i=1:numBodies
        for j=1:length(bodyTxt{i})
            if ~isempty(find(strcmp(stemDict(:,1), bodyTxt{i}(j))))
                bodyTxt{i}(j) = stemDict(find(strcmp(stemDict(:,1), bodyTxt{i}(j))), 2);
            end
        end
        
        if (mod(i, 100)==0)
            disp(['Current Iter: ', num2str(i), '/', num2str(numBodies)]);
        end
    end
    
    save([baseDir, 'bodyTxt_fin.mat'], 'bodyTxt');
    
else
    % Count the number of each term appears in the set of documents
    TFIDF = cell(numBodies,1);
    for i=1:numBodies
%         txtOnlyIdcs = find(~cellfun(@isempty, regexp(bodyTxt{i}, '[A-z,a-z]')));
%         bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
%         txtOnlyIdcs = find(cellfun(@isempty, strfind(bodyTxt{i}, '<')));
%         bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
%         txtOnlyIdcs = find(cellfun(@isempty, strfind(bodyTxt{i}, '>')));
%         bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
%         txtOnlyIdcs = find(cellfun(@isempty, strfind(bodyTxt{i}, '/')));
%         bodyTxt{i} = bodyTxt{i}(txtOnlyIdcs, :);
        
        TFIDF(i) = {tabulate(bodyTxt{i})};
        TFIDF(i) = {[TFIDF{i}, cell(size(TFIDF{i},1), 1)]};
        keep_idx = find(cell2mat(TFIDF{i}(:,2))>1);
        TFIDF{i} = TFIDF{i}(keep_idx,:);
        tfMax = max(cell2mat(TFIDF{i}(:,2)));
        TFIDF{i}(:,4) = num2cell(0.5 + ((0.5 * cell2mat(TFIDF{i}(:,2))) / tfMax));
        if (mod(i, 1000)==0)
            disp(['Generate TF mat: ', num2str(i), '/', num2str(numBodies)]);
        end
    end
    
    numTerms = 0;
    for i=1:numBodies
        numTerms = numTerms + size(TFIDF{i},1);
    end
    
    % Generate IDF matrix
    termMap = cell(numTerms,1);
    idx = 1;
    for i=1:numBodies
        termMap(idx:idx+size(TFIDF{i},1)-1) = TFIDF{i}(:,1);
        idx = idx + size(TFIDF{i},1);
        if (mod(i, 1000)==0)
            disp(['Generate termSet: ', num2str(i), '/', num2str(numBodies)]);
        end
    end
    
    IDF = tabulate(termMap);
    txtOnlyIdcs = find(~cellfun(@isempty, regexp(IDF(:,1), '[A-z,a-z]')));
    IDF = IDF(txtOnlyIdcs, :);
    IDF = [IDF, cell(size(IDF,1), 2)];
    IDF(:,4) = num2cell(log10((length(IDF)./cell2mat(IDF(:,2)))));
    keep_idx = find(cell2mat(IDF(:,2))>1);
    IDF = IDF(keep_idx,:);
    [~,s_idx] = sort(cell2mat(IDF(:,2)));
    IDF = IDF(s_idx,:);
    
    % Calculate TF-IDF
    for i=1:numBodies
        [~, iTFIDF, iIDF] = intersect(TFIDF{i}(:,1), IDF(:,1));
        TFIDF{i}(iTFIDF,5) = num2cell(cell2mat(TFIDF{i}(iTFIDF,4)) .* cell2mat(IDF(iIDF,4)));
        
        keep_idx = find(~cellfun(@isempty, TFIDF{i}(:,5)));
        TFIDF{i} = TFIDF{i}(keep_idx,:);
        [~,s_idx] = sort(cell2mat(TFIDF{i}(:,5)), 'descend');
        TFIDF{i} = TFIDF{i}(s_idx,:);
        if (mod(i, 1000)==0)
            disp(['Generate TFIDF mat: ', num2str(i), '/', num2str(numBodies)]);
        end
    end
end

a = 1 