function textPreProc

inputDir = 'C:\Users\Young Suk Cho\Dropbox\CSE 5243 Datamining\CSE5243-master\reuters\';

% Read all the sgm file names  (Import findFiles)
[files,paths] = findFiles(inputDir, '.sgm', 1);
[files,s_idx] = sort(files);
paths = paths(s_idx);

% Check whether there is already a dictionary or not. If there is, then
% read it.
baseDir = 'C:\Temp\CSE_5243\';
isDictComplete = 0;
if exist([baseDir, 'stemDict_fin.mat'], 'file')
    disp('Load existing Dictionary');
    load([baseDir, 'stemDict_fin.mat']);
    isDictComplete = 1;
end
isBodyComplete = 0;
if exist([baseDir, 'bodyTxt_fin.mat'], 'file')
    disp('Load existing Body text');
    load([baseDir, 'bodyTxt_fin.mat']);
    isBodyComplete = 1;
end
isBodyFeat = 0;
if exist([baseDir, 'bdyFeat_fin.mat'], 'file')
    disp('Load existing body features');
    load([baseDir, 'bdyFeat_fin.mat']);
    isBodyFeat = 1;
end
isTpcFeat = 0;
if exist([baseDir, 'tpcFeat_fin.mat'], 'file')
    disp('Load existing topic features');
    load([baseDir, 'tpcFeat_fin.mat']);
    isTpcFeat = 1;
end
isTFIDF = 0;
if exist([baseDir, 'TFIDF_fin.mat'], 'file')
    disp('Load existing TFIDF');
    load([baseDir, 'TFIDF_fin.mat']);
    isTFIDF = 1;
end
isPlc = 0;
if exist([baseDir, 'plcFeat_fin.mat'], 'file')
    disp('Load existing place feature');
    load([baseDir, 'plcFeat_fin.mat']);
    isPlc = 1;
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
bdy_idcs = cell(length(paths),3);
tpc_idcs = cell(length(paths),3);
plc_idcs = cell(length(paths),3);
ttl_idcs = cell(length(paths),3);

% Locate all the indices for each tag
stTags = {'newid'; '<title>'; '<topics>'; '<body>'; '<places>'};
edTags = {''; '</title>'; '</topics>'; '</body>'; '</places>'};
disp('Extracting indices ...')
[idcSets, nums] = extIdcs(rawTexts, stTags, edTags);

numBodies = nums(find(~cellfun(@isempty, strfind(stTags, '<body>'))));
numTopics = nums(find(~cellfun(@isempty, strfind(stTags, '<topics>'))));
numTitles = nums(find(~cellfun(@isempty, strfind(stTags, '<title>'))));

ttl_idcs = idcSets{find(~cellfun(@isempty, strfind(stTags, '<title>')))};
tpc_idcs = idcSets{find(~cellfun(@isempty, strfind(stTags, '<topics>')))};
bdy_idcs = idcSets{find(~cellfun(@isempty, strfind(stTags, '<body>')))};
plc_idcs = idcSets{find(~cellfun(@isempty, strfind(stTags, '<places>')))};
newid_idcs = idcSets{find(~cellfun(@isempty, strfind(stTags, 'newid')))};

% Stem for body text
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
            if isDictComplete == 1
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

% Extract topic text
isTpcComplete = 0;
if isTpcComplete == 0
    tpcTxt = cell(numTopics,1);    
    term_idx = 0;
    
    for i=1:length(paths)
        for j=1:length(tpc_idcs{i})

            term_idx = term_idx + 1;

            tpcTxt{term_idx} = char(strread(rawTexts{i}(tpc_idcs{i,1}(j)+length('<title>') ...
            :tpc_idcs{i,2}(j)-1),'%s','delimiter',' '));
            tpcTxt{term_idx} = strrep(tpcTxt{term_idx}, '><d>', '');
            tpcTxt{term_idx} = strrep(tpcTxt{term_idx}, '</d>', '');
            tpcTxt{term_idx} = strrep(tpcTxt{term_idx}, '</d', ' ');
            tpcTxt{term_idx} = strtrim(strrep(tpcTxt{term_idx}, '>', ''));
 
            if (mod(term_idx, 100)==0)
                disp(['Current Iter: ', num2str(term_idx), '/', num2str(numBodies)]);
            end
        end
    end
end

% Extract place text
isPlcComplete = 0;
if isPlcComplete == 0
    plcTxt = cell(numTopics,1);    
    term_idx = 0;
    
    for i=1:length(paths)
        for j=1:length(plc_idcs{i})

            term_idx = term_idx + 1;

            plcTxt{term_idx} = char(strread(rawTexts{i}(plc_idcs{i,1}(j)+length('<places>') ...
            :plc_idcs{i,2}(j)-1),'%s','delimiter',' '));
            plcTxt{term_idx} = strrep(plcTxt{term_idx}, '<d>', ' ');
            plcTxt{term_idx} = strrep(plcTxt{term_idx}, '</d>', '');
            plcTxt{term_idx} = strrep(plcTxt{term_idx}, '</d', ' ');
            plcTxt{term_idx} = strtrim(strrep(plcTxt{term_idx}, '>', ''));
 
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
    if isTFIDF == 1
        % Count the number of each term appears in the set of documents
        TFIDF = cell(numBodies,2);
        for i=1:numBodies        
            TFIDF(i) = {tabulate(bodyTxt{i})};
            TFIDF(i) = {[TFIDF{i}, cell(size(TFIDF{i},1), 1)]};
            keep_idx = find(cell2mat(TFIDF{i}(:,2)));
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
        keep_idx = find(cell2mat(IDF(:,2))>0);
        IDF = IDF(keep_idx,:);
        [~,s_idx] = sort(cell2mat(IDF(:,2)));
        IDF = IDF(s_idx,:);

        % Calculate TF-IDF
        for i=1:numBodies
            [~, iTFIDF, iIDF] = intersect(TFIDF{i}(:,1), IDF(:,1));
            TFIDF{i}(iTFIDF,5) = num2cell(cell2mat(TFIDF{i}(iTFIDF,4)) .* cell2mat(IDF(iIDF,4)));

            keep_idx = find(~cellfun(@isempty, TFIDF{i}(:,5)));
            TFIDF{i} = TFIDF{i}(keep_idx,:);
            keep_idx = find(cellfun(@length, TFIDF{i,1}(:,1))>1);
            TFIDF{i} = TFIDF{i}(keep_idx,:);
            [~,s_idx] = sort(cell2mat(TFIDF{i}(:,5)), 'descend');
            TFIDF{i} = TFIDF{i}(s_idx,:);
          
            if (mod(i, 1000)==0)
                disp(['Generate TFIDF mat: ', num2str(i), '/', num2str(numBodies)]);
            end
        end
        
        maxPrevIdx = 0;
        idxCnt = 1;
        TFnewIdx = zeros(length(plcTxt),1);
        for i=1:length(bdy_idcs)
            for j=1:length(bdy_idcs{i,3})
                if ~isempty(cell2mat(bdy_idcs{i,3}(j)))
                    TFnewIdx(idxCnt) = maxPrevIdx + cell2mat(bdy_idcs{i,3}(j));
                end
                idxCnt = idxCnt + 1;
            end
            maxPrevIdx = maxPrevIdx + length(bdy_idcs{i});
        end
        TFnewIdx = TFnewIdx';
        
        save([baseDir, 'TFIDF_fin.mat'], 'TFIDF', 'TFnewIdx');
    end
end

% Generate feature vector for topics
        tpcVec = tabulate(tpcTxt);
        allTpcs = strjoin(tpcVec(:,1));
        tpcVec = strread(allTpcs,'%s','delimiter',' ');    
        tpcVec = unique(tpcVec);

% Composite the feature vector
if isTpcFeat == 1
    maxWordLen = max(cellfun(@length, strfind(tpcTxt, ' ')));
    tpcMat = cell(length(tpcTxt), maxWordLen);
    for i=1:length(tpcTxt)
        tpcWVec = strread(tpcTxt{i}, '%s', 'delimiter', ' ')';
        tpcMat(i,1:length(tpcWVec)) = tpcWVec;
    end

    cntVec = zeros(length(plcTxt),length(tpcVec));
    rowfVec = length(plcTxt);
    for i=1:rowfVec
        for j=1:length(tpcVec)
            chkVec = find(strcmp(tpcMat(i,:), tpcVec{j}));
            if ~isempty(chkVec)
                cntVec(i,j) = 1;
            end
        end
        if (mod(i, 1000)==0)
            disp(['Generate topic feature vector: ', num2str(i), '/', num2str(rowfVec)]);
        end
    end
    
    tpcFeatMat = cntVec;
    tpcVectLabel = tpcVec;
    [tpc_is,tpc_js,tpc_vals] = find(tpcFeatMat);
    save([baseDir, 'tpcFeat_fin.mat'], 'tpc_is' ,'tpc_js', 'tpcVectLabel', 'cntVec');
end

if isBodyFeat == 1
    [row, ~] = find(cntVec);
    row = unique(row);
    idcs = [1:length(cntVec)]';
%     emptRow = setdiff(idcs,row);
    emptRow = [1:length(cntVec)]';
    bdy_loc = zeros(length(emptRow),2);
    bdy_loc(:,1) = floor(emptRow./1000)+1;
    bdy_loc(:,2) = mod(emptRow, 1000);
    mod0Idx = find(bdy_loc(:,2)==0);
    bdy_loc(mod0Idx,2) = 1000;
    bdy_loc(mod0Idx,1) = bdy_loc(mod0Idx,1) - 1;

    max_len = 0;
    for i=1:length(TFIDF)
        max_len = max(max_len,length(TFIDF{i}));
    end
    
    featCand = cell(length(emptRow)*max_len,2);
    cnt = 1;
    for i=1:length(emptRow)
        if TFnewIdx(emptRow(i))>0
            if ~isempty(TFIDF{TFnewIdx(emptRow(i))})
                for j=1:size(TFIDF{TFnewIdx(emptRow(i))},1)
                    featCand{cnt,1} = char(TFIDF{TFnewIdx(emptRow(i))}(j,1));
                    featCand{cnt,2} = cell2mat(TFIDF{TFnewIdx(emptRow(i))}(j,5));
                    cnt = cnt + 1;
                end
            end
        end
    end
    [keep_idx,~] = find(~cellfun(@isempty, featCand));
    featCand = featCand(keep_idx,:);
    [~,s_idx] = sort(cell2mat(featCand(:,2)), 'descend');
    featCand = featCand(s_idx,:);
    [~,row] = unique(featCand(:,1), 'stable');
    featCand = featCand(row,:);
    
    bdyFeatMat = zeros(length(plcTxt),length(featCand));
    for i=1:length(featCand)
        for j=1:length(TFIDF)
           matchInJ = find(strcmp(TFIDF{j}(:,1), featCand{i}));
           if ~isempty(matchInJ)
               bdyFeatMat(find(TFnewIdx == j),i) = cell2mat(TFIDF{j}(matchInJ,2));
           end
        end
        if (mod(i, 100)==0)
            disp(['Generate body feature vector: ', num2str(i), '/', num2str(length(featCand))]);
        end
    end

    colSums = sum(bdyFeatMat)';
    [~, s_idx] = sort(colSums, 'descend');
    featCand = featCand(s_idx,:);
    bdyFeatMat = bdyFeatMat(:,s_idx);
%     featCand(1:19) = [];
%     bdyFeatMat(:,1:19) = [];
%     featCand(801:end) = [];
%     bdyFeatMat(:,801:end) = [];
   
    bdyVectLabel = featCand;
    [bdy_is,bdy_js,bdy_vals] = find(bdyFeatMat);
    save([baseDir, 'bdyFeat_fin.mat'], 'bdy_is', 'bdy_js', 'bdy_vals', 'bdyVectLabel', '-v7.3');
end

if isPlc == 1
    plcFeatMat = plcTxt;
    save([baseDir, 'plcFeat_fin.mat'], 'plcFeatMat');
end

% fVec = cell(length(plcTxt),length(tpcVec)+length(featCand)+2);
% docIDs = [1:length(plcTxt)]';
% labels = ['doc_id', tpcVec', featCand'];
% 
% fVec(:,1) = num2cell(docIDs);
% fVec(:,end) = plcTxt;
% fVec(:,2:end-1) = num2cell(zeros(length(plcTxt),length(tpcVec)));
% fVec(:,2:length(tpcVec)+1) = num2cell(cntVec);
% fVec(:,length(tpcVec)+2:end-1) = num2cell(bdyFeatMat);
% 
% fileID = fopen([baseDir, 'fVect_fin.txt'],'w');
% for i=1:length(fVec)
%     for j=1:size(fVec,2)
%         fprintf(fileID, '%s\t', num2str(cell2mat(fVec(i,j))));
%     end
%     fprintf(fileID, '\n');
% end
% fclose(fileID);

disp('Done')