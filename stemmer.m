function [strData, stemDict, isTerminate] = stemmer(url1, strData, iterNo, stemDict)

if ~exist('url1', 'var') || isempty(url1), url1 = ...
    sprintf('http://www.thefreedictionary.com/'); end
if ~exist('strData', 'var') || isempty(strData), strData = {'1010'}; end
if ~exist('stemDict', 'var') || isempty(stemDict), stemDict = {'',''}; end

baseDir = 'C:\Temp\CSE_5243\';
dictPart1 = stemDict;
isTerminate = 0;
url2 = 'https://en.wiktionary.org/wiki/';

% Check every word in the input stream
stemDict = cell(length(strData), 2);
for i=1:length(stemDict)
    
    % Check if the query word is already added in the DB.
    if ~isempty(find(strcmp(dictPart1(:,1), strData{i})))
        strData{i} = char(dictPart1(min(find(strcmp(dictPart1(:,1), strData{i}))),2));
        continue;
    end
    
    % Check if the word doesn't need to be queried
    alphaOnly = regexp(strData{i}, '[A-Z,a-z]');
    if length(alphaOnly) < length(strData{i})
        stemDict{i,1} = strData{i};
        stemDict{i,2} = strData{i};
        continue;
    end
    
    % Read the dictionary page source for each given word
    [sourcefile, status] = urlread(sprintf([url1, strData{i}]));  
    if (isempty(sourcefile))
        disp('Access for site 1 is blocked!!')
        iterNo
        isTerminate = 1;
        return;
    end

    stStemIdcs_p1 = strfind(sourcefile, 'tense of ')';
    stStemIdcs_p2 = strfind(sourcefile, 'Past participle of ')';
    if (length(stStemIdcs_p1)+length(stStemIdcs_p2)==0)
        [sourcefile, status] = urlread(sprintf([url2, strData{i}]));  
        if (isempty(sourcefile))
            disp('Access for site 2 is blocked!!')
            iterNo
            isTerminate = 1;
            return;
        end
        stStemIdcs = strfind(sourcefile, '<a href="/wiki/')';
        edStemIdcs = (strfind(sourcefile, '#English" title')-1)';
    else
        stStemIdcs = union(stStemIdcs_p1, stStemIdcs_p2);
        edStemIdcs = strfind(sourcefile, '">')';
    end
    
    num_elem = min(length(stStemIdcs), length(edStemIdcs));
    tgtSrc = cell(num_elem,1);

    % Find minimum distance </span> from each <span .. dbox-bold>, retreive
    % only relevant parts and extract only necessary words from them
    for j=1:num_elem
        if num_elem == length(stStemIdcs) 
            edStemIdcs(j) = min(edStemIdcs(find(stStemIdcs(j) < edStemIdcs)));
        else
            stStemIdcs(j) = max(stStemIdcs(find(stStemIdcs < edStemIdcs(j))));
        end
        tgtSrc{j} = sourcefile(stStemIdcs(j):edStemIdcs(j));
        extIdc = regexp(tgtSrc{j}, '[\w\s\<\>\/]');
        tgtSrc{j} = tgtSrc{j}(extIdc);
        if (length(stStemIdcs_p1)+length(stStemIdcs_p2)==0)
            tgtSrc{j} = strrep(tgtSrc{j}, '<a href/wiki/', '');
        else
            if length(stStemIdcs_p1) > 0
                tgtSrc{j} = strrep(tgtSrc{j}, 'tense of  <a href', '');
            else
                tgtSrc{j} = strrep(tgtSrc{j}, 'Past participle of <a href', '');
            end
        end
        
        tgtSrc{j} = strtrim(tgtSrc{j});
    end
    
    tgtSrc = unique(tgtSrc);
    if (length(tgtSrc)>1)
        correct = 0;
        while correct==0;
            [~, minIdx] = min(cellfun('length', tgtSrc));
            if strcmp(strData{i}(1), tgtSrc{minIdx}(1)) == 1
                if length(strfind(tgtSrc{minIdx},' ')>0)
                    tgtSrc = strData{i};
                else
                    tgtSrc = tgtSrc{minIdx};
                end
                correct = 1;
            else
                tgtSrc(minIdx) = [];
            end
        end
    else
        tgtSrc = char(tgtSrc);
    end
    stemDict{i,1} = strData{i};
    stemDict{i,2} = lower(tgtSrc);
%     pause(10);
end

if (~isempty(dictPart1))
    stemDict = [dictPart1;stemDict];
    emptyCells = cellfun('isempty', stemDict); 
    stemDict(all(emptyCells,2),:) = [];
end

if mod(iterNo, 10)==0
    sd = stemDict;
    [~,idx]=unique(strcat(sd(:,1),sd(:,2)), 'rows');
    stemDict = sd(idx,:);
    save([baseDir, 'stemDictionary.mat'], 'stemDict');
end

a = 1; 