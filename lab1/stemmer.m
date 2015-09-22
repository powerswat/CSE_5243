function [strData, stemDict] = stemmer(url, strData, iterNo, stemDict)

if ~exist('url', 'var') || isempty(url), url = ...
    sprintf('http://dictionary.reference.com/browse/'); end
if ~exist('strData', 'var') || isempty(strData), strData = {'eat'}; end
if ~exist('stemDict', 'var') || isempty(stemDict), stemDict = {'',''}; end

baseDir = 'C:\Temp\CSE_5243\';
dictPart1 = stemDict;
isTerminate = 0;

% Check every word in the input stream
stemDict = cell(length(strData), 2);
for i=1:length(stemDict)
    
    % Check if the query word is already added in the DB.
    if ~isempty(find(strcmp(dictPart1(:,1), strData{i}))) 
        strData{i} = char(dictPart1(min(find(strcmp(dictPart1(:,1), strData{i}))),2));
        continue;
    elseif ~isempty(find(strcmp(stemDict(:,1), strData{i})))
        strData{i} = char(dictPart1(min(find(strcmp(stemDict(:,1), strData{i}))),2));
        continue;
    end
    
    % Check if the word doesn't need to be queried
    alphaOnly = regexp(strData{i}, '[A-Z,a-z]');
    if length(alphaOnly)< length(strData{i})
        stemDict{i,1} = strData{i};
        stemDict{i,2} = strData{i};
        continue;
    end
    
    % Read the dictionary page source for each given word
    [sourcefile, ~] = urlread(sprintf([url, strData{i}, '?s=t']));  
    if (isempty(sourcefile))
        stemDict{i,1} = strData{i};
        stemDict{i,2} = strData{i};
        continue;
    end

    stStemIdcs = strfind(sourcefile, '<span class="me" data-syllable=')';
    tmpEdStemIdcs = strfind(sourcefile, '</span>')';
    edStemIdcs = zeros(length(stStemIdcs),1);
    tgtSrc = cell(length(stStemIdcs),1);

    % Find minimum distance </span> from each <span .. dbox-bold>, retreive
    % only relevant parts and extract only necessary words from them
    for j=1:length(stStemIdcs)
        edStemIdcs(j) = min(tmpEdStemIdcs(find(stStemIdcs(j) < tmpEdStemIdcs)));
        tgtSrc{j} = sourcefile(stStemIdcs(j):edStemIdcs(j));
        extIdc = regexp(tgtSrc{j}, '[\w\s\<\>\/]');
        tgtSrc{j} = tgtSrc{j}(extIdc);
        tgtSrc{j} = strrep(tgtSrc{j}, '<', '');
        stSubIdx = strfind(tgtSrc{j}, 'span classme');
        edSubIdx = strfind(tgtSrc{j}, '>');
        tgtSrc{j}(stSubIdx:edSubIdx) = [];
        tgtSrc{j} = strtrim(tgtSrc{j});
    end
    
    tgtSrc = unique(tgtSrc);
    if (length(tgtSrc)>1)
        [~, minIdx] = min(cellfun('length', tgtSrc));
        if length(strfind(tgtSrc{minIdx},' ')>0)
            tgtSrc = strData{i};
        else
            tgtSrc = tgtSrc{minIdx};
        end
    else
        tgtSrc = char(tgtSrc);
    end
    stemDict{i,1} = strData{i};
    stemDict{i,2} = lower(tgtSrc);
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