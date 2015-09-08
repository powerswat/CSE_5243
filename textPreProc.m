function textPreProc

inputDir = 'C:\Users\Young Suk Cho\Dropbox\CSE 5243 Datamining\CSE5243-master\reuters';

% Read all the sgm file names  (Import findFiles)
[files,paths] = findFiles(inputDir, '.sgm', 1);

% Check whether there is already a dictionary or not. If there is, then
% read it.
baseDir = 'C:\Temp\CSE_5243\';
if exist([baseDir, 'stemDictionary.mat'], 'file')
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
bodyTxt = cell(numBodies,1);
term_idx = 0;
for i=1:length(paths)
    for j=1:length(bdy_idcs{i})
%         text = rawTexts{i}(bdy_idcs{i,1}(j)+length('<body>'):bdy_idcs{i,2}(j)-1);
%         while ~isempty(text)
%             [bodyTxt{term_idx}, text] = strtok(text);
%             term_idx = term_idx + 1;
%         end
        term_idx = term_idx + 1;

        bodyTxt{term_idx} = strsplit(rawTexts{i}(bdy_idcs{i,1}(j)+length('<body>') ...
        :bdy_idcs{i,2}(j)-1))';
    
%         if term_idx < 8100
%             continue
%         end
%     
%         % Stemming needed here!!!!!!!!!!!!
%         [bodyTxt{term_idx}, stemDict] = stemmer('', bodyTxt{term_idx}, term_idx, dictPart1);
%         bodyTxt{term_idx}(find(cellfun(@isempty, bodyTxt{term_idx}))) = [];
%     
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
a = 1 