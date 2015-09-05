function textPreProc

inputDir = 'C:\Users\Young Suk Cho\Dropbox\CSE 5243 Datamining\CSE5243-master\reuters';

% Read all the sgm file names  (Import findFiles)
[files,paths] = findFiles(inputDir, '.sgm', 1);

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
for i=1:length(paths)
    for j=1:length(bdy_idcs{i})
        bodyTxt{(i-1)*j+j} = strsplit(rawTexts{i}(bdy_idcs{i,1}(j)+length('<body>') ...
        :bdy_idcs{i,2}(j)-1));
    
        % Stemming needed here!!!!!!!!!!!!
    
        [uqWds,~,wdIdx]=unique(bodyTxt{(i-1)*j+j});
        nUqWds = length(uqWds);
        counts = hist(wdIdx,1:nUqWds);
        plot(counts)
        text([1:nUqWds], counts, uqWds)
        drawnow;
    end
end
a = 1 