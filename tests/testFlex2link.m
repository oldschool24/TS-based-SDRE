function testFlex2link(testConfigPath)
    testConfig = jsondecode(fileread(testConfigPath));
    expPath = fullfile('../runs', testConfig.expName);
    trainConfigPath = fullfile(expPath, 'trainData.json');
    trainConfig = jsondecode(fileread(trainConfigPath));

    [~, testConfigName, ~] = fileparts(testConfigPath);
    folderName = createFolder(expPath, testConfigName);
    folderPath = fullfile(expPath, folderName);
    copyfile(testConfigPath, fullfile(folderPath, 'test.json'))

    q1Range = testConfig.q1Range;
    q2Range = testConfig.q2Range;
    z1Range = testConfig.z1Range;
    z2Range = testConfig.z2Range;
    q = testConfig.q;  % q, r - parameters of the control criterion
    r = testConfig.r;
    known = testConfig.known;
    modelPath = fullfile(expPath, 'model.mat');
    dt = trainConfig.dt;
    isWrap = trainConfig.isWrap;
    T = testConfig.T;
    nExamples = testConfig.nExamples;
    isAnalyze = testConfig.isAnalyze;
    stopType = testConfig.stopType;
    xRange = testConfig.xRange;


    addpath('../')
    if isempty(z1Range) || isempty(z2Range)
        [q1, q2] = ndgrid(q1Range, q2Range);
        z1 = zeros(size(q1));
        z2 = zeros(size(q1));
        for k=1:numel(q1)
            G = rhsG([q1(k) q2(k)]);
            z1(k) = -G(1);
            z2(k) = -G(2);
        end
    else
        [q1, q2, z1, z2] = ndgrid(q1Range, q2Range, z1Range, z2Range);
    end
    nTests = numel(q1);
    idxExamples = randsample(nTests, nExamples);
    Q = q * eye(8);
    R = r * eye(2);
    if isAnalyze
        analyzeDir = fullfile(expPath, "testAnalysis");
        if ~exist(analyzeDir, 'dir')
            mkdir(analyzeDir)
        end
    end
    xRange(1, 1:2) = -inf;
    xRange(2, 1:2) = inf;
    tsX = [];
    predX = [];
    f_true = [];
    f_pred = []; 
    B_true = []; 
    B_pred = [];

    testStats = zeros(nTests, 8+7+8);
    warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
    warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    
    WaitMessage = utils.parfor_wait(nTests, 'Waitbar', true); %append('Testing. Export path is ', expPath), true);

    tic
    parfor k=1:nTests
%     idxExamples = 1;
%     for k=idxExamples:idxExamples
        warning('off', 'all')    
        x0 = [q1(k); q2(k); z1(k); z2(k); 0; 0; 0; 0]; 
        if ismember(k, idxExamples)
            imgDir = fullfile( ...
                folderPath, ['example_' num2str(k)]);
            if ~exist(imgDir, "dir")
                mkdir(imgDir)
            end
        else
            imgDir = '';
        end
        simStats = mainSim(modelPath, 'flex2link', dt, T, x0, Q, R, ...
                           stopType, xRange, @ode15s, isWrap, ...
                           imgDir, known, isAnalyze);
        testStats(k, :) = [
            x0', simStats.tsCriterion, simStats.sdreCriterion, ...
            simStats.tsTime, simStats.sdreTime, ...
            simStats.tsWallTime, simStats.sdreWallTime, ...
            simStats.insideEpsTube, simStats.stopPoint];
        if isAnalyze
            tsX = [tsX; simStats.tsX]; 
            predX = [predX; simStats.predX]; 
            f_true = [f_true; simStats.f_true]; 
            f_pred = [f_pred; simStats.f_pred]; 
            B_true = [B_true; simStats.B_true]; 
            B_pred = [B_pred; simStats.B_pred]; 
        end

        WaitMessage.Send;
    end

    WaitMessage.Destroy;

    testStats = table2map(testStats);
    save(fullfile(folderPath, testConfigName), 'testStats')
    printResults(testStats)

    if isAnalyze
        utils.wrapperStats('flex2link', tsX, predX, f_true, f_pred, ...
                           B_true, B_pred, analyzeDir)
    end
end
    
function res = createFolder(expPath, folderName)
    if ~exist(fullfile(expPath, folderName), 'dir')
        % Folder doesn't exist, create it
        mkdir(fullfile(expPath, folderName));
        res = folderName;
    else
        % Folder already exists, find the highest numbered folder
        existingFolders = dir([fullfile(expPath, folderName) '*']);
        highestNumber = 1;
        for k = 1:length(existingFolders)
            folder = existingFolders(k).name;
            % Extract the number from the folder name
            [~, name, ext] = fileparts(folder);
            if strcmp(ext, '') && startsWith(name, folderName + "_")
                numberStr = extractAfter(name, length(folderName) + 1);
                number = str2double(numberStr);
                if ~isnan(number) && number > highestNumber
                    highestNumber = number;
                end
            end
        end
        % Create the new folder with the incremented number
        newFolderName = [folderName '_' num2str(highestNumber + 1)];
        mkdir(fullfile(expPath, newFolderName));
        res = newFolderName;
    end
end

function G = rhsG(q)
    G = [79.38*sin(q(1)) + 11.074*sin(q(1)+q(2)); 
         11.074*sin(q(1)+q(2))];
end

function map = updateMap(map, key, idx, value)
    arr = map(key);
    arr(idx, :) = value;
    map(key) = arr;
end

function map=table2map(table)
    [nTests, ~] = size(table);
    map = containers.Map( ...
        {'x0', 'tsCriterion', 'sdreCriterion', 'tsTime', 'sdreTime', ...
        'tsWallTime', 'sdreWallTime', 'insideEpsTube', 'stopPoint'}, ...
        {zeros(nTests, 8), zeros(nTests, 1), zeros(nTests, 1), ...
        zeros(nTests, 1), zeros(nTests, 1), zeros(nTests, 1), ...
        zeros(nTests, 1), zeros(nTests, 1), zeros(nTests, 8)});
    
    for iTest=1:nTests
        map = updateMap(map, 'x0', iTest, table(iTest, 1:8));
        map = updateMap(map, 'tsCriterion', iTest, table(iTest, 9));
        map = updateMap(map, 'sdreCriterion', iTest, table(iTest, 10));
        map = updateMap(map, 'tsTime', iTest, table(iTest, 11));
        map = updateMap(map, 'sdreTime', iTest, table(iTest, 12));
        map = updateMap(map, 'tsWallTime', iTest, table(iTest, 13));
        map = updateMap(map, 'sdreWallTime', iTest, table(iTest, 14));
        map = updateMap(map, 'insideEpsTube', iTest, table(iTest, 15));
        map = updateMap(map, 'stopPoint', iTest, table(iTest, 16:end));
    end
end
