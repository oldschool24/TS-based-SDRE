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

    addpath('../')
    [q1, q2, z1, z2] = ndgrid(q1Range, q2Range, z1Range, z2Range);
    nTests = numel(q1);
    idxExamples = randsample(nTests, nExamples);
    Q = q * eye(8);
    R = r * eye(2);
    if isAnalyze
        tsX = [];
        predX = [];
        f_true = [];
        f_pred = []; 
        B_true = []; 
        B_pred = [];
        analyzeDir = fullfile(expPath, "testAnalysis");
        if ~exist(analyzeDir, 'dir')
            mkdir(analyzeDir)
        end
    end

    criterion = zeros(nTests, 8+6+8);
    warning('off', 'fuzzy:general:warnEvalfis_NoRuleFired')
    warning('off', 'fuzzy:general:diagEvalfis_OutOfRangeInput')
    tic
    for k=1:nTests
        x0 = [q1(k); q2(k); z1(k); z2(k); 0; 0; 0; 0];
        if ismember(k, idxExamples)
            imgDir = fullfile( ...
                expPath, ['example_' testConfigName '(' num2str(k) ')']);
            if ~exist(imgDir, "dir")
                mkdir(imgDir)
            end
        else
            imgDir = '';
        end
        simStats = mainSim(modelPath, 'flex2link', dt, T, x0, Q, R, ...
                           @ode15s, isWrap, imgDir, known, isAnalyze);
        mainStats = [x0', simStats.tsCriterion, simStats.sdreCriterion, ...
                     simStats.tsTime, simStats.sdreTime, ...
                     simStats.tsWallTime, simStats.sdreWallTime, ...
                     simStats.stopPoint];
        for iStats=1:22
            criterion(k, iStats) = mainStats(iStats);
        end
        if isAnalyze
            tsX = [tsX; simStats.tsX]; 
            predX = [predX; simStats.predX]; 
            f_true = [f_true; simStats.f_true]; 
            f_pred = [f_pred; simStats.f_pred]; 
            B_true = [B_true; simStats.B_true]; 
            B_pred = [B_pred; simStats.B_pred]; 
        end
    end
    
    save(fullfile(expPath, testConfigName), 'criterion')
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
