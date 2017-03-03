%webGeneration;
%Webs:
% 1: Bahia San Quentin - Baja California
% 2: Carpinteria Salt Marsh - SoCal (north west of LA)
% 3: Punta de la Banda - Baja California
% 4: Flensburg Fjord - Germany
% 5: Otago Harbor - New Zealand
% 6: Sylt Tidal Basin - Germany/Denmark
propertiesAll = [propertiesCell{1,1};
    propertiesCell{2,1};
    propertiesCell{3,1};
    propertiesCell{4,1};
    propertiesCell{5,1};
    propertiesCell{6,1}];

paraAll = [propertiesCell{1,2};
    propertiesCell{2,2};
    propertiesCell{3,2};
    propertiesCell{4,2};
    propertiesCell{5,2};
    propertiesCell{6,2}];

classAll = [speciesTypeCell{1};
    speciesTypeCell{2};
    speciesTypeCell{3};
    speciesTypeCell{4};
    speciesTypeCell{5};
    speciesTypeCell{6}];

classAll = classAll + 3*paraAll;

nBSQ = length(propertiesCell{1,2});
nParBSQ = sum(propertiesCell{1,2});
nCSM = length(propertiesCell{2,2});
nParCSM = sum(propertiesCell{2,2});
nPDB = length(propertiesCell{3,2});
nParPDB = sum(propertiesCell{3,2});
nFLEN = length(propertiesCell{4,2});
nParFLEN = sum(propertiesCell{4,2});
nOTAG = length(propertiesCell{5,2});
nParOTAG = sum(propertiesCell{5,2});
nSYLT = length(propertiesCell{6,2});
nParSYLT = sum(propertiesCell{6,2});

nSpecies = [nBSQ, nCSM, nPDB, nFLEN,nOTAG,nSYLT];
nPar = [nParBSQ, nParCSM, nParPDB, nParFLEN,nParOTAG,nParSYLT];

%1
BSQList = 1:nBSQ;
%2
CSMList = (1:nCSM) + nBSQ;
%3
PDBList = (1:nPDB) + nBSQ + nCSM;
%4
FLENList = (1:nFLEN) + nBSQ + nCSM + nPDB;
%5
OTAGList = (1:nOTAG) + nBSQ + nCSM + nPDB + nFLEN;
%6
SYLTList = (1:nSYLT) + nBSQ + nCSM + nPDB + nFLEN + nOTAG;

trainCell = cell(7,1);
testCell = cell(7,1);

longList = 1:(nBSQ + nCSM + nPDB + nFLEN + nOTAG + nSYLT);



trainCell{2} = [BSQList, CSMList,PDBList];
testCell{2} = [FLENList, OTAGList, SYLTList];

trainCell{3} = FLENList;
testCell{3} = SYLTList;

trainCell{4} = [BSQList, FLENList];
testCell{4} = OTAGList;

trainCell{5} = [BSQList, OTAGList];
testCell{5} = FLENList;

trainCell{6} = [FLENList, OTAGList];
testCell{6} = BSQList;

trainCell{7} = [BSQList, CSMList];
testCell{7} = PDBList;

costs  = 1:10:100;
N = 100;
nCosts = length(costs);
truePositives = zeros(N,nCosts);
falsePositives = zeros(N,nCosts);
falseNegatives = zeros(N,nCosts);
trueNegatives = zeros(N,nCosts);
nParaAll = zeros(N,length(costs));
%Could try bootstrapping.



S = 156;
C = 0.0889;
trial = 1;
nCost = 0;
%                1         2    3         4             5           6               7           8               9           10      11
varNames = {'clustCoef','gen','vul','meanVulPrey','meanImpPrey','meanGenPred','meanImpPred','minSPToBasal','numConnBasal','SWTL','inLoop'};
predictorsUsed = [1;
    2;
    3;
    4;
    ...5;
    6;
    ...7;
    ...8;
    9;
    10;
    11
    ];
N=1;
trees = cell(N,1);
errors = zeros(N,8);
cost =1 ;
for ii=1
    for kk = 1:N
        %trainCell{ii} = datasample(1:length(paraAll(classAll~=0)),ceil(.75*length(paraAll(classAll~=0))),'replace',false);
        %binTrain1 = false(length(longList(classAll~=0)),1);
        %binTrain1(trainCell{1}) = true;
        %testCell{1} = longList(~binTrain1);
        trees{kk} = fitctree(propertiesAll(classAll~=0,predictorsUsed),classAll(classAll~=0)...
            ...,'Cost',[0 1;cost 0]...
            ,'PruneCriterion','impurity'...
            ,'MaxNumSplits',7 ...,
            ,'PredictorNames',varNames(predictorsUsed));
        preds = predict(trees{ii},propertiesAll(testCell{ii},predictorsUsed));
        truth = paraAll(testCell{ii});
        n = length(preds);
        
        err = sum(preds~=truth)/n;
        prior = sum(truth)/n;
        truePos = sum(preds&truth)/sum(truth);
        falsePos = sum(preds&(~truth))/sum(~truth);
        trueNeg = sum((~preds)&(~truth))/sum(~truth);
        falseNeg = sum((~preds)&truth)/sum(truth);
        
        close all hidden
        
        errors(kk,:) = [err,prior,truePos,falsePos,falseNeg,trueNeg,sum(truth),sum(preds)];
    end
    %{
    for trial = 1:N
        speciesList = 1:S;
        [res,con] = NicheModel_nk(S,C);
        propertiesTest = calculateLocalProperties(res,con);
        
        
        [preds, post, ~] = predict(nb,propertiesTest(:,predictorsUsed));
        [postSort, idxSort] = sort(post(:,2));
        speciesList = speciesList(idxSort);
        nParaAll(trial,nCost) = sum(preds);
        
    end
    %}
end
%The idea is that I want to have about the same false negatives as false
%positives so that I can get about the right number of parasites.
