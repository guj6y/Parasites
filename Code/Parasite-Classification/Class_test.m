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

nBSQ = length(propertiesCell{1,2});
nCSM = length(propertiesCell{2,2});
nPDB = length(propertiesCell{3,2});
nFLEN = length(propertiesCell{4,2});
nOTAG = length(propertiesCell{5,2});
nSYLT = length(propertiesCell{6,2});

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

train{1} = datasample(1:length(paraAll),ceil(.75*length(paraAll)),'replace',false);
binTrain1 = false(1,length(longList));
binTrain1(train{1}) = true;
test{1} = longList(~binTrain1);

train{2} = [BSQList, CSMList,PDBList];
test{2} = [FLENList, OTAGList, SYLTList];

train{3} = FLENList;
test{3} = SYLTList;

train{4} = [BSQList, FLENList];
test{4} = OTAGList;

train{5} = [BSQList, OTAGList];
test{5} = FLENList;

train{6} = [FLENList, OTAGList];
test{6} = BSQList;

train{7} = [BSQList, CSMList];
test{7} = PDBList;

costs  = 0.5:0.5:5.5;
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
                  5;
                  6;
                  7;
                  ...8;
                  ...9;
                  10;
                  ...11
                  ];

for cost=costs
    nCost = nCost+1;
    %train = datasample(1:length(paraAll),ceil(.75*length(paraAll)),'replace',false);
    n = length(paraAll);
    %train = (1:(nBSQ+nPDB+nCSM));
    binTrain = false(n,1);
    binTest = binTrain;
    binTrain(train) = true;
    %binTest(train) = true;
    %binTest = ~binTest;
    binTest(test) = true;
    nb = fitctree(propertiesAll(binTrain,predictorsUsed),paraAll(binTrain),...
        'Cost',[0 1;cost 0]...
        ,'PruneCriterion','impurity'...
        ,'MaxNumSplits',4 ...,
        ,'PredictorNames',varNames(predictorsUsed));
    %preds = predict(tree,propertiesAll(binTest,1:11));
    
    
    %
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
%mean(falsePositives)-mean(falseNegatives)
%sum(paraAll(binTest))
%length(preds)-sum(paraAll(binTest))
%sum([truePositives;falsePositives])
%sum([trueNegatives;falseNegatives])