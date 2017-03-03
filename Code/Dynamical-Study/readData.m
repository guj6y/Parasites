%Reading saved data from parasiteNumericalExperiments to analyze results.

S = 50;
C = 0.1;
SList = 1:S;

%The abundance level web directory
dataDirS = sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/NMWebs/S%03u',S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));


codeDir = pwd;

nFPar = 10;
fParAll = linspace(0,1,nFPar);
nWeb = 100;
nParAll = zeros(nWeb,nFPar);

model = 0;
persistenceAll0 = zeros(nWeb,nFPar);
persistenceCon0 = zeros(nWeb,nFPar);
persistencePar0 = zeros(nWeb,nFPar);
persistenceFree0 = zeros(nWeb,nFPar);
totBiomass0 = zeros(nWeb,nFPar);
resultsDir = sprintf('%s/Results%u',dataDirSC,model);

for ii = 1:nWeb
    
    nodeData = csvread(sprintf('%s/nodeData%04u.csv',dataDirSC,ii-1));
    basal = nodeData(:,7) > 0;
    finals = csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1));
    means = csvread(sprintf('%s/means%04u.csv',resultsDir,ii-1));
    stds = csvread(sprintf('%s/stds%04u.csv',resultsDir,ii-1));
    
    nParAll(ii,:) = round(50-sum(basal))*fParAll;
    
    persistenceAll0(ii,:) = sum(finals>0);
    persistenceCon0(ii,:) = sum(finals(~basal,:)>0);
    totBiomass0 = sum(finals>0);
    
end

model = 1;
persistenceAll1 = zeros(nWeb,nFPar);
persistenceCon1 = zeros(nWeb,nFPar);
persistencePar1 = zeros(nWeb,nFPar);
persistenceFree1 = zeros(nWeb,nFPar);
totBiomass1 = zeros(nWeb,nFPar);
resultsDir = sprintf('%s/Results%u',dataDirSC,model);

for ii = 1:nWeb
    
    nodeData = csvread(sprintf('%s/nodeData%04u.csv',dataDirSC,ii-1));
    basal = nodeData(:,7) > 0;
    finals = csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1));
    means = csvread(sprintf('%s/means%04u.csv',resultsDir,ii-1));
    stds = csvread(sprintf('%s/stds%04u.csv',resultsDir,ii-1));
    
    
    persistenceAll1(ii,:) = sum(finals>0);
    persistenceCon1(ii,:) = sum(finals(~basal,:)>0);
    totBiomass1 = sum(finals>0);
    
end

model = 2;
persistenceAll2 = zeros(nWeb,nFPar);
persistenceCon2 = zeros(nWeb,nFPar);
persistencePar2 = zeros(nWeb,nFPar);
persistenceFree2 = zeros(nWeb,nFPar);
totBiomass2 = zeros(nWeb,nFPar);
resultsDir = sprintf('%s/Results%u',dataDirSC,model);

for ii = 1:nWeb
    
    nodeData = csvread(sprintf('%s/nodeData%04u.csv',dataDirSC,ii-1));
    basal = nodeData(:,7) > 0;
    finals = csvread(sprintf('%s/finals%04u.csv',resultsDir,ii-1));
    means = csvread(sprintf('%s/means%04u.csv',resultsDir,ii-1));
    stds = csvread(sprintf('%s/stds%04u.csv',resultsDir,ii-1));
    
    
    persistenceAll2(ii,:) = sum(finals>0);
    persistenceCon2(ii,:) = sum(finals(~basal,:)>0);
    totBiomass2 = sum(finals>0);
    
end
close all
figure;
plot(fParAll,mean(persistenceAll0)/50,'b*'...
    ,fParAll,mean(persistenceAll1)/50,'ro'...
    ,fParAll,mean(persistenceAll2)/50,'gs')
title('Fraction of persistence')
xlabel('Fraction of Parasites')
ylabel('% persistence')
legend('Model 0','Model 1','Model 2')

figure;
plot(fParAll,mean(persistenceCon0)/50,'b*'...
    ,fParAll,mean(persistenceCon1)/50,'ro'...
    ,fParAll,mean(persistenceCon2)/50,'gs')
title('Fraction of Consumers persistence')
xlabel('Fraction of Parasites')
ylabel('% persistence')
legend('Model 0','Model 1','Model 2')