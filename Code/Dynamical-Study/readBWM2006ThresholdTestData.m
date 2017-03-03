clear all
close all

S = 40;
sList = 1:S;
C = .15;
nBasal = 5;
ax = .314;

nWeb = 91;

SList = 1:S;
nZ = 8;
ZAll = logspace(-2,5,nZ);
%The abundance level web directory
dataDir =...
    sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/Brose2006Webs');

codeDir = pwd;

%The abundance level web directory
dataDirS = sprintf('%s/S%03u',dataDir,S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));







nThresh = 7;
threshAll = logspace(-30,-6,7);





plotThis = zeros(nThresh,nZ);
persistences = zeros(nWeb,nZ);
persistencesDeep = zeros(nWeb,nZ);
stdPersistences = zeros(nThresh,nZ);
for jj = 1:nThresh

extctThresh = threshAll(jj);
resultsDir = sprintf('%s/ResultsThresh%u',dataDir,extctThresh);
resultsDirDeep = sprintf('%s/ResultsTf1000000',dataDir);
tic
extctDataAll = zeros(0,5);
for ii = 1:nWeb
    
   
        meanDeep = csvread(sprintf('%s/means%04u.csv',resultsDirDeep,ii-1));
        means = csvread(sprintf('%s/means%04uTh%u.csv',resultsDir,ii-1,jj));
        persistences(ii,:) = mean(means>0);
        persistencesDeep(ii,:) = mean(meanDeep>0);
    
end
plotThisDeep(jj,:) = mean(persistencesDeep);
stdDeep(jj,:) = mean(persistencesDeep);
plotThis(jj,:) = mean(persistences);
stdPersistences(jj,:) = std(persistences);
end

%
h1 = figure;

tcrit = tinv(.975,nWeb-1);
err = tcrit*stdPersistences/sqrt(nWeb);
errDeep = tcrit*stdDeep/sqrt(nWeb);

hold on
plot(log10(ZAll),plotThis(1:5,:)h,'-o')
h1e = errorbar(log10([ZAll;ZAll]'),plotThis(6:7,:)',err(6:7,:)','-o');
h2e = plot(log10(ZAll),plotThisDeep(1,:),'k-x');
hold off
title(sprintf('Persistence vs log(Z) for different Extinction Thresholds\nTf=5000'));
xlabel('log_{10}(Z)')
ylabel('Fraction of Species with B_i > B_e')

legend('B_e = 10^{-30}',...
    'B_e = 10^{-26}',...
    'B_e = 10^{-22}',...
    'B_e = 10^{-18}',...
    'B_e = 10^{-14}',...
    'B_e = 10^{-10}',...
    'B_e = 10^{-6}',...
    'B_e = 10^{-30},Tf=10^6',...
    'Location','NorthWest');
saveas(h1,'BWM2006ThresholdGraph.jpg')

%
% h2 = figure;
% for ii = 1:8
%     subplot(2,4,ii);
%     hist(persistenceAtSnapshots(:,ii,8))
%     title(sprintf('Persistence for Z = %u',ii))
% end
gscatterGroupsNames = {'Z = 10^{-2}';...
    'z = 10^{-1}';...
    'Z = 10^0';...
    'Z = 10^1';...
    'Z = 10^2';...
    'Z = 10^3';...
    'Z = 10^4';...
    'Z = 10^5'};
gscatterGroups = gscatterGroupsNames(extctDataAll(:,4));

 markerColors = [1.0 0.0 0.0; %Red
    0.3 0.0 0.5; %Indigo
    1.0 0.8 0.0; %Yellow
    1.0 0.5 0.0; %Orange
    0.0 0.0 1.0; %Blue
    0.0 1.0 0.0; %Green
    1.0 0.0 1.0; %Violet
    0.0 0.0 0.0]; %Black
h3 = figure;
gscatter(log10(extctDataAll(:,3)),log10(extctDataAll(:,1)),...
    gscatterGroups,markerColors,'.',6)
title(sprintf('Extinction Time and Metabolic Rate grouped by consumer-resource Body Size Ratio'))
xlabel('log_{10}(x_i)')
ylabel('log_{10}(T_{extct})')
saveas(h3,'figure1','pdf')

h4 = figure;
k40Cutoff = extctDataAll(:,1)<=40000;
gscatter(log10(extctDataAll(k40Cutoff,3)),log10(extctDataAll(k40Cutoff,1)),gscatterGroups(k40Cutoff))
%}