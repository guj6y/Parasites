clear all
close all

S = 40;
sList = 1:S;
C = .15;
nBasal = 5;
ax = .314;



SList = 1:S;

%The abundance level web directory
dataDir =...
    sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/Brose2006Webs');
resultsDir = sprintf('%s/ResultsTf1000000',dataDir);


codeDir = pwd;

nWeb = 91;

nZ = 8;
nSnapshots = 8;
ZAll = logspace(-2,5,nZ);


nThresh = 7;
threshAll = logspace(-30,-6,7);



%The abundance level web directory
dataDirS = sprintf('%s/S%03u',dataDir,S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));

persistenceAtSnapshots = zeros(nWeb,nZ,nSnapshots);
totalExtinctionsAtSnapshots = zeros(nWeb,nZ,nSnapshots);

tic
extctDataAll = zeros(0,5);
for ii = 1:nWeb
    
    count = 0;
    
    for Z = 1:nZ
        
        
        
        snapshots = csvread(sprintf('%s/snapshots%04uZ%u.csv',resultsDir,ii-1,Z));
        persistenceAtSnapshots(ii,Z,:) = mean(snapshots>1e-30);
        totalExtinctionsAtSnapshot = sum(snapshots(:,end)<1e-30);
        
        extctData = csvread(sprintf('%s/extinctions%04uZ%u.csv',resultsDir,ii-1,Z));
        extctIds = extctData(:,2);
        extinctionsIiZ = zeros(totalExtinctionsAtSnapshot,5);
        count = 1;
        while numel(extctIds)>0
            extinctionsIiZ(count,:) = [extctData(1,:),Z,ii];
            extctData(extctIds==extctIds(1),:) = [];
            extctIds(extctIds==extctIds(1)) = [];
            count = count+1;
        end
        
        extctDataAll = [extctDataAll;extinctionsIiZ];    
    end
    
end


%h
h1 = figure;
stdPersistences = zeros(nSnapshots,nZ);
meanPersistences = zeros(nSnapshots,nZ);



meanPersistences(:) = mean(persistenceAtSnapshots);
stdPersistences(:) = std(persistenceAtSnapshots);
tcrit = tinv(.975,nWeb-1);
err = tcrit*stdPersistences/sqrt(nWeb);
x = repmat(log10(ZAll)',1,8);
hold on
plot(x(:,1:7),meanPersistences(:,1:7),'-o')
h1e = errorbar(x(:,8),meanPersistences(:,8),err(:,8),'-o');
hold off
title('Persistence vs log(Z) for different Final Times')
xlabel('log_{10}(Z)')
ylabel('Fraction of Species with B_i > 10^{-30}')
legend('Tf = 5000',...
    'Tf = 10000',...
    'Tf = 25000',...
    'Tf = 50000',...
    'Tf = 100000',...
    'Tf = 250000',...
    'Tf = 500000',...
    'Tf = 1000000',...
    'Location','SouthEast');
saveas(h1,'BWM2006DeepTimeGraph.jpg')

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