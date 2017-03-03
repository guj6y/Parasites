clear all
close all

S = 40;
sList = 1:S;
C = .15;
nBasal = 5;
nZ = 8;
h=1.2;
ZAll = logspace(-2,5,nZ);
ax = .314;

spTypes = {'Inv','EVe'};


SList = 1:S;

%The abundance level web directory
dataDir =...
    sprintf('/Volumes/Buster/Research/Parasitism/Dynamics/JulySimulations');
resultsDir = sprintf('%s/ResultsTf10000',dataDir);


codeDir = pwd;

nWeb = 100;
Tfs = [2000,5000,10000,20000,40000];
nT = length(Tfs);

nZ = 3;
ZAll = logspace(1,3,nZ);







%The abundance level web directory
dataDirS = sprintf('%s/S%03u',dataDir,S);

%The abundance,connectance level web directory
dataDirSC = sprintf('%s/C%03u',dataDirS,round(100*C));


try
    load('persistenceAllNullModel.mat');
    load('ExtctDataNullModel.mat');
catch
    tic
    extctData = zeros(0,5);
    persistenceAll = cell(2,1);
    [persistenceAll{:}] = deal(zeros(nWeb,nZ));
    for ii = 1:nWeb
        for type = 1:2
        count = 0;
        
        linkList = csvread(sprintf('%s/web%04u.csv',dataDirSC,ii-1));
        finals = csvread(sprintf('%s/finals%s%04u.csv',resultsDir,spTypes{type},ii-1));
        means = csvread(sprintf('%s/means%s%04u.csv',resultsDir,spTypes{type},ii-1));
        stds = csvread(sprintf('%s/stds%s%04u.csv',resultsDir,spTypes{type},ii-1));
        
        res = linkList(:,1);
        con = linkList(:,2);
        
        properties = calculateLocalProperties(res,con);
        swtl = properties(:,10);
        
        fprintf('Loading Z = ')
        for jj = 1:nZ
            sol = load(sprintf('%s/sol%s%04uZ%u.mat',resultsDir,spTypes{type},ii-1,jj));
            curTime = toc;
            clc
            fprintf('Loading data for %s web %u, Z = %.1e\n',spTypes{type},ii,ZAll(jj))
            fprintf('Approximately %.2f minutes have elapsed.\n',curTime/60)
            numLoaded = (ii-1)*6 + jj -1;
            fprintf('%u out of 300 solution structures loaded.\n',numLoaded)
            %     x = 0:40000;
            %     y = zeros(40,40001);
            %     curX = 0;
            
            for	kk = 1:sol.n;
                count = count+1;
                curPart = sprintf('part%u',kk);
                tExtct = sol.(curPart).xe;
                indxExtct = sol.(curPart).ye<1e-30;
                xExtct = .314*(ZAll(jj).^(swtl(indxExtct))).^(-0.25);
                sExtct = SList(indxExtct)';
                catchMult = ones(size(sExtct));
                extctData = [extctData;
                    ii*catchMult, log10(ZAll(jj))*catchMult, sExtct, xExtct,tExtct*catchMult];
                
                kkX0 = ceil(sol.(curPart).x(1));
                kkXEnd = floor(sol.(curPart).x(end));
                idx0 = kkX0 + 1;
                idxEnd = kkXEnd + 1;
                
            end
            persistenceAll{type}(ii,jj) = mean(sol.(curPart).y(:,end)>1e-30);


            
            
        end
        end
    end
    save('persistenceAllNullModel','persistenceAll');
    save('extctDataNullModel','extctData');
end

Z1 = extctData(:,2) == -2;
Z2 = extctData(:,2) == -1;
Z3 = extctData(:,2) == 0;
Z4 = extctData(:,2) == 1;
Z5 = extctData(:,2) == 2;
Z6 = extctData(:,2) == 3;
Z7 = extctData(:,2) == 4;
Z8 = extctData(:,2) == 5;
figure;
hold on
for ii = 1:2
    
    plot(log10(ZAll),mean(persistenceAll{ii}),'-o')
    
end
hold off
% markerColors = [1.0 0.0 0.0; %Red
%     1.0 0.5 0.0; %Orange
%     1.0 0.8 0.0; %Yellow
%     0.0 1.0 0.0; %Green
%     0.0 0.0 1.0; %Blue
%     0.3 0.0 0.5; %Indigo
%     1.0 0.0 1.0; %Violet
%     0.0 0.0 0.0]; %Black

gscatterGroupsNames = {'Z = 10^{-2}';...
    'z = 10^{-1}';...
    'Z = 10^0';...
    'Z = 10^1';...
    'Z = 10^2';...
    'Z = 10^3';...
    'Z = 10^4';...
    'Z = 10^5'};

gscatterGroups = gscatterGroupsNames(extctData(:,2)+3);

h1 = figure;
hold on
gscatter(log10(extctData(:,4)),log10(extctData(:,5)),...
    gscatterGroups)
title('Extinction Time and Metabolic Rate grouped by consumer-resource Body Size Ratio')
xlabel('log_{10}(x_i)')
ylabel('log_{10}(T_{extct})')

saveas(h1,'metabolicRateAndExtinctionTimeImpreciseLowIC.pdf');

[persistenceZ5Sorted, idxZ5]= sort(persistenceAll{1}(:,3));
web25 = idxZ5(25);
web50 = idxZ5(50);
web75 = idxZ5(75);

sol25 = load(sprintf('%s/sol%s%04uZ%u.mat',resultsDir,spTypes{1},web25-1,3));
sol50 = load(sprintf('%s/sol%s%04uZ%u.mat',resultsDir,spTypes{1},web50-1,3));
sol75 = load(sprintf('%s/sol%s%04uZ%u.mat',resultsDir,spTypes{1},web75-1,3));
Z = 100;
count = 0;
x25 = [];
y25 = zeros(40,0);
extctSp25 = [];
linkList25 = csvread(sprintf('%s/web%04u.csv',dataDirSC,web25));
properties25 = calculateLocalProperties(linkList25(:,1),linkList25(:,2));
swtl25 = properties25(:,10);
xi25 = .314*(Z.^swtl25).^(-.25);
for	kk = 1:sol25.n;
    count = count+1;
    curPart = sprintf('part%u',kk);
    x25 = [x25 sol25.(curPart).x(2:end)];
    y25 = [y25 sol25.(curPart).y(:,2:end)];
    extctSp25 = [extctSp25; sList(sol25.(curPart).ye<1e-30)'];
end

x50 = [];
y50 = zeros(40,0);
extctSp50 = [];
linkList50 = csvread(sprintf('%s/web%04u.csv',dataDirSC,web50));
properties50 = calculateLocalProperties(linkList50(:,1),linkList50(:,2));
swtl50 = properties50(:,10);
xi50 = .314*(Z.^swtl50).^(-.25);

for	kk = 1:sol50.n;
    count = count+1;
    curPart = sprintf('part%u',kk);
    x50 = [x50 sol50.(curPart).x(2:end)];
    y50 = [y50 sol50.(curPart).y(:,2:end)];
    
    extctSp50 = [extctSp50; sList(sol50.(curPart).ye<1e-30)'];
end

x75 = [];
y75 = zeros(40,0);
extctSp75 = [];
linkList75 = csvread(sprintf('%s/web%04u.csv',dataDirSC,web75));
properties75 = calculateLocalProperties(linkList75(:,1),linkList75(:,2));
swtl75 = properties75(:,10);
xi75 = .314*(Z.^swtl75).^(-.25);
for	kk = 1:sol75.n;
    count = count+1;
    curPart = sprintf('part%u',kk);
    x75 = [x75 sol75.(curPart).x(2:end)];
    y75 = [y75 sol75.(curPart).y(:,2:end)];
    extctSp75 = [extctSp75; sList(sol75.(curPart).ye<1e-30)'];
    
end

webCheck = 50;
web = idxZ5(webCheck);

sol = load(sprintf('%s/sol%s%04uZ%u.mat',resultsDir,spTypes{1},web-1,3));
x = [];
y = zeros(40,0);
extctSp = [];
linkList = csvread(sprintf('%s/web%04u.csv',dataDirSC,web));
properties = calculateLocalProperties(linkList(:,1),linkList(:,2));
swtl = properties(:,10);
xi = .314*(Z.^swtl).^(-.25);
for	kk = 1:sol.n;
    count = count+1;
    curPart = sprintf('part%u',kk);
    x = [x sol.(curPart).x(2:end)];
    y = [y sol.(curPart).y(:,2:end)];
    
end
h2 = figure;
plot(x,log10(y))
title('Example of Dynamics')
xlabel('Time')
ylabel('log_{10}(Biomass)')
saveas(h2,'ExampleDynamics50_Z8ImpreciseLowIC.pdf')
extctSp = sList(sol.(curPart).y(:,end)<1e-30);

numExtct = length(extctSp);

extctCompare = zeros(numExtct,2);
figure;
hold on
for ii = 1:numExtct
    
    
    yExtctii = y(extctSp(ii),:);
    yExtctii = yExtctii(yExtctii>0);
    xExtctii = x(yExtctii>0);
    xii = xi(extctSp(ii));
    extctCompare(ii,1) = xii;
    extctCompare(ii,2) = (log10(yExtctii(end)) - log10(yExtctii(floor(.75*end))))/(xExtctii(end) - xExtctii(floor(.75*end)));
    
    %plot(xExtctii,log10(yExtctii),xExtctii,log10(yExtctii(end)) - xii*(xExtctii-xExtctii(end)))
end

%plot(x,log10(y))
