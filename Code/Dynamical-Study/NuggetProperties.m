%This analyzes web nuggets.
clear all
close all
Z = 100;

dataDir = '~/Desktop/JulySimulations/WebNuggets';

OGWebDir = sprintf('%s/OGWebs',dataDir);
ZWebDir = sprintf('%s/Z10Nuggets',dataDir);


cutoff = 30;

%Choosing the final S value.
S = 40;


OGIndex = csvread(sprintf('%s/webIndex.csv',OGWebDir),1,0);


deadWebCount = 0;


for Z = [10]
    count = 0;
    
    genAll = [];
    vulAll = [];
    genAll0 = [];
    vulAll0 = [];
    extctAll = [];
    cAll = [];
    rAll = [];
    nAll = [];
    
    medianPATLPreyDiffAll= [];
    medianSWTLPreyDiffAll = [];
    medianPATLPredDiffAll = [];
    medianSWTLPredDiffAll = [];
    
    meanGenPredAll0 = [];
    meanVulPreyAll0 = [];
    
    maxVulPreyAll = [];
    minVulPreyAll = [];
    
    maxGenPredAll = [];
    minGenPredAll = [];
    
    meanImpPreyAll0 = [];
    meanImpPredAll0 = [];
    spToBasalAll = [];
    nBasalConAll = [];
    
    genAll1 = [];
    vulAll1 = [];
    
    
    meanImpPredAll1 = [];
    meanImpPreyAll1 = [];
    
    meanGenPredAll1 = [];
    meanVulPreyAll1 = [];
    
    maxGenPredAll1 = [];
    minGenPredAll1 = [];
    
    maxVulPreyAll1 = [];
    minVulPreyAll1 = [];
    
    spToBasalAll1 = [];
    nBasalCon = [];
    
    meanBiomassPrey = zeros(80,0);
    extntAll = zeros(80,0);
    extctOrderAll = [];
    BAll = [];
    patlAll = [];
    fBasalAll = [];
    persistenceAll = [];
    
    totB  = [];
    badWebs = [];
    
    patlDiff = [];
    
    patlAll0 = [];
    patlAll = [];
    simAll = [];
    maxSimPredAll = [];
    maxSimPreyAll = [];
    ZWebDir = sprintf('%s/Z%uNuggets',dataDir,Z);
    ZWebIndex = csvread(sprintf('%s/webIndex.csv',ZWebDir),1,0);
    webNumbersZ = ZWebIndex(:,1);
    sims = zeros(size(webNumbersZ));
    sims0 = zeros(size(webNumbersZ));
    
    
    deadWebs = zeros(size(webNumbersZ));
    
    for web = webNumbersZ'
        %Testing attractivity of final states.
        web
        count = count + 1;
        S0 = ZWebIndex(count,2);
        S = ZWebIndex(count,5);
        
        sListOld = 1:S0;
        sListNew = 1:S;
        sRename = sListOld;
        nExtinct =  S0-S;
        
        OGNodeData = csvread(sprintf('%s/nodeData%06u.csv',OGWebDir,web));
        nodeData = csvread(sprintf('%s/nodeData%06u.csv',ZWebDir,web));
        linkData = csvread(sprintf('%s/linkData%06u.csv',ZWebDir,web));
        
        
        
        gr0 = OGNodeData(:,4);
        n0 = OGNodeData(:,1);
        c0 = OGNodeData(:,2);
        r0 = OGNodeData(:,3);
        
        alive = isinf(nodeData(:,5));
        
        n = n0(alive);
        c = c0(alive);
        r = r0(alive);
        
        res0 = linkData(:,1);
        con0 = linkData(:,2);
        L = sum(~linkData(:,3));
        S = numel(n);
        adjMx = sparse(res0,con0,1,S0,S0);
        if L/S^2 < .01
            badWebs = [badWebs;count];
            continue
        end
        try
            [res,con,~,~,~,~] = AddTrophicConsumers(n,c,r,0,L/S^2);
        catch
            badWebs = [badWebs;count];
            continue
        end
        extctOrder = nodeData(:,5);
        relExtctOrder = extctOrder/max(extctOrder(isfinite(extctOrder)));
        extctOrderAll = [extctOrderAll; extctOrder];
        
        localProperties = calculateLocalProperties(res,con,S);
        localProperties0 = calculateLocalProperties(res0,con0,80);
        
        sim0 = calculateSimilarity(res0,con0);
        simStructure = sim0.*adjMx;
        maxSimPred = max(simStructure')';
        maxSimPrey = max(simStructure)';
        
        maxSimPredAll = [maxSimPredAll;maxSimPred];
        maxSimPreyAll = [maxSimPreyAll;maxSimPrey];
        simAll = [simAll;max(sim0)'];
        
        swtl = localProperties0(:,10);
        sptb = localProperties0(:,8)+1;
        patl = 2*swtl - sptb;
        
        B = nodeData(:,4);
        BAll = [BAll;B/sum(B)];
        
        A = sparse(res0,con0,1,80,80);
        
        totalPreyBiomass(:,count) = A'*B.^1.2;
        extntAll(:,count) = B>0;
        
        totB = [totB;sum(B)];
        patlAll = [patlAll;patl];
        
        fBasal = sum(swtl==1)/80;
        persistence = sum(B>0)/80;
        
        fBasalAll = [fBasalAll;fBasal];
        persistenceAll = [persistenceAll;persistence];
        
        %         extctList = nodeData(:,end);
        %
        %         disconnected = false;
        
        vul = localProperties(:,3)*(length(res)/S);
        gen = localProperties(:,2)*(length(res)/S);
        
        vulAll = [vulAll;vul];
        genAll = [genAll;gen];
        
        vul0 = localProperties0(:,3)*(length(res0)/80);
        gen0 = localProperties0(:,2)*(length(res0)/80);
        
        
        
        meanVulPrey0 = localProperties0(:,4);
        meanVulPrey0(isnan(meanVulPrey0)) = 0;
        
        meanImpPrey0 = localProperties0(:,5);
        
        meanGenPred0 = localProperties0(:,6);
        meanGenPred0(isnan(meanGenPred0)) = 0;
        
        meanImpPred0 = localProperties0(:,7);
        
        spToBasal0 = localProperties0(:,8);
        
        nBasalCon0 = localProperties(:,9);
        
        maxGenPred0 = localProperties0(:,12);
        minGenPred0 = localProperties0(:,13);
        maxVulPrey0 = localProperties0(:,14);
        minVulPrey0 = localProperties0(:,15);
        
        median_patl_prey_diff0 = localProperties0(:,16);
        median_patl_pred_diff0 = localProperties0(:,17);
        median_swtl_prey_diff0 = localProperties0(:,18);
        median_swtl_pred_diff0 = localProperties0(:,19);
        
        medianPATLPreyDiffAll= [medianPATLPreyDiffAll; median_patl_prey_diff0];
        medianSWTLPreyDiffAll = [medianSWTLPreyDiffAll; median_patl_pred_diff0];
        medianPATLPredDiffAll = [medianPATLPredDiffAll; median_swtl_prey_diff0];
        medianSWTLPredDiffAll = [medianSWTLPredDiffAll; median_swtl_pred_diff0];
        
        
        gen0 = localProperties0(:,2)*(length(res0)/80);
        vul0 = localProperties0(:,3)*(length(res0)/80);
        
        
        
        genAll0 = [genAll0;gen0];
        vulAll0 = [vulAll0;vul0];
        
        
        meanImpPredAll0 = [meanImpPredAll0; meanImpPred0];
        meanImpPreyAll0 = [meanImpPreyAll0; meanImpPrey0];
        
        meanGenPredAll0 = [meanGenPredAll0; meanGenPred0];
        meanVulPreyAll0 = [meanVulPreyAll0; meanVulPrey0];
        
        maxGenPredAll = [maxGenPredAll;maxGenPred0];
        minGenPredAll = [minGenPredAll;minGenPred0];
        
        maxVulPreyAll = [maxVulPreyAll;maxVulPrey0];
        minVulPreyAll = [minVulPreyAll;minVulPrey0];
        
        spToBasalAll = [spToBasalAll; spToBasal0];
        
        res=res0;
        con=con0;
        notExtinct = true(80,1);
        for ii = 1:nExtinct
            extctSp = find(extctOrder==ii);
            
            if sum(extctOrder==ii)>1
                fprintf('uhoh\n')
            end
            notExtinct(extctSp) = false;
            
            deadLinks = (res==extctSp)|(con==extctSp);
            res(deadLinks) = [];
            con(deadLinks) = [];
            
            
            localProperties0 = calculateLocalProperties(res,con,80);
            
            gen0(notExtinct) = localProperties0(notExtinct,2)*...
                (length(res)/(80));
            vul0(notExtinct) = localProperties0(notExtinct,3)*...
                (length(res)/(80));
            
            meanVulPrey0(notExtinct) = localProperties0(notExtinct,4);
            meanVulPrey0(isnan(meanVulPrey0)) = 0;
            
            meanImpPrey0(notExtinct) = localProperties0(notExtinct,5);
            
            meanGenPred0(notExtinct) = localProperties0(notExtinct,6);
            meanGenPred0(isnan(meanGenPred0)) = 0;
            
            meanImpPred0(notExtinct) = localProperties0(notExtinct,7);
            
            spToBasal0(notExtinct) = localProperties0(notExtinct,8);
            
            nBasalCon0(notExtinct) = localProperties0(notExtinct,9);
            
            maxGenPred0(notExtinct) = localProperties0(notExtinct,12);
            minGenPred0(notExtinct) = localProperties0(notExtinct,13);
            maxVulPrey0(notExtinct) = localProperties0(notExtinct,14);
            minVulPrey0(notExtinct) = localProperties0(notExtinct,15);
            
        end
        nBasalCon0(patl==1) = nan;
        genAll1 = [genAll1;gen0];
        vulAll1 = [vulAll1;vul0];
        
        
        meanImpPredAll1 = [meanImpPredAll1; meanImpPred0];
        meanImpPreyAll1 = [meanImpPreyAll1; meanImpPrey0];
        
        meanGenPredAll1 = [meanGenPredAll1; meanGenPred0];
        meanVulPreyAll1 = [meanVulPreyAll1; meanVulPrey0];
        
        maxGenPredAll1 = [maxGenPredAll1;maxGenPred0];
        minGenPredAll1 = [minGenPredAll;minGenPred0];
        
        maxVulPreyAll1 = [maxVulPreyAll1;maxVulPrey0];
        minVulPreyAll1 = [minVulPreyAll1;minVulPrey0];
        
        spToBasalAll1 = [spToBasalAll1; spToBasal0];
        
        nBasalConAll = [nBasalConAll;nBasalCon0];
    end
    figure;
    extntGood = extntAll;
    extntGood(:,badWebs) = [];
    histogram(mean(extntGood),'norm','pdf');
    if Z == 100
        extntAll100 = extntAll;
    end
    
end

totalPreyBiomass(:,badWebs) = [];
totalPreyBiomass = reshape(totalPreyBiomass,[],1);

maxVulPreyAll(isnan(maxVulPreyAll)) = 0;
minVulPreyAll(isnan(minVulPreyAll)) = 0;

maxGenPredAll(isnan(maxGenPredAll)) = 0;
minGenPredAll(isnan(minGenPredAll)) = 0;

extntAll100Vec = reshape(extntAll100,[],1);
extntAllVec = reshape(extntAll,[],1);
patlAllNonBasal = patlAll(patlAll>1);
extntAllVecNonBasal = extntAllVec(patlAll>1);

allProperties = full([spToBasalAll,maxSimPredAll,maxSimPreyAll,...
    maxGenPredAll,minGenPredAll,maxVulPreyAll,minVulPreyAll,...
    genAll0,vulAll0,simAll,medianPATLPredDiffAll,medianPATLPreyDiffAll]);

kFree = [1 2];
kPara = [-3 -4];
mAllk1 = 10.^(patlAllNonBasal-1);
mAllk2 = 100.^(patlAllNonBasal-1);
mAll_11 = 10^(kPara(1) - kFree(1))*10.^(kFree(1)*(patlAllNonBasal-1));
mAll_21 = 10^(kPara(1) - kFree(2))*10.^(kFree(2)*(patlAllNonBasal-1));
mAll_12 = 10^(kPara(2) - kFree(1))*10.^(kFree(1)*(patlAllNonBasal-1));
mAll_22 = 10^(kPara(2) - kFree(2))*10.^(kFree(2)*(patlAllNonBasal-1));


xAll_k1 = .314*(mAllk1).^-(0.25);
xAll_k2 = .314*(mAllk2).^-(0.25);
xAll_11 = .314*(mAll_11).^-(0.25);
xAll_21 = .314*(mAll_21).^-(0.25);
xAll_12 = .314*(mAll_12).^-(0.25);
xAll_22 = .314*(mAll_22).^-(0.25);
normalizationChoice = 'cdf';
figure;
hold on
histogram(log10(xAll_k1),'DisplayStyle','stairs'...
    ,'normalization',normalizationChoice)
histogram(log10(xAll_11),'DisplayStyle','stairs'...
    ,'normalization',normalizationChoice)
histogram(log10(xAll_12),'DisplayStyle','stairs'...
    ,'normalization',normalizationChoice)

% histogram(log10(xAll_k2),'DisplayStyle','stairs'...
%     ,'normalization',normalizationChoice)
% histogram(log10(xAll_21),'DisplayStyle','stairs'...
%     ,'normalization',normalizationChoice)
% histogram(log10(xAll_22),'DisplayStyle','stairs'...
%     ,'normalization',normalizationChoice)

legend('All Free Livers, Z=10','All Parasites, Z=10^{-3}','All Parasites, Z = 10^{-4}')
xlabel('log_{10}(x_i)')


figure;
hold on
histogram((xAll_k1),'DisplayStyle','stairs'...
    ,'normalization',normalizationChoice)
histogram((xAll_11),'DisplayStyle','stairs'...
    ,'normalization',normalizationChoice)
histogram((xAll_12),'DisplayStyle','stairs'...
    ,'normalization',normalizationChoice)

% histogram((xAllk2),'DisplayStyle','stairs'...
%     ,'normalization',normalizationChoice)
% histogram((xAll_21),'DisplayStyle','stairs'...
%     ,'normalization',normalizationChoice)
% histogram((xAll22),'DisplayStyle','stairs'...
%     ,'normalization',normalizationChoice)


legend('All Free Livers, Z=10','All Parasites, Z=10^{-3}','All Parasites, Z = 10^{-4}')
xlabel('x_i')

figure;
hold on
kVec = -5:5;
xMxAll = zeros(length(xAll_k1),length(kVec));
normalizationChoice = 'pdf';
for ii = 1:length(kVec)
    
    xMxAll(:,ii) = .314*(10.^(kVec(ii)*(patlAllNonBasal-1))).^(-0.25);
    histogram(log10(xMxAll(:,ii)),'DisplayStyle','stairs'...
        ,'normalization',normalizationChoice);
end
legend(cellstr(num2str(kVec')))
% vulCat = discretize(vulAll,linspace(0,5,41),linspace(1/16,79/16,40));
% meanGenCat = grpstats(genAll,vulCat);
% meanci = grpstats(genAll,vulCat,'meanci');
% err = meanci(:,2) - meanGenCat;
% grpstats(genAll,vulCat,.05)
% %errorbar(linspace(0,max(vulAll),19),meanGenCat,err)
% hold on
%
% vulCat0 = discretize(vulAll0,linspace(0,5,41),linspace(1/16,79/16,40));
% meanGenCat0 = grpstats(genAll0,vulCat0);
% meanci0 = grpstats(genAll0,vulCat0,'meanci');
% err0 = meanci0(:,2) - meanGenCat0;
% %errorbar(linspace(0,max(vulAll0),19),meanGenCat0,err0,'r')
% grpstats(genAll0,vulCat0,.05)
% axis([0 numel(unique(vulCat0)) 0 2])
% hold off

%plot(patlAll,BAll,'k.')

