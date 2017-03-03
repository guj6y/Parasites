%Calculates P-value for differences in mean between parasites and
%non-freeliving consumers for all local properties. 

%Controls whether or not the webs generated include the observed
%concomittant links.
%concomittantWebs = false;

 %webGeneration;



%Bunch of cell higgery-jiggery here.  This code does t-tests for all the
%properties.  Not grouping webs, so don't need to worry too much about
%independence.
%
%close all
N = 13;
meanPropFreeCon = zeros(6,13);
meanPropPara = zeros(6,13);
varPropFreeCon = zeros(6,13);
varPropPara = zeros(6,13);
sePropPool = zeros(6,13);
dfProp = zeros(6,13);
freeCutOff = -1;
paraCutOff = 1;
%

shortWebName = {'bahia','carp','punta','flens','otago','sylt'};

genAllPara = cell(6,1);
vulAllPara = cell(6,1);

genAllFree = cell(6,1);
vulAllFree = cell(6,1);


for ii = 1:6
    para = propertiesCell{ii,2};
    basal = binAutoCell{ii};
    
    genAllPara{ii} = propertiesCell{ii,1}(para,2);
    vulAllPara{ii} = propertiesCell{ii,1}(para,3);
    
    genAllFree{ii} = propertiesCell{ii,1}(~(para|basal),2);
    vulAllFree{ii} = propertiesCell{ii,1}(~(para|basal),3);
    
    meanPropFreeCon(ii,:) = mean(propertiesCell{ii,1}(...
        ~(para|basal),:));
    meanPropPara(ii,:) = mean(propertiesCell{ii,1}(para,:));
    varPropFreeCon(ii,:) = var(propertiesCell{ii,1}(...
        ~(para|basal),:));
    varPropPara(ii,:) = var(propertiesCell{ii,1}(para,:));
    sePropPool(ii,:) = sqrt(varPropFreeCon(ii,:)/totFreeCon(ii) + ...
        varPropPara(ii,:)/totPara(ii));
    dfProp(ii,:) = (sePropPool(ii,:).^4)./...
    (((varPropFreeCon(ii,:)./totFreeCon(ii)).^2./(totFreeCon(ii)-1))+...
    ((varPropPara(ii,:)./totPara(ii)).^2./(totPara(ii)-1)));
    

    gen_vulHist = figure('Position',[0 0 1440 900]);
    sp1 = subplot(2,2,1);
    minGen = min([genAllFree{ii};genAllPara{ii}]);
    maxGen = max([genAllFree{ii};genAllPara{ii}]);
    histogram(genAllFree{ii},10,...
        'BinLimits',[minGen maxGen]...
        ...,'Normalization','pdf'...
        )
    avg = mean(genAllFree{ii});
    line([avg avg],get(sp1,'YLim'),'color','r')
    title(sprintf('Generality of Free-Liviers, %s',shortWebName{ii}))
    xlabel(sprintf(...
        'Normalized Generality, mean = %.2f',avg))
    
    sp2 = subplot(2,2,2);
    minVul = min([vulAllFree{ii};vulAllPara{ii}]);
    maxVul = max([vulAllFree{ii};vulAllPara{ii}]);
    histogram(vulAllFree{ii}...
        ,'BinLimits',[minVul maxVul]...
        ...,'Normalization','pdf'...
        )
    avg = mean(vulAllFree{ii});
    line([avg avg],get(sp2,'YLim'),'color','r')
    title(sprintf('Vulnerability of Free-Liviers, %s',shortWebName{ii}))
    xlabel(sprintf(...
        'Normalized Vulnerality, mean = %.2f',avg))
    
    sp3 = subplot(2,2,3);
    histogram(genAllPara{ii}...
        ,'BinLimits',[minGen maxGen]...
        ...,'Normalization','pdf'...
        )
    avg = mean(genAllPara{ii});
    line([avg avg],get(sp3,'YLim'),'color','r')
    title(sprintf('Generality of Parasites, %s',shortWebName{ii}))
    xlabel(sprintf(...
        'Normalized Generality, mean = %.2f',avg))
    sp4 = subplot(2,2,4);
    histogram(vulAllPara{ii}...
        ,'BinLimits',[minVul maxVul]...
        ...,'Normalization','pdf'...
        )
    title(sprintf('Vulnerability of Parasites, %s',shortWebName{ii}))
     avg = mean(vulAllPara{ii});
    line([avg avg],get(sp4,'YLim'),'color','r')
    xlabel(sprintf(...
        'Normalized Vulnerality, mean = %.2f',avg))
    
    name = sprintf('gen-vul_para-freeHists_%s',shortWebName{ii});
    
    saveas(gen_vulHist,name,'jpg')
    close(gen_vulHist);
end

tProp = (meanPropFreeCon - meanPropPara)./sePropPool;
tProp = tProp(:,1:11);
dfProp = dfProp(:,1:11);
sePropPool = sePropPool(:,1:11);


 
pProp = 2*tcdf(-abs(tProp),dfProp);

varNames = {'clustCoef','gen','vul','$<$vulPrey$>$','$<$impPrey$>$','$<$genPred$>$','$<$ImpPred$>$','minTL','basalCon','SWTL','inLoop'};
fileName = 'pValsAllProps.textab';
fid = fopen(fileName,'w');
fprintf(fid,'webName&clustCoef&gen&vul&meanVulPrey&meanImpPrey&meanGenPred&meanImpPred&minSPToBasal&numConnBasal&SWTL&inLoop\\\\\n');
dlmwrite(fileName,pProp.*sign(tProp),'delimiter','&','-append')
fclose(fid);
shortWebName = shortWebName';
pVals = pProp.*sign(tProp);


pValsCell = cell(size(pVals) + [0 1]);
pValsCell(:,1) = deal(shortWebName);
pValsCell(:,2:end) = deal(num2cell(pVals));
fileName = '../../Presentations/Presentation-4-11-16/pValsAllProps2.textab';
pValsCell = pValsCell';

sigs = zeros(size(pValsCell));
sigs(2:end,:) = cellfun(@(x) lt(abs(x),0.05),pValsCell(2:end,:));
insigs = ~ sigs;
insigs(1,:) = false;
pValsCell(sigs>0) = deal(cellfun(@(x) sprintf('\\textcolor{red}{%.2g}',x),pValsCell(sigs>0),'UniformOutput',false));
pValsCell(insigs>0) = deal(cellfun(@(x) sprintf('%.2g',x),pValsCell(insigs>0),'UniformOutput',false));
fid = fopen(fileName,'w');
fprintf(fid,'web&clustCoef&gen&vul&meanVulPrey&meanImpPrey&\\\\\n\\hline\n');
fprintf(fid,'%s&%s&%s&%s&%s&%s&\\\\\n',pValsCell{1:6,:});
fprintf(fid,'\\hline\n\\hline\n web&meanGenPred&meanImpPred&minSPToBasal&numConnBasal&SWTL&inLoop\\\\\n\\hline\n');
fprintf(fid,'%s&%s&%s&%s&%s&%s&%s\\\\\n',pValsCell{[1,7:end],:});
fclose(fid);

%Looking at correlations in the data
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

%Take non0basal species
corr_sp = classAll>0;

corrMapPears = figure;
set(gcf, 'Visible', 'off')
propCorrPear = corr([propertiesAll(corr_sp,1:11),paraAll(corr_sp)]);
imagesc(propCorrPear,[-1,1]);
xlabels = {'Clust',...
           'Gen',...
           'Vul',...
           '<Vul Prey>',...
           '<ImpPrey>',...
           '<GenPred>',...
           '<ImpPred>',...
           'minTL',...
           'basalCon',...
           'SWTL',...
           'inLoop',...
           'para'};
hold on
plot([11.5,11.5],[0,13],'k-','LineWidth',2)
plot([0,13],[11.5,11.5],'k-','LineWidth',2)
textStrings = num2str(propCorrPear(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:12);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
colormap(jet)
colorbar
title('Heat Map of Pearson''s (linear) Correlation Matrix')
ax = gca;
ax.XTickLabel = xlabels;
ax.XTick = 1:12;
ax.YTick = 1:12;
ax.XTickLabelRotation=45;
ax.YTickLabel = xlabels;
print(corrMapPears,'-djpeg','../../Presentations/Presentation-4-11-16/corMapPear.jpg');
close(corrMapPears);

corrMapSpear = figure;
set(gcf,'Visible','off')
propCorrSpear = corr([propertiesAll(corr_sp,1:11),paraAll(corr_sp)],...
    'type','spearman');

imagesc(propCorrSpear,[-1,1]);
hold on
plot([11.5,11.5],[0,13],'k-','LineWidth',2)
plot([0,13],[11.5,11.5],'k-','LineWidth',2)
textStrings = num2str(propCorrSpear(:),'%0.2f');  %# Create strings from the matrix values
textStrings = strtrim(cellstr(textStrings));  %# Remove any space padding
[x,y] = meshgrid(1:12);   %# Create x and y coordinates for the strings
hStrings = text(x(:),y(:),textStrings(:),...      %# Plot the strings
                'HorizontalAlignment','center');
colormap(jet)
colorbar
title('Heat Map of Spearman''s (rank) Correlation Matrix')
ax = gca;
ax.XTickLabel = xlabels;
ax.XTick = 1:12;
ax.YTick = 1:12;
ax.XTickLabelRotation=45;
ax.YTickLabel = xlabels;
print(corrMapSpear,'-djpeg','../../Presentations/Presentation-4-11-16/corMapSpear.jpg');
close(corrMapSpear);
hold off

propKey = cell(12,2);
propNums = 1:12;
propNumsCell = num2cell(propNums)';
propKey(:,1) = propNumsCell;
propKey(1:11,2) = varNames;
propKey(12,2) = {'para'};
propKey = propKey';
fileName = '../../Presentations/Presentation-4-11-16/variableKey.textab';
fid = fopen(fileName,'w');
fprintf(fid,'Number&Name\\\\\n\\hline\n');
fprintf(fid,'%u&%s\\\\\n',propKey{:});
fclose(fid);
%}

% close all



%Examining body SIze ratios & body size distributions
bodySizeRatios = cell(6,1);
bodySizesRes = cell(6,1);
bodySizesCon = cell(6,1);

endoByCon{ii} = cell(6,1);
ectoInvertByCon{ii} = cell(6,1);
ectoVertByCon{ii} = cell(6,1);
paraByCon{ii} = cell(6,1);
paraByRes{ii} = cell(6,1);
paraLink = cell(6,1);
bodySizes = cell(6,1);
swtls = cell(6,1);
bodySizesEctoInvert{ii} = cell(6,1);
bodySizesEctoVert = cell(6,1);
bodySizesEndo = cell(6,1);
bodySizesPara = cell(6,1);

swtlEctoInvert{ii} = cell(6,1);
swtlEctoVert = cell(6,1);
swtlEndo = cell(6,1);
swtlPara = cell(6,1);

binEctoInvertByCon = cell(6,1);
binEndoByCon = cell(6,1);
binEctoVertByCon = cell(6,1);
binParaLink = cell(6,1);

bodySizeSWTL = figure;
bodyRatioSWTL = figure;

dotColors0 = ['b';'k';'g';'r'];
dotStyles0 = ['o';'x';'+';'*';'s';'d'];
dotStyles = reshape(repmat(dotStyles0,1,4)',24,1);
dotColors = repmat(dotColors0,6,1);

dotChars = [dotStyles dotColors];

bodySizesEV = [];
bodySizesEI = [];
bodySizesEN = [];
bodySizesPA = [];

swtlsEV = [];
swtlsEI = [];
swtlsEN = [];
swtlsPA = [];

bodySizeRatiosEV = [];
bodySizeRatiosEI = [];
bodySizeRatiosEN = [];
bodySizeRatiosPA = [];

for ii = 1:6
    
    nEctoInvert = sum(binEctoInvertCell{ii});
    nEctoVert = sum(binEctoVertCell{ii});
    nEndo = sum(binEndoCell{ii});
    nPara = sum(binParaCell{ii});
    
    linkList = linkListCell{ii};
    res = linkList(:,1);
    con = linkList(:,2);
    
    bodySizes{ii} = propertiesCell{ii,1}(:,12);
    swtls{ii} = propertiesCell{ii,1}(:,10);
    bodySizesRes{ii} = bodySizes{ii}(res);
    bodySizesCon{ii} = bodySizes{ii}(con);
    
    bodySizesEctoInvert{ii} = bodySizes{ii}(binEctoInvertCell{ii});
    bodySizesEctoVert{ii} = bodySizes{ii}(binEctoVertCell{ii});
    bodySizesEndo{ii} = bodySizes{ii}(binEndoCell{ii});
    bodySizesPara{ii} = bodySizes{ii}(binParaCell{ii});
    
    swtlEctoInvert{ii} = swtls{ii}(binEctoInvertCell{ii});
    swtlEctoVert{ii} = swtls{ii}(binEctoVertCell{ii});
    swtlEndo{ii} = swtls{ii}(binEndoCell{ii});
    swtlPara{ii} = swtls{ii}(binParaCell{ii});
    
    swtlEctoInvert{ii} = swtlEctoInvert{ii}(isfinite(bodySizesEctoInvert{ii}));
    swtlEctoVert{ii} = swtlEctoVert{ii}(isfinite(bodySizesEctoVert{ii}));
    swtlEndo{ii} = swtlEndo{ii}(isfinite(bodySizesEndo{ii}));
    swtlPara{ii} = swtlPara{ii}(isfinite(bodySizesPara{ii}));
    
    bodySizesEctoInvert{ii} = bodySizesEctoInvert{ii}(isfinite(bodySizesEctoInvert{ii}));
    bodySizesEctoVert{ii} = bodySizesEctoVert{ii}(isfinite(bodySizesEctoVert{ii}));
    bodySizesEndo{ii} = bodySizesEndo{ii}(isfinite(bodySizesEndo{ii}));
    bodySizesPara{ii} = bodySizesPara{ii}(isfinite(bodySizesPara{ii})); 
    
    bodySizeRatios{ii} = bodySizesCon{ii}./bodySizesRes{ii};

    paraByCon{ii} = propertiesCell{ii,2}(con);
    paraByRes{ii} = propertiesCell{ii,2}(res);
    
    binEndoByCon{ii} = binEndoCell{ii}(con);
    binEctoInvertByCon{ii} = binEctoInvertCell{ii}(con);
    binEctoVertByCon{ii} = binEctoVertCell{ii}(con);
    binParaLink{ii} = paraByCon{ii}&(~paraByRes{ii});
     
    endoByCon{ii} = con(binEndoCell{ii}(con));
    ectoInvertByCon{ii} = con(binEctoInvertByCon{ii});
    ectoVertByCon{ii} = con(binEctoVertByCon{ii});
    paraLink{ii} = con(binParaLink{ii});
    
    figure(bodySizeSWTL)
    hold on
    plot(swtlEctoInvert{ii},log10(bodySizesEctoInvert{ii}),dotChars(4*(ii-1)+1,:))
    plot(swtlEctoVert{ii},log10(bodySizesEctoVert{ii}),dotChars(4*(ii-1)+2,:))
    plot(swtlEndo{ii},log10(bodySizesEndo{ii}),dotChars(4*(ii-1)+3,:))
    plot(swtlPara{ii},log10(bodySizesPara{ii}),dotChars(4*(ii-1)+4,:))
    
    figure(bodyRatioSWTL)
    hold on
    plot(swtls{ii}(ectoInvertByCon{ii}),log10(bodySizeRatios{ii}(binEctoInvertCell{ii}(con))),dotChars(4*(ii-1)+1,:))
    plot(swtls{ii}(ectoVertByCon{ii}),log10(bodySizeRatios{ii}(binEctoVertCell{ii}(con))),dotChars(4*(ii-1)+3,:))
    plot(swtls{ii}(endoByCon{ii}),log10(bodySizeRatios{ii}(binEndoCell{ii}(con))),dotChars(4*(ii-1)+4,:))
    plot(swtls{ii}(paraLink{ii}),log10(bodySizeRatios{ii}(paraLink{ii})),dotChars(4*(ii-1)+2,:))
    
    bodySizesEV = [bodySizesEV; bodySizesEctoInvert{ii}];
    bodySizesEI = [bodySizesEI; bodySizesEctoVert{ii}];
    bodySizesEN = [bodySizesEN; bodySizesEndo{ii}];
    bodySizesPA = [bodySizesPA; bodySizesPara{ii}];
    
    swtlsEV = [swtlsEV; swtlEctoInvert{ii}];
    swtlsEI = [swtlsEI; swtlEctoVert{ii}];
    swtlsEN = [swtlsEN; swtlEndo{ii}];
    swtlsPA = [swtlsPA; swtlPara{ii}];
    
    if ii < 4
    bodySizeRatiosEI = [bodySizeRatiosEI; bodySizeRatios{ii}(binEctoInvertCell{ii}(con))];
    bodySizeRatiosEV = [bodySizeRatiosEV; bodySizeRatios{ii}(binEctoVertCell{ii}(con))];
    bodySizeRatiosEN = [bodySizeRatiosEN; bodySizeRatios{ii}(binEndoCell{ii}(con))];
    bodySizeRatiosPA = [bodySizeRatiosPA; bodySizeRatios{ii}(paraByCon{ii}&(~paraByRes{ii}))];
    end
end

%Ignorning mosquitos here.  This study is not so much about parasites as it
%is about the effect of body size ratios.
bodySizeRatiosFree = [bodySizeRatiosEI(log10(bodySizeRatiosEI)>freeCutOff);
                      bodySizeRatiosEV;
                      bodySizeRatiosEN;
                      bodySizeRatiosPA(log10(bodySizeRatiosPA)>paraCutOff)];

%All the small things that eat big things.
bodySizeRatiosPara_Mosq = [bodySizeRatiosPA(log10(bodySizeRatiosPA)<paraCutOff);
                         bodySizeRatiosEI(log10(bodySizeRatiosEI)<freeCutOff)];

propertiesAll
                     
xEV = [swtlsEV ones(size(swtlsEV))];
xEI = [swtlsEI ones(size(swtlsEI))];
xEN = [swtlsEN ones(size(swtlsEN))];
xPA = [swtlsPA ones(size(swtlsPA))];

bEV = xEV\log10(bodySizesEV);
bEI = xEI\log10(bodySizesEI);
bEN = xEN\log10(bodySizesEN);
bPA = xPA\log10(bodySizesPA);

rhoEV = corr(log10(bodySizesEV),xEV*bEV);
rhoEN = corr(log10(bodySizesEN),xEN*bEN);
rhoEI = corr(log10(bodySizesEI),xEI*bEI);
rhoPA = corr(log10(bodySizesPA),xPA*bPA);

rSquaredEV = rhoEV^2;
rSquaredEN = rhoEN^2;
rSquaredEI = rhoEI^2;
rSquaredPA = rhoPA^2;

figure(bodySizeSWTL)
hold on 
plot(swtlsEI,xEI*bEI,'k-')
plot(swtlsEV,xEV*bEV,'b-')
[swtlsEN,idx] = sort(swtlsEN);
plot(swtlsEN,xEN(idx,:)*bEN,'g-')
plot(swtlsPA,xPA*bPA,'r-')
xlabel('Short-Weighted Trophic Level')
ylabel('log_{10}[Body Size]')
title('Body Size vs. SWTL')
legend('Ect. Inv.','Ect. Vert.','Endo.','para','Location','southwest')
axis([2 6 -7 4])
print(bodySizeSWTL,'-djpeg','../../Presentations/Presentation-4-11-16/bstl.jpg')

figure(bodyRatioSWTL)
xlabel('Short-Weighted Trophic Level')
ylabel('log_{10}[consumer/resource bodySizes]')
title('Body Size vs. SWTL (by link)')
legend('Ect. Inv.','para','Ect. Vert.','Endo.','Location','southwest')
axis([2 6 -10 10])
print(bodyRatioSWTL,'-djpeg','../../Presentations/Presentation-4-11-16/bsrtl.jpg')

BMR = figure;
%set(gcf,'visible','off')
title('Body Mass Ratios by Link')
subplot(2,2,1)
hist(log10(bodySizeRatiosEI))
title('Ectotherm Invertebrate Consumers')
ylabel('Frequency')
subplot(2,2,2)
hist(log10(bodySizeRatiosEV))
title('Ectotherm Vertebrate Consumers')
subplot(2,2,3)
hist(log10(bodySizeRatiosEN))
title('Endotherm Consumers')
xlabel('log_{10} Consumer Resource Body mass Ratio')
ylabel('Frequency')
subplot(2,2,4)
hist(log10(bodySizeRatiosPA))
title('Parasitic Consumers (no hyper-parasitism)')
xlabel('log_{10} Consumer Resource Body mass Ratio')
%Problem? parasites plot includes links froma free-living stages
print(BMR,'-djpeg','../../Presentations/Presentation-4-11-16/BMR.jpg');
%close(BMR);

BMRAll = figure;
hist(log10([bodySizeRatiosEI;bodySizeRatiosEV;bodySizeRatiosEN;bodySizeRatiosPA]))
title('BodyMassRatios for All!')
xlabel('log_{10}(BMR)')
%{
medEI0 = median(log10(bodySizeRatiosEI),'omitnan');
medEI = median(log10(bodySizeRatiosEI(log10(bodySizeRatiosEI)>-3)),'omitnan');
medEV = median(log10(bodySizeRatiosEV),'omitnan');
medEN = median(log10(bodySizeRatiosEN),'omitnan');
medPA = median(log10(bodySizeRatiosPA),'omitnan');

meanEI = mean(log10(bodySizeRatiosEI(log10(bodySizeRatiosEI)>-3)),'omitnan');
meanEI0 = median(log10(bodySizeRatiosEI),'omitnan');
meanEV = mean(log10(bodySizeRatiosEV),'omitnan');
meanEN = mean(log10(bodySizeRatiosEN),'omitnan');
meanPA = mean(log10(bodySizeRatiosPA),'omitnan');

iqrEI = iqr(log10(bodySizeRatiosEI(log10(bodySizeRatiosEI)>-3)));
iqrEV = iqr(log10(bodySizeRatiosEV));
iqrEN = iqr(log10(bodySizeRatiosEN));
iqrPA = iqr(log10(bodySizeRatiosPA));

stdEI = std(log10(bodySizeRatiosEI(log10(bodySizeRatiosEI)>-3)),'omitnan');
stdEV = std(log10(bodySizeRatiosEV),'omitnan');
stdEN = std(log10(bodySizeRatiosEN),'omitnan');
stdPA = std(log10(bodySizeRatiosPA),'omitnan');
%}
medEI0 = median(log(bodySizeRatiosEI),'omitnan');
medEI = median(log(bodySizeRatiosEI(log10(bodySizeRatiosEI)>freeCutOff)),'omitnan');
medEV = median(log(bodySizeRatiosEV),'omitnan');
medEN = median(log(bodySizeRatiosEN),'omitnan');
medPA = median(log(bodySizeRatiosPA),'omitnan');

medFree = median(log(bodySizeRatiosFree),'omitnan');
medPA_Mosq = median(log(bodySizeRatiosPara_Mosq),'omitnan');

meanEI = mean(log(bodySizeRatiosEI(log10(bodySizeRatiosEI)>freeCutOff)),'omitnan');
meanEI0 = median(log(bodySizeRatiosEI),'omitnan');
meanEV = mean(log(bodySizeRatiosEV),'omitnan');
meanEN = mean(log(bodySizeRatiosEN),'omitnan');
meanPA = mean(log(bodySizeRatiosPA),'omitnan');

meanFree = mean(log(bodySizeRatiosFree),'omitnan');
meanPA_Mosq = mean(log(bodySizeRatiosPara_Mosq),'omitnan');

iqrEI = iqr(log(bodySizeRatiosEI(log10(bodySizeRatiosEI)>freeCutOff)));
iqrEV = iqr(log(bodySizeRatiosEV));
iqrEN = iqr(log(bodySizeRatiosEN));
iqrPA = iqr(log(bodySizeRatiosPA));

iqrFree = iqr(log(bodySizeRatiosFree));
iqrPA_Mosq = iqr(log(bodySizeRatiosPara_Mosq));

stdEI = std(log(bodySizeRatiosEI(log10(bodySizeRatiosEI)>freeCutOff)),'omitnan');
stdEV = std(log(bodySizeRatiosEV),'omitnan');
stdEN = std(log(bodySizeRatiosEN),'omitnan');
stdPA = std(log(bodySizeRatiosPA),'omitnan');

stdFree = std(log(bodySizeRatiosFree),'omitnan');
stdPA_Mosq = std(log(bodySizeRatiosPara_Mosq),'omitnan');
%}
%newNamesCell{1}(linkListCell{1}((log10(bodySizeRatios{1})<-3)&binEctoInvertCell{1}(linkListCell{1}(:,2),1),:))
