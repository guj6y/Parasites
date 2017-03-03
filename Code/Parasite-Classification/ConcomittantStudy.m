%This script attempts to find patterns with concomittant links
%
% concomittant links are this:
% p: Parasite
% h: host of parasite
% c: consumer of host
%
%         _o°°°°°°°°°o¸_
%        (P)?¬     ???(C)
%         ?   \   /    ?
%              \_/
%              (H)
%               ¯
%(not ovaries)
%The parasite (p) has a host (H) that is eaten by a consumer (c).  in some
%situations, by eating H, C also eats some P. (the dotted line).  This
%doesn't always happen, so this code tries to figure out roughly when it
%does.



hosts = cell(6,1);
%webGeneration
names = {'BSQ','CSM','EPB'};
stats = cell(6,1);
trees = cell(6,1);
data = cell(6,1);
paraList = cell(6,1);

concomittant = cell(6,1);
predation = cell(6,1);
linkListConcomittant = cell(6,1);

for ii = 1:6
    para = propertiesCell{ii,2};
    S = totSpecies(ii);
    spList = 1:S;
    paraList{ii} = spList(para);
    
    linkListConcomittant{ii} = [];
    res = linkListCell{ii}(:,1);
    con = linkListCell{ii}(:,2);
    linkTypes = linkTypesCell{ii};
    concLinks = (linkTypes==14) | (linkTypes==15) ;
    hosts{ii} = unique(res(linkParaCell{ii}>0));
    fracConcs1 = []; 
    fracConcs2 = [];
    microLinks = linkTypes == 3;
    bsh1 = [];
    bsc = [];
    swTLh1 = [];
    swTLc = [];
    
    
    bsh2 = [];
    bsp = [];
    swTLh2 = [];
    swTLp = [];
    
    swTLh3 = [];
    swTLp2 = [];
    swTLc2 = [];
    
    
    bsh3 = [];
    bsp2 = [];
    bsc2 = [];
    
    
    swTLs = propertiesCell{ii,1}(:,10);
    bs = propertiesCell{ii,1}(:,12);
    
    
    %"For each host"
    for jj = hosts{ii}'
        if sum(jj==paraList{ii})
            continue
        else
            predOfjj = con((res==jj)&(~linkParaCell{ii}));
            paraOfjj = con((res==jj)&(linkParaCell{ii}));
            %THe following are looking for patterns int he fractions of
            %parasites of a host that each predator of that host consumes:
            
            %For each parasite-host pair, what are the odds that a parasite
            %of that host is concomittantly eaten by a predator of that
            %host.  But if that parasite also eats a different host that is
            %a prey of that same predator, we wouldn't know whether one or
            %both were the source of the existing concomittant link.
            
            %Additionaly, these data should probably be taken with grain of
            %salt: the evidence for all of these links is 'inferred'.
            
            %This is looking by predators of each host (percent of parasites of
            %each host that a predator of that host consumes
            %         for kk = predOfjj'
            %             swTLh1 = [swTLh1;swTLs(jj)];
            %             swTLc = [swTLc;swTLs(kk)];
            %             bsh1 = [bsh1;bs(jj)];
            %             bsc = [bsc;bs(kk)];
            %             actualConcLinks = 0;
            %             for ll = paraOfjj'
            %                 actualConcLinks = actualConcLinks+ ...
            %                     sum((con==kk)&concLinks&(res==ll));
            %             end
            %             possibleConcLinks = length(paraOfjj);
            %             fracConcs1 = [fracConcs1; actualConcLinks/possibleConcLinks];
            %         end
            %         %This is looking by parasites of each host (percent of predators of
            %         %each host that concomittantly consume a particular parasite)
            %         for kk = paraOfjj'
            %             swTLh2 = [swTLh2;swTLs(jj)];
            %             swTLp = [swTLp;swTLs(kk)];
            %             bsh2 = [bsh2;bs(jj)];
            %             bsp = [bsp;bs(kk)];
            %             actualConcLinks = 0;
            %             for ll = predOfjj'
            %                 actualConcLinks = actualConcLinks + ...
            %                     sum((con==ll)&concLinks&(res==kk));
            %             end
            %             possibleConcLinks = length(predOfjj);
            %
            %             fracConcs2 = [fracConcs2; actualConcLinks/possibleConcLinks];
            %         end
            
            %Fully separating everything
            
            for kk = paraOfjj'
                for ll = predOfjj'
                    swTLh3 = [swTLh3; swTLs(jj)];
                    swTLp2 = [swTLp2; swTLs(kk)];
                    swTLc2 = [swTLc2; swTLs(ll)];
                    bsh3 = [bsh3; bs(jj)];
                    bsp2 = [bsp2; bs(kk)];
                    bsc2 = [bsc2; bs(ll)];
                    concomittant{ii} = [concomittant{ii}; ...
                        sum((con==ll)&(res==kk)&concLinks)];
                    predation{ii} = [predation{ii} ...
                        sum((con==ll)&(res==kk)&~concLinks)];
                    linkListConcomittant{ii} = [linkListConcomittant{ii};
                        kk ll jj];
                    %NB:  Here I am keeping track of the host, so
                    %linkListConcomittant is going to wind up being much
                    %larger than it needs to be.  It may be interesting
                    %and/or important to actually track through which hosts
                    %the parasite is incidentally devoured.  Especially if
                    %(as I suspect may be the most  sound method of
                    %modeling parasites in general) anyone down the line
                    %winds up dividing the parasite species into
                    %populations within each of its hosts(along with a
                    %free-living subset).  In this way the parasites might
                    %be able to be thought of as a meta-population.  that
                    %scneario opens up a lot of very interesting options
                    %for theoretical work, especially viz. Chesson's
                    %storage effect.  The parasites are effectively
                    %changing who their predators are by inhabiting
                    %different hosts( or being free-living).  
                    
                    %For now(probably as far as NK will take it) I am not
                    %separating the different incidental links.
                end
            end
            
            
        end
    end
    sum(concomittant{ii}&predation{ii}')
        data{ii} = [swTLh3 swTLp2 swTLc2 log(bsh3) log(bsp2) log(bsc2)];
        data{ii} = data{ii}(~predation{ii},:);
        concomittant{ii} = concomittant{ii}(~predation{ii});
        linkListConcomittant{ii} = linkListConcomittant{ii}(~predation{ii},:);
        
        
        %Compactifying the data
        linkList = [0 0];
        for jj= length(linkListConcomittant{ii}):-1:1
            p = linkListConcomittant{ii}(jj,1);
            c = linkListConcomittant{ii}(jj,2);
            if (sum((linkList(:,1)==p)&(linkList(:,2)==c)))
                linkListConcomittant{ii}(jj,:) = [];
                data{ii}(jj,:) = [];
                concomittant{ii}(jj,:) = [];
            else
                linkList = [linkList;
                    p, c];
            end
            
        end
        data{ii}(:,[1, 4]) = [];
        %    [~,~,stats{ii}] = mnrfit([log(bsp2),log(bsh3),log(bsc2)],categorical(concomittant));
        %     figure
        %     hold on
        %     title(sprintf('%s',names{ii}))
        %     plot(stats{ii}.beta(1) + stats{ii}.beta(2)*log(bsp2)+stats{ii}.beta(3)*log(bsh3)+stats{ii}.beta(4)*log(bsc2),concomittant,'k.')
        %     t = linspace(-5,10,100);
        %     plot(t,1./(1+exp(t)))
        %     legend('data',sprintf('\\beta_1 = %.3f \\pm %.3f\n \\beta_2 = %.3f \\pm %.3f',stats{ii}.beta(1),1.96*stats{ii}.se(1),stats{ii}.beta(2),1.96*stats{ii}.se(2)))
        
        

end

variableNames = {'TL p','TL c','BS p','BS c'};
%Lumping all the data
xData = [data{1};
    data{2};
    data{3}];

yData = [concomittant{1};
    concomittant{2};
    concomittant{3}];

nans = isnan(sum(xData,2));

xData = xData(~nans,:);
yData = yData(~nans);

[~,order1] = sort(xData(:,1));
[~,order2] = sort(xData(:,2));
[~,order3] = sort(xData(:,3));
[~,order4] = sort(xData(:,4));

xDataO = [order1 order2 order3 order4]/length(order1);

cvPart = cvpartition (yData,'holdout',0.3);

xTrain = xData(training(cvPart),:);
yTrain = yData(training(cvPart),:);
xTest = xData(test(cvPart),:);
yTest = yData(test(cvPart),:);

xTrainO = xDataO(training(cvPart),:);
xTestO = xDataO(test(cvPart),:);

%Using all predictors
tree = fitctree(xTrain,yTrain,'PredictorNames',variableNames);
treeO = fitctree(xTrainO,yTrain,'PredictorNames',variableNames);

% % THe best pruning level?
% [E,SE,NLEAF,BESTLEVEL] = cvloss(tree,'Subtrees','All');
%
% cvtree = crossval(tree); %%%This takes some time to run.
% cvloss = kfoldLoss(cvtree);
%
% %Credit to the mathworks online documentation for this code; it takes a hot
% %second to run because of the cross validation.
% leafs = logspace(1,2,10);
% N = numel(leafs);
% err = zeros(N,1);
%
% for n=1:N
%     t = fitctree(xData,yData,'CrossVal','On',...
%         'KFold',3,...
%         'MinLeafSize',leafs(n));
%     err(n) = kfoldLoss(t);
% end
% figure
% plot(leafs,err);
% xlabel('Min Leaf Size');
% ylabel('cross-validated error');
%
%
% splits = linspace(1,50,50);
% N = numel(splits);
% err = zeros(N,1);
%
% for n=1:N
%     t = fitctree(xData,yData,'CrossVal','On',...
%         'KFold',3,...
%         'MaxNumSplits',splits(n));
%     err(n) = kfoldLoss(t);
% end
% figure
% plot(splits,err);
% xlabel('Max Num Splits ');
% ylabel('cross-validated error');
%
% %Trying some ensemble methods.  NOt sure why RUSboost acts so strangely.
% %In the end I think these just aren't worth it for this problem.
%
% adaBoostEnsemble = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,'tree');
% adaBoostCVEnsemble = fitensemble(xData,yData,'AdaBoostM1',1000,'tree',...
%     'type','classification','kfold',5);
% figure;
% plot(loss(adaBoostEnsemble,xTest,yTest,'Mode','cumulative'));
% hold on;
% plot(kfoldLoss(adaBoostCVEnsemble,'mode','cumulative'),'r.');
% hold off;
% xlabel('Number of trees');
% ylabel('Classificaiton error');
% legend('Test','Cross-validation','Location','NE')
% title('Adaboost algorithm')
%
%
% temp = templateTree('MaxNumSplits',2);
% adaBoostEnsemble2 = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,temp);
% temp = templateTree('MaxNumSplits',3);
% adaBoostEnsemble3 = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,temp);
% temp = templateTree('MaxNumSplits',4);
% adaBoostEnsemble4 = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,temp);
% temp = templateTree('MaxNumSplits',5);
% adaBoostEnsemble5 = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,temp);
% temp = templateTree('MaxNumSplits',10);
% adaBoostEnsemble10 = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,temp);
% temp = templateTree('MaxNumSplits',20);
% adaBoostEnsemble20 = fitensemble(xTrain,yTrain,'AdaBoostM1',1000,temp);
%
%
%
% figure;
% hold on;
% plot(loss(adaBoostEnsemble,xTest,yTest,'mode','cumulative'),'r.');
% plot(loss(adaBoostEnsemble2,xTest,yTest,'mode','cumulative'),'b.');
% plot(loss(adaBoostEnsemble3,xTest,yTest,'mode','cumulative'),'g.');
% plot(loss(adaBoostEnsemble4,xTest,yTest,'mode','cumulative'),'c.');
% plot(loss(adaBoostEnsemble5,xTest,yTest,'mode','cumulative'),'m.');
% plot(loss(adaBoostEnsemble10,xTest,yTest,'mode','cumulative'),'k.');
% plot(loss(adaBoostEnsemble20,xTest,yTest,'mode','cumulative'),'ro');
% hold off;
% xlabel('Number of trees');
% ylabel('Classificaiton error');
% legend('Stump','2','3','4','5','10','20')
% title('Adaboost algorithm')
%
%
% % rusBoostEnsemble = fitensemble(xTrain,yTrain,'RUSBoost',1000,'tree');
% % rusBoostCVEnsemble = fitensemble(xData,yData,'RUSBoostM1',1000,'Tree',...
% %     'type','classification','kfold',5);
% % figure;
% % plot(loss(rusBoostEnsemble,xTest,yTest,'Mode','cumulative'));
% % hold on;
% % plot(kfoldLoss(rusBoostCVEnsemble,'mode','cumulatie'),'r.');
% % hold off;
% % xlabel('Number of trees');
% % ylabel('Classificaiton error');
% % figure;
% % plot(loss(rusBoostEnsemble,xTest,yTest,'mode','cumulative'));
% % xlabel('Number of trees');
% %  ylabel('Classificaiton error');
% %  grid on;
% % legend('Test','Cross-validation','Location','NE')
% % title('RUSboost algorithm')