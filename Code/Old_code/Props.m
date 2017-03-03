%function [fracs] = Web_Properties(res,con,save_glob_save_loc)
% Function that calculates properties of food webs.
%TOP,BAS,INT,CAN,HERB,VUL,GEN
%Initially tried to avoid constructing adjacency matrix, but seems to work
%fine if represented as sparse matrix.

%Program by Nick Kappler, Initial version written 29-10-2015
%I also allow for saving global (save_glob) and local (save_loc) properties
%to file.  I tried to predict somewhat useful local properties for future
%statistical analysis of dynamical results.

%number of links:
L = length(res);

%number of species:
S = max(max([res,con]));

Slist = 1:S;

%Connectance (directed, allowing for self-loops (i.e. canibalism).
C = L/S^2;

%Mathematician's, "column eats row" format:
A = sparse(res,con,true,S,S);
Asym = A+A';

%resource, consumer list without cannibalism
Cann_links = res==con;
nCanns = sum(Cann_links);
Cann = nCanns/S;
Canns = zeros(S,1);
Canns(res(Cann_links))=1;

%Non cannibalistic links
res0 = res(~Cann_links);
con0 = con(~Cann_links);

%Number non-cannibalistic links
L0 = L-nCanns;

%Counting prdators and prey
nPrey = zeros(S,1);
nPred = zeros(S,1);

predation_on_basal = zeros(L0,1);
Bas = false(S,1);
Top = false(S,1);

for ii = Slist
    %counting the number of prey (number of times species ii appears in
    %the consumer list)
    nPrey(ii) = sum(ii == con0);
    
    %counting the number of predators (number of times species ii
    %appears in the resource list)
    nPred(ii) = sum(ii == res0);
    
    if nPrey(ii) == 0
        Bas(ii) = true;
        predation_on_basal(res==ii)=1;
    end
    Top(ii) = nPred(ii)==0;
end

Herbs = zeros(S,1);

for ii = Slist(~Bas) %loop over consuemrs
    %~con0(~predation_on_basal) is the vector of consumers for all 
    %non-cannibalistic, non herbivorous feeding links.
    if sum(ii==con0(~predation_on_basal))==0 %if it only has herbivorous 
        Herbs(ii) = 1;                         %feeding links
    end
end

nHerbs = sum(Herbs);
nBas = sum(Bas);
nTop = sum(Top);

fBas = nBas/S;
fTop = nTop/S;
fHerb = nHerbs/S;

%Including cannibalism in number of predators and prey
nPred = nPred + Canns;
nPrey = nPrey + Canns;

EVul = mean(nPred);
EGen = mean(nPrey);
sdVul = std(nPred,1);
sdGen = std(nPrey,1);

sdLink = std(nPred+nPrey,1);

VulGenCor = corrcoef(nPrey,nPred);

%Similarity Coefficients:
%Predator Similarity, Prey Similarity, link similarity
linkSim = zeros(S,S);
preySim = zeros(S,S);
predSim = zeros(S,S);

for ii = 1:(S-1);
    for jj = (1+ii):S
        
        linkSim(ii,jj) = nnz(Asym(:,ii)&Asym(:,jj))...
            /nnz(Asym(:,ii)|Asym(:,jj));
        
        preySim(ii,jj) = nnz(A(ii,:)&A(jj,:))/nnz(A(ii,:)|A(jj,:));
        
        predSim(ii,jj) = nnz(A(:,ii)&A(:,jj))/nnz(A(:,ii)|A(:,jj));

        
    end 
end

%Loops, chains: numLoops, lengLoops, frac Chains
nLoops = 0;
nChains = 0;
ChainLengs = [];
LoopLengs = [];
inLoop = zeros(S,1);

