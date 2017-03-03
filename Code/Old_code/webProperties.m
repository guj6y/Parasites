%----------------------------------------------------
%% Program by: Perrine Tonin and Coralie Picoche
%% Ref: Dunne et al., PLOS Biology, 2008.
%%      Williams and Martinez, Journal of Animal Ecology, 2008.
%%      Williams and Martinez, Nature, 2000.
%%---------------------------------------------------
%% This function takes a niche web matrix in input
%% (rows eat columns) and calculates the 17 properties
%% of a niche web:
%% Top      Int     Bas     Herb    Can     Omn
%% Loop     ChLen   ChSD    ChNum   TL
%% MaxSim   VulSD   GenSD   LinkSD
%% Path     Clust
%%---------------------------------------------------

function [nicheproperties, herbsp, topsp]= webProperties(nicheweb)

nichewebsize = length(nicheweb);  % number of species
numberlinks  = length(find(nicheweb==1)); % L: number of links in the foodweb


%% Top, Bas, Int and Can coefficients
%%-------------------------------------
    basalsp = find(sum(nicheweb,2)==0); % indices of basal species
    topsp = find(sum(nicheweb,1)==0);   % indices of top species
    cansp = find(diag(nicheweb)==1);    % indices of top species

    Bas=length(basalsp)/nichewebsize;   % percentage of basal species
    Top=length(topsp)/nichewebsize;     % percentage of top species
    Int=1-Top-Bas;         % percentage of intermediate species
    Can=length(cansp)/nichewebsize;     % percentage of cannibal species
%%-------------------------------------


%% Herb coefficient
%% = species eating only at basal level
%%-------------------------------------
    eatingbasalsp=find(sum(nicheweb(:,basalsp),2)>0);    % indices of all species eating basal
    nonbasalsp=find(sum(nicheweb,2)>0);                  % indices of non basal
    nonbasalpreys=nicheweb(eatingbasalsp,nonbasalsp);    % matrix rows=eating basal species / col=non basal species
    herbsp=eatingbasalsp(find(sum(nonbasalpreys,2)==0)); % indicies of herbivore species
    Herb=length(herbsp)/nichewebsize;
%%-------------------------------------


%% MaxSim coefficient
%% mean of the maximum similarity between species
%%-------------------------------------
    Sij=zeros(nichewebsize); % matrix of the trophic similarities
    for i=1:nichewebsize
        for j=1:nichewebsize % or j=1:i-1 ???
            if i~=j
                compred=find(nicheweb(:,i)==1 & nicheweb(:,j)==1); % predators in common
                comprey=find(nicheweb(i,:)==1 & nicheweb(j,:)==1); % preys in common
                totpred=find(nicheweb(:,i)==1 | nicheweb(:,j)==1); % total of predators
                totprey=find(nicheweb(i,:)==1 | nicheweb(j,:)==1); % total of preys
                Sij(i,j)=(length(compred)+length(comprey))/(length(totpred)+length(totprey));
            end
        end
    end
    MaxSim=sum(max(Sij,[],2))/nichewebsize;
%%-------------------------------------


%% VulSD GenSD and LinkSD coefficients
%% variability of the vulnerability (ie number of predators),
%% the generality (ie number of preys) and connectivity (total)
%%-------------------------------------
    Vul=nichewebsize/numberlinks.*sum(nicheweb);
    Gen=nichewebsize/numberlinks.*sum(nicheweb,2)';
    VulSD=std(Vul);
    GenSD=std(Gen);
    LinkSD=std(Vul+Gen);
%%-------------------------------------


%% Clust coefficient
%%-------------------------------------
    Clustvec=zeros(1,nichewebsize);
    nichewithoutcan=nicheweb-eye(nichewebsize); % to prevent counting cannibal species as their own neighbour --> should we do that?
    for v=1:nichewebsize
        neighbours=find(nichewithoutcan(v,:)'==1 | nichewithoutcan(:,v)==1);   % indices of species v neighbours (prey & pred)
        kv=length(neighbours);                                  % number of neighbours of v's species
        betweenneighbours=nicheweb(neighbours,neighbours);      % submatrix of only v's species neighbours
        Clustvec(v)=length(find(betweenneighbours==1))/kv^2;    % Clustering coeff for species v
    end
    Clust=mean(Clustvec);
%%-------------------------------------

%% Clust coefficient
%% mean shortest path between two species
%% Floyd-Warshall algorithm
%%-------------------------------------
niche_nonoriented=max(nicheweb,nicheweb');                  % symetrization of the nicheweb
clustmat=-10^6.*(niche_nonoriented-1)+niche_nonoriented;    % assign high values where there is no link
for k=1:nichewebsize
    for i=1:nichewebsize
        for j=1:nichewebsize
            clustmat(i,j)=min(clustmat(i,j),clustmat(i,k)+clustmat(k,j));
        end
    end
end
clustmat=triu(clustmat); % to keap only pairs of species (clust(i,j)=clust(j,i))
Clust=mean(clustmat(clustmat~=0));
%%-------------------------------------

%% Omn --> fraction of species which consum at least two species and have
%% food chains of different lengths

%% Loop --> fraction of species involved in loops (longer than 1 species loops = cannibalism)

%% ChLen --> Cf Watts & Strogatz 1998
%% ChSD

%% ChNum --> Log of number of food chains (without loops)

%% TL

nicheproperties=[Top Int Bas Herb Can 0 0 0 0 0 0 MaxSim VulSD GenSD LinkSD Path Clust];
    %% Top      Int     Bas     Herb    Can     Omn
%% Loop     ChLen   ChSD    ChNum   TL
%% MaxSim   VulSD   GenSD   LinkSD
%% Path     Clust