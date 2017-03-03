function [res,con,n,c,r,noobs] = AddTrophicConsumers(n,c,r,N,C0)

%Add trophic consumers to a niche web.  Attempt to maintain the current
%connectance.  This also needs the original connectance you were aiming for
%when you generated the web since that tells us how the existing ranges are
%distributed.

%Defining the expected width of parasitic ranges.  Adjusting for the fact
%that each parasite is guaranteed to have a host.

%Went back and forth on whether or not I should include a link list.
%Decided to keep it.  I need to exclude basal species consuming new species
%in this.  I don't think there is necessarily a clear benefit to passing or
%not passing it in.

S1 = length(n);

%lower border of niche range for every prey:
preymins = c - r/2;

%upper border of niche range for every predator:
preymaxs = c + r/2;

%fills the empty matrix with niche ranges:
n_mx = n*ones(1,S1);
%matrix with the lowest points of ranges in every column:
preymins_mx = ones(S1,1)*preymins';
%same, with highest:
preymaxs_mx = ones(S1,1)*preymaxs';

%Construct the web adjacency matrix; if species in the row is in the
%diet range of the species in the column, it gets eaten (if species i
%is eaten by species j, set (i,j) = 1).

web_mx=((n_mx>=preymins_mx)&(n_mx<=preymaxs_mx));


if N <= 0
    [res,con] = find(web_mx);
    noobs = [];
    return
end
    

OGbasal = (sum(web_mx)==0)';
nOGBasal = sum(OGbasal);

SFree = S1-nOGBasal;


%DOn't try to maintain current connectance, use the same parameters to
%generate the new species.  But note that the new species must be
%consumers, so the desired diet width needs to be smaller to account for
%that.

%Number of links to original free livers
L_OGConsumers = sum(sum(web_mx(:,~OGbasal)));

%Expected number of new links given expected coverage of niche space and
%additional species added
L_NewOGConsumers = L_OGConsumers + C0*SFree*N;

%THis is the generality we are going to try to match when we add new
%species; the species of the web nugget have a certain distribution of
%(absolute) generalities and I want the parasites to fit into that
%distribution, or at least have the same mean.  Will be worthwhile to check
%out the distribution of generalities of the original webs.
genOGConsumers = L_NewOGConsumers/SFree;

wNew = genOGConsumers/(S1+N);
S = S1;
Sp = 1:S1;
%%%Testing for connectivity of web first; adding consumers means web must
%%%be connected so this isn't necessary later on.  But, I'm thinking I do
%%%want ot make sure that the generalities of parasites is the same.
weak_comp = graphconncomp(sparse(web_mx),'Directed',1,'Weak',1);
web_connected = (weak_comp == 1);

basal = sum(web_mx)==0;

%make sure that all species have a directed connection from a basal
%species; only if connected.

if web_connected
    unconnected_species = (1:S)';
    
    
    %check that all energy flow starts at a basal
    connected_species = [];
    
    for kk = Sp(basal)
        connected_species = walk(kk,connected_species,web_mx);
    end
    
    unconnected_species(connected_species) = [];
    if isempty(unconnected_species)
        web_connected = 1;
    else
        web_connected = 0;
    end
end

if ~web_connected
    error('The web is not connected. \\n')
end
%This is a little funny since we aren't actually aiming for an overall web
%connectance, just that the mean generality of these consumers be the same
%as the consumers in the original web.  They will probably have slightly
%different vulnerabilities because we need to force the existing basal
%species to remain basal; the original web structure must be a stable core
%for the food web.  So, the check for web ok-ness includes connectivity and
%a match in generalities.


%I'm thinking it might be impossible, especially with smaller webs that end
%up very sparse, to add any number of species and still maintain
%the connectance.  So, this is just in case that happens.
if wNew <0
    error('Error: It''s not gonna work!.')
end


%error on the connectance we are trying to hitt
err = 0.05;

%validation of the niche web:
ok_n=0;

S = S1 + N;


Sp = 1:S;
cNoob0 = zeros(N,1);

tries = 1000;

while ((tries>0) && (ok_n==0))
    
    tries=tries-1;
    %assign noob species niche values according to a uniform distribution.
    %Slight bias introduced by making these species have higher niche
    %values (especially if low niche values got all destroyed by the
    %filtering process.  Not ideal, but probably necessary.  Also, might
    %have problems if filtered webs are not strongly connnected... SHould
    %exclude those.
    n1 = min(n);
    
    nNoob = rand(N,1)*(1-n1) + n1;
    
    %parameters for beta distribution:
    alpha = 1;
    beta = alpha/wNew/2 - alpha;
    
    yNoob = betarnd(alpha,beta,N,1);
    
    rNoob = yNoob.*nNoob;
    
    %The guaranteed resources of the noobs; but I don't want this to be the
    %center.  Think this is the most straighforward way; if the diet is
    %small, the possible centers could be several disjoint intervals >.<
    for jj = 1:N
        cNoob0(jj) = datasample(n(n<=nNoob(jj)),1);
    end
    
    %Now, we add the little bump:
    %This just ensures that the diet is within the niche space, and that
    %the center doesn't go above the niche value.
    bumpMins = max([-rNoob/2,rNoob/2-cNoob0],[],2);
    bumpMaxs = min([rNoob/2,1-(nNoob+rNoob/2),nNoob-cNoob0],[],2);
    
    bump = bumpMins + rand(N,1).*(bumpMaxs-bumpMins);
    cNoob = cNoob0 + bump;
    
    %sort the noobs:
    [nNoobSort, Indx] = sort(nNoob);
    %np_new: niche values in ascending order Indx: indices of species in
    %descending order
    
    rNoobSort = rNoob(Indx);                                                       %NK: Maybe not how I'd do this
    cNoobSort = cNoob(Indx);
    
    %Put the noobs in with free-livers
    n_new = [n;nNoobSort];
    c_new = [c;cNoobSort];
    r_new = [r;rNoobSort];
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %lower border of niche range for every prey:
    preymins = c_new - r_new/2;
    
    %upper border of niche range for every predator:
    preymaxs = c_new + r_new/2;
    
    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,S);
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(S,1)*preymins';
    %same, with highest:
    preymaxs_mx = ones(S,1)*preymaxs';
    
    %Construct the web adjacency matrix; if species in the row is in the
    %diet range of the species in the column, it gets eaten (if species i
    %is eaten by species j, set (i,j) = 1).
    
    web_mx=((n_mx>=preymins_mx)&(n_mx<=preymaxs_mx));
    
    %Make sure basal species stay basal species.
    basal = zeros(S,1);
    basal(1:S1) = OGbasal>0;
    
    web_mx(:,basal>0) = 0;
    
   
    
    %Check the generalities (web is connected by default).  But don't worry
    %about these if you are only adding a few species!  
    %10% of original species is chosen arbitrarily.
    
    if (web_connected)
        if N > ceil(.1*S1)
            
            genAll = sum(web_mx);
            genOG = genAll(1:S1);
            genOGCon = mean(genOG(~OGbasal));
            genNew = mean(genAll((S1+1):end));
            %Within Tolerance?
            if abs(genOGCon-genNew)/genOGCon > err
                ok_n=0;
            else
                ok_n=1;
            end
            
        else
            ok_n = 1;
        end
    end
end

if tries == 0
    error('Could not add species in 1000 tries.')
else
    
    [res,con] = find(web_mx);
    n = n_new;
    c = c_new;
    r = r_new;
    noobs = false(S,1);
    noobs((S1+1):end) = true;
    r(OGbasal) = 0;
    c(OGbasal) = 0;
end

end