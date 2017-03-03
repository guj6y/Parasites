%-----------------------------------------------
% Program by: Rosalyn Rael

% Barbara Bauer added a part removing species not connected to a basal
% and changed the output to a 'row eats column' type matrix

% Coralie Picoche added the parasites and changed the 'removing part' 
% in order not to change the connectance

% Nick Kappler made significant modifications:
% Coding modifications (significantly improved efficiency of code)
% -added separate function for finding incidental predation to assist 
% profiling
% -Changed how parasites were handled; instead of being sorted into the
% niche range right away they are left at the end of the niche vectors (
% and center and diet vectors) then sorted at the end.  
%-Made range of parasite niche unrestricted and made parasite diets omit
%other parasites.

% Ref: Williams and Martinez, Nature, 2000.
% Ref: Warren et al. 2009
%---------------------------------------------------
% This function produces a niche model food web with
% Input: number of species, connectance, number of parasites, 
% connectance of parasites, fraction of free-living (??)
% Output: adjacency matrix called 'nicheweb' 
%    A(i,j) = 1 if i eats j.(row eats column)                              
%    
%--------------------------------------------------
%
 function [nicheweb,web_ip,Index_par,n,c,r,C_web]= ...
     NicheModelWithParasites_nk2(num_species, connectance, num_parasites,...
     connectance_parasites)                                                %NK: Make ouput a link list (or options)?  anyway, Make option to save as 
 %stream = RandStream.getGlobalStream;
 %reset(stream);                                                                           %.gdf (or text) file.
%------------------------------------------------------------------
% do a loop until the right properties are found in the web. This part
% only builds the normal niche model

tries=10000;                                                               %NK: Also consider implementation using sparse matrices.
%error on the connectance for free-livers
error_n=0.05;
%error on the connectance for parasites
error_p=0.05;
nicheweb=0;

%nicheweb with incidental predation to calculate corresponding properties
%web_ip=0; 

%index of the parasites
Index_par=[];

n=[]; %niche values
c=[]; %center values
r=[]; %range values

%validation of the plain niche web:
ok_n=0; 

%validation of the fact that every parasites has a host. 
%0 means that it is false by default
ok_par_host=0; 

while ((num_species>0) && (tries>0) && (ok_n==0))
    tries=tries-1;
    ok_n=1;
    %assign niche values from a uniform distribution
    n = rand(num_species,1);  
    
    
    %designate range for each species

    %parameters for beta distribution:
    alpha = 1;
    beta = (1-2*connectance)/(2*connectance); 
    r = betarnd(alpha,beta,num_species,1); 
    
    %vector of ranges:
    r = r.*n;  
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % set center of range, uniformly distributed in [r_i/2,n_i]; 
    c = min(1-r./2,rand(num_species,1).*(n-r./2)+r./2);

    %sort everything:
    [n_new, Indx] = sort(n);                                                
    %n_new: niche values in ascending order
    %Indx: indices of species in descending order 
    %(-> 1 is the index of the smallest niche range, 10 is the index of the
    %largest)

    %the smallest r to highest index species 
    r_new = r(Indx);                                                       %NK: Maybe not how I'd do this
    c_new = c(Indx); 

    r_new(1) = 0; %change the r of highest index species to 0
    %so we have a basal species in every web
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    
    
    %lower border of niche range for every prey:
    preymins = c_new - r_new/2;
    
    %upper border of niche range for every predator:
    preymaxs = c_new + r_new/2;

    %fills the empty matrix with niche ranges:
    n_mx = n_new*ones(1,num_species); 
    %matrix with the lowest points of ranges in every column:
    preymins_mx = ones(num_species,1)*preymins'; 
    %same, with highest:
    preymaxs_mx = ones(num_species,1)*preymaxs';

    %Construct the web adjacency matrix;
    %if species in the row is in the diet range of the species in the 
    %column, it gets eaten (if species i is eaten by species j, set (i,j) =
    %1).
    
    web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==2*ones(num_species));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  if there is an isolated species                                    %
    %  or something not connected to a basal                              %
    %  or a disconnected group                                            %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %1. removal if isolated
    
    %transposing the matrix (in the following codes we need that format)
    nicheweb = web_mx'; 

    nichewebsize=size(nicheweb,2);
    
    %find zeros
    testmx1 = (nicheweb==zeros(nichewebsize)); 
    
    %columns with only zero (they are not eaten by anyone)
    zerocols = find(sum(testmx1)==nichewebsize); 
    
    %rows with only zero (they don't eat anyone)
    zerorows = find(sum(testmx1,2)==nichewebsize);
    
    %both: index of isolated species
    qq = intersect(zerorows,zerocols); 
   
    if numel(qq)~=0
        ok_n=0;
    else                                                                   %NK: They are redrawing every variable for all species if something goes 
        %Save the web                                                      %wrong.  Is there a more efficient way to do this?  Can we detect groups  
        web_mx=nicheweb;                                                   %as we go?I.e. tracking the number of groups as we add diets? 
                                                                           % Concerns are: 1.efficiency, 2.elegance, 3.redability and fcv 
        %2. removal if species not connected to a basal                    %4.intuitiveness
        nichewebsize=size(nicheweb,2);
        basalsp=find(sum(nicheweb,2)==0);

        %identify list1
        list1=basalsp;
        list2=[];
        
        %This is almost like finding spanning trees from basal species.
        while numel(list1)>0                                               
        %look for connections; loop over 'basal species' - those in list1
            for whichspec=1:length(list1)
                %Determine which species eat the 'basal species'
                whichcol=list1(whichspec);
                %add those to list2
                list2=[list2; find(nicheweb(:,whichcol)==1)];
                list2=unique(list2);
            end
 
        %make species connected to list1 zeros (isolates already removed)
            nicheweb(list1,:)= zeros(length(list1),nichewebsize); 
            nicheweb(:,list1) = zeros(nichewebsize,length(list1));
    
        %make list2 the new list1
            list1=list2;
            list2=[];
        end

        %find indices of species not connected to basal
        list3=find(sum(nicheweb)>0);
        list3=[list3 find(sum(nicheweb,2)>0)'];
        list3=unique(list3);
        
        %putting the web with no isolates and not connected to basals to 
        %nicheweb
        nicheweb = web_mx; 
        
        %find if there are disconnected groups
        web_connected = isConnected(nicheweb);
        
        %if there are species not connected to basal, delete the web
        if numel(list3)>0
            ok_n=0;
            
        %if there are disconnected groups, delete the web
        elseif ~web_connected
            ok_n=0;
        
        %if all species are connected to a basal species, and there are no
        %disconnected groups, check the connectance.
        elseif web_connected
            % indices of links (which element in web_mx?)
            links = sum(web_mx>0);  
            % Actual connectance
            C_web = sum(links)/(num_species^2);  
            if (abs(C_web-connectance)*1.0/connectance) > error_n
                ok_n=0;
            else
                
                ok_n=1;
            end
        end  
        
        
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Now we can add the parasites, following the same scheme

%web_without_parasites=nicheweb;
too_high_par = 0;
too_low_par = 0;
tries=10000;
ok_p=0;

fraction_fp = .6;
fraction_pp = .3;
np_min = 0.1;
np_max = np_min+fraction_fp;
c_pmin = np_min+(np_max-np_min)*(1-fraction_pp);
%Primary loop for adding parasites.  As above, the parasites are drawn en
%masse.  Unlike above, there is a secondary loop that ensures all parasites
%have a host before checking the connectance and deciding whether or not to
%throw out the drawn parasites.
while ((num_parasites>0) && (ok_n==1) && (tries>0) && (ok_p==0))
    tries=tries-1;
    ok_p=1;
    n_par=rand(num_parasites,1)*(np_max-np_min)+np_min;

%%%%%%%%%%%%    designate range for each parasite

%parameters for beta distribution:
    alpha = 1;
    param=connectance_parasites/(1-np_min-fraction_fp/2)/(1-c_pmin);
    betap = (1-param)/param;
    r_p = betarnd(alpha,betap,num_parasites,1); %vector of ranges
    
    %Scale the ranges: lower niche value = more generalist.
    
    r_p = r_p.*(1-n_par);
%Rules for scaling r. If above isn't working, can try this.
%{
%first rule : if n_p increases, r_p increases (the bigger you are, the more
%generalist you are, Williams & Martinez 2000)
    if rule_r==1
        r_p = r_p.*n_par; 
    elseif rule_r==2  
%seconde rule : if n_p increases, r_p decreases (the bigger you are, the
%more specialist you are (Warren et al 2010))
%means that r can't be r_min and r can't be r_max either
        r_p=r_p.*(np_max-n_par+np_min);

%third rule : r is independant from n_p
    else
        r_p=r_p.*(np_max+np_min)/2;
    end
%}
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    set center for each parasites, according to the rules
c_p=c_pmin+r_p./2+rand(num_parasites,1).*(1-r_p/2-n_par);
%Rules for c. If above isn't working, can try this.
%{
    %first rule : allows a lot of hyperparasitism : c_p=>np_min+rp_2
    if rule_c==1
        c_p=np_min+r_p./2+rand(num_parasites,1).*(1-r_p-np_min);
    
    %second rule : allows no hyperparasitism : cp=>np_max+rp_2    
    elseif rule_c==2    
        c_p=np_max+r_p./2+rand(num_parasites,1).*(1-r_p-np_max);
    
    %third rule : c_p can be anywhere
    else
        c_p=r_p./2+rand(num_parasites,1).*(1-r_p);
    end;
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%    Construct the web adjacency matrix

    n=[n_new ; n_par];
    c=[c_new ; c_p];
    r=[r_new ; r_p];
    Index_par = (num_species+1):(num_species+num_parasites);
    %[n,Indx] = sort(n); 
    %n: niche values in ascending order
    %Indx: indexes of species in descending order 
    %(-> 1 is the index of the smallest niche range, 10 is the index of the
    %largest)

    %r = r(Indx); %the smallest r to highest index species
    %c = c(Indx); 
    %lower border of niche range for every prey:
    preymins = c - r/2; 
    %upper border of niche range for every predator:
    preymaxs = c + r/2; 

    %fills the empty matrix with niche ranges
    n_mx = n*ones(1,num_species+num_parasites);  
    %matrix with the lowest points of ranges in every column
    preymins_mx = ones(num_species+num_parasites,1)*preymins';
    %same, with highest
    preymaxs_mx = ones(num_species+num_parasites,1)*preymaxs'; 

    
    %if species in the row is in the niche range of the species in the 
    %column, it gets eaten
    web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==...
        2*ones(num_species+num_parasites));                                %NK: The reason for restricting the diets of parasites seems to be the 
                                                                           %adherence to this method of calculating the adjacency matrix.  Surely
                                                                           %there must be a way to simply skip over the parasitic rows/columns.  to
                                                                           %Avoid hyperparasitism/predation on parasites.
                                                                           
                                                                           
   %Removing predation on parasites and hyperparasitism; per warren et al.
   %web_mx(Index_par,:) = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %if there is a parasite that does not eat a host, it is replaced by
    %another parasite

    %Finding the Indices of the parasites in the parasites+freelivers
    %nichespace
    
    %[Index_par] = Find_Index_par_nk(n,n_par);
    %{
    Index_par = [];
   
    for i=1:length(n_par)
        j=1;
        trouve=0;
        while (trouve==0) && (j<length(n)+1)
            if n_par(i)==n(j)
                trouve=1;
                Index_par=[Index_par j];
            end;
            j=j+1;
        end;
    end;
    %}
    
    %transposing the matrix (in the following codes we need that format)
    nicheweb = web_mx'; 
    nicheweb_par=nicheweb(Index_par, 1:num_species);
    
    %looking for parasites that have no host
    testmx1_par = (nicheweb_par==...
        zeros(num_parasites,num_species));
    zerorows_par = find(sum(testmx1_par,2)==num_species+num_parasites);
    
    
    if numel(zerorows_par)~=0
        %some parasites have no host.
        ok_par_host=0;
    else
        %all parasites have a host
        ok_par_host=1;
    end;    
    tries_bis=10000;
    
    %This loop re-draws all parasites that don't have a host.
    
    while (tries_bis >0) && (ok_par_host==0)
        tries_bis=tries_bis-1;
        
         n_new_par=rand(numel(zerorows_par),1).*(np_max-np_min)+np_min;
         r_new_par = betarnd(alpha,betap,numel(zerorows_par),1)...
             .*(1-n_new_par);
         c_new_par=c_pmin+r_new_par./2 ...
         +rand(numel(zerorows_par),1).*(1-r_new_par/2-n_new_par);
        
        %Choosing diet ranges for re-drawn parasites
        %{
        if rule_r==1
            r_new_par = r_new_par.*n_new_par; 
        elseif rule_r==2  
            r_new_par=r_new_par.*(np_max-n_new_par+np_min);
        else
            r_new_par=r_new_par.*(np_max+np_min)/2;
        end;
        %}
        %Choosing diet centers for re-drawn parasites
        %{
        if rule_c==1
            c_new_par=np_min+r_new_par/2+rand(numel(zerorows_par),1)...
                .*(1-r_new_par-np_min);
        elseif rule_c==2
            c_new_par=np_max+r_new_par/2+rand(numel(zerorows_par),1)...
                .*(1-r_new_par-np_max);
        else
            c_new_par=...
                r_new_par/2+rand(numel(zerorows_par),1).*(1-r_new_par);
        end;
        %}
        
        %replacing parasites that need to be redrawn in parasite vectors
        n_par(zerorows_par)=n_new_par;
        c_p(zerorows_par)=c_new_par;
        r_p(zerorows_par)=r_new_par;
    
        %adding parasites to accepted free-liver vectors
        n=[n_new ; n_par];
        c=[c_new ; c_p];
        r=[r_new ; r_p];

        %Sorting species smalles to highest (parasites and free livers)
        %[n,Indx] = sort(n); 
        %r = r(Indx); 
        %c = c(Indx); 

        %lower border of niche range for every prey
        preymins = c - r/2; 
        %upper border of niche range for every predator
        preymaxs = c + r/2; 
    
        %fills the empty matrix with niche ranges
        n_mx = n*ones(1,num_species+num_parasites); 
        
        %matrix with the lowest points of ranges in every column
        preymins_mx = ones(num_species+num_parasites,1)*preymins'; 
        
        %same, with highest
        preymaxs_mx = ones(num_species+num_parasites,1)*preymaxs'; 
        
        %Building adjacency matrix
        web_mx=((n_mx>=preymins_mx)+(n_mx<=preymaxs_mx)==...
            2*ones(num_species+num_parasites));
        
        %Removing predation on parasites and hyperparasitism; per warren et al.
        %web_mx(Index_par,:) = 0;
   
        %Checking for parasites that have no hosts...
        
        %Retrieving indices of parasites
        Index_par = (num_species+1):(num_species+num_parasites);
        %[Index_par] = Find_Index_par_nk(n,n_par);
    
        %{
        Index_par = [];

        for i=1:length(n_par)
            j=1;
            trouve=0;
            while (trouve==0) && (j<length(n)+1)
                if n_par(i)==n(j)
                    trouve=1;
                    Index_par=[Index_par j];
                end;
                j=j+1;
            end;
        end;
        %}

        %transposing the matrix (in the following codes we need that format) 
        nicheweb = web_mx'; 
        nicheweb_par=nicheweb(Index_par, 1:num_species);
        
        %Parsaites with no host
        testmx1_par = (nicheweb_par==...
        zeros(num_parasites,num_species));
        zerorows_par = find(sum(testmx1_par,2)==num_species+num_parasites); 
        if numel(zerorows_par)~=0
            ok_par_host=0;
        else
            ok_par_host=1;
        end;
    end;
    web=nicheweb;
    %web_ip=nicheweb;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%At this point, all parasites have hosts

    % Adding incidental predation : warning : it is calulated another way
    % in the functional response : web_ip is only used for the properties !
    
    %i is the prey, j is the predator, %k is the parasite of i that gets 
    %incidentally eaten by the predator j
    %{
    for i=1:nichewebsize                                                   %NK: Warren et al. suggest only adding a (random) percentage of these.
        for j=1:nichewebsize                                               %NK: This loops N^2 * P times (.4s!!!).  Better way to do this?
            for k=Index_par 
                if web_ip(j,i)==1 && web_ip(k,i)==1
                    web_ip(j,k)=1;
                end
            end    
        end
    end
    %}
    
    %[web_ip] = Find_Incidental_Predation_nk(web_ip,nichewebsize,Index_par);
    %[web_ip] = Find_Incidental_Predation_nk2(web_ip,Index_par);
    links=sum(web(Index_par,1:num_species)>0);
    % Actual connectance   
    C_par_web = sum(links)/(num_parasites*num_species);
    %
    if abs(C_par_web-connectance_parasites)*1.0/connectance_parasites>error_p
        ok_p=0;
        if C_par_web-connectance_parasites >0
            too_high_par = too_high_par+1;
        else
            too_low_par = too_low_par+1;
        end
    else
        ok_p=1;
    end
    %}
    if ok_p==1
    
        links=sum(web(:,:)>0);
    
        % Actual connectance of the whole web:
        C_web = sum(links)/((num_parasites+num_species)...
            *(num_parasites+num_species));  
        %{
        connectance_tot=0.11;
        if abs(C_web-connectance_tot)*1.0/connectance_tot>error_p
            ok_p=0;
            if C_web-connectance >0
                too_high = too_high+1;
            else
                too_low = too_low+1;
            end
        else
            ok_p=1;
        end
        %}
    end
   
end
%Catching any failure.
if (num_parasites==0 && ok_n==0)                                           %NK: How about we do this earlier.. avoid superfluous loop iterations
    nicheweb=0;
elseif (num_parasites>0 && (ok_n==0 || ok_p==0 || ok_par_host==0))        
    nicheweb=0;
elseif (num_species==0 && (ok_p==0 || ok_par_host==0))
    nicheweb=0;
end
%[~,Indx] = sort(n);
%nicheweb = nicheweb(Indx,:);
%nicheweb = nicheweb(:,Indx);
                                                                           %About 150s for Carpinteria
end

