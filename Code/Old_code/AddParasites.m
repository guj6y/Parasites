%-----------------------------------------------

 function [nicheweb,Index_par,n,c,r]= ...
     AddParasites(nicheweb,n,c,r,num_parasites,connectances,varargin)  
num_species = length(n);
%------------------------------------------------------------------
if nargin == 6
    %must specify to have 2 axis version
    %version decides one or two axes.
    version =  1;
elseif nargin==7
    version = varargin{4};
elseif nargin == 8
    version = varargin{4};
    error_p = varargin{5};
end

Cff = connectances(1);
Cpf = connectances(2);

if numel(connectances) ==3
    Cpp = connectances(3);
elseif numel(connectances) == 4
    Cpp = connectances(3);
    Cfp = connectances(4);
end

%Cff is connectance between free-livers
%Cfp is connectance of predation on parasites
%Cpp is connectances of hyperparasitism
%Cpf is connectance of parasitism

%error on the connectances for parasites
error_p=0.05;

%nicheweb with incidental predation to calculate corresponding properties
%web_ip=0; 

%index of the parasites
Index_par=[];

%validation of the plain niche web:
ok_n=1; 

%validation of the fact that every parasites has a host. 
%0 means that it is false by default
ok_par_host=0; 

%Now we can add the parasites, following the same scheme

%web_without_parasites=nicheweb;
too_high_par = 0;
too_low_par = 0;
tries=10000;
ok_p=0;

%setting parameters for verison 1.  need to think about importance of
%restricting np_min and np_max.

%version 1 = one axis;

%match two connectances; Cff, Cpf
if (version == 1) && (numel(connectances) == 2)
    fprintf('a\n')
    np_min = 0;
    np_max = 1;
    c_pmin = 0;
    P=1;
end

%also match Cpp
if (version == 1) && (numel(connectances) == 3)
    fprintf('b\n')
    np_min = 0;
    np_max = 1;
    c_pmin = np_min+(np_max-np_min)*(1-Cpp/Cpf);
    P=1;
end

%also match Cfp; lose intervality
if (version == 1) && (numel(connectances) == 4)
    fprintf('c\n')
    np_min = 0;
    np_max = 1;
    c_pmin = np_min+(np_max-np_min)*(1-Cpp/Cpf);
    %probability of accepting predation on a parasite
    P = Cfp/Cff;
end
    
%Primary loop for adding parasites.  As above, the parasites are drawn en
%masse.  Unlike above, there is a secondary loop that ensures all parasites
%have a host before checking the connectance and deciding whether or not to
%throw out the drawn parasites.
n_new = n;
c_new = c;
r_new = r;

while ((num_parasites>0) && (ok_n==1) && (tries>0) && (ok_p==0))
    tries=tries-1;
    ok_p=1;
    n_par=rand(num_parasites,1)*(np_max-np_min)+np_min;

%%%%%%%%%%%%    designate range for each parasite
    
    
%parameters for beta distribution:
    alpha = 1;
    param=2*Cpf/(1-np_min)/(1-c_pmin);
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
        2*ones(num_species+num_parasites));                                

    %accept only P% of 'free liver <- parasite' links
    P_mx = rand(num_parasites,num_species);
    
    P_mx = P_mx<=P;
    
    web_mx(Index_par,1:num_species) = ...
        (P_mx + web_mx(Index_par,1:num_species)) == 2;
    
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
        
            %accept only P% of 'free liver <- parasite' links
        P_mx = rand(num_parasites,num_species);
    
        P_mx = P_mx<=P;
    
        web_mx(Index_par,1:num_species) = ...
            (P_mx + web_mx(Index_par,1:num_species)) == 2;
        
   
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
   
    
    %[web_ip] = Find_Incidental_Predation_nk(web_ip,nichewebsize,Index_par);
    %[web_ip] = Find_Incidental_Predation_nk2(web_ip,Index_par);
    links=sum(web(Index_par,1:num_species)>0);
    % Actual connectance   
    C_par_web = sum(links)/(num_parasites*num_species);
    %
    if abs(C_par_web-Cpf)*1.0/Cpf>error_p
        ok_p=0;
        if C_par_web-Cpf >0
            too_high_par = too_high_par+1;
        else
            too_low_par = too_low_par+1;
            C_par_web
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
tries

end

