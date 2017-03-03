function [num_chn,num_scc_end,num_scc_in,num_scc_and,chn_lens,down_paths,loop_num] =  count_chains(res,con)

%Function for counting the number of paths from a basal species to a top
%species; but not all webs have top species!  In testing I found that these
%tended to have strongly connected top communities - i.e. species mutually
%preying on each other (and usually even more complicated sets of
%interactions between 4 or 5 species).  In the end I decided to count 4
%things:
%
%1. Number of paths from a basal to a top (discounting cannibalism) species
%2. Number of paths ending in a strongly connected community of 'apex'
%predators.
%3. Number of paths passing through a strongly connected sub-community.
%4.  Number of paths passing through and ending at a strongly connected
%sub-community.
%
%This seemed to be a reasonable thing to do.  We will see what the
%dynamical results say; it could be that adding parasaites (probably,
%creating many more strongly connected subcommuities) will make this whole
%thing metric somewhat pointless.
%
%
% Cannibals count as the end of a food chain.
%
% May also want to think aboutthe statistics of the strongly connected
% components - "dynamics of strongly connected components of consumer
% resource food web models"
%
% For diagnostics on a graph:
%global h
%h = view(biograph(A));
%h.Nodes(i).color = [1 0 0] change sthe color of node i to red.  I used
%this to watch the algorithm go; was able to convince myself it was working
%as I expected.


%No  cannibals
res0 = res(res~=con);
con0 = con(con~=res);
res = res0;
con = con0;

%Number of species:
N = max([res;con]);
sp_list = 1:N;

%initializing some things...
%Form adjacency matrix for SCC detection (or for diagnostics);
A = sparse(res,con,1,N,N);

%Number of downstream paths of each type from this node.
down_paths = zeros(N,4)-inf;


[n_comp, comps] = graphconncomp(A,'Directed',true,'Weak',false);

num_in_loops = 0;
sccs = zeros(n_comp,1);

%Collapsing the strongly connected components (a strongly connected
%component could consist of a single node >.<

res_new = zeros(size(res));
con_new = zeros(size(con));

for kk = 1:comps
    num_in_ii = sum(kk==comps);
    if num_in_ii> 1 %If the SCC is nontrivial
        %Everyone in the SCC is in a loop
        num_in_loops = num_in_loops + num_in_ii;
        %'new' node ii is a scc - want to track this properly.
        sccs(kk) = true;
    end
    %Re-labeling.
    for jj = sp_list(kk==comps);
        res_new(res==jj) = kk;
        con_new(con==jj) = kk;
    end
end
% Getting rid of redundant links.  might be possible with 'unique', but I
% don't trust it.
A_new = sparse(res_new,con_new,1,n_comp);
[res,con] = find(A_new);

N = n_comp;
L = length(res);

%Logical if each species is a top species
tops = sum(A,2)==0;

%Recording the state of the current path
in_path = zeros(N,1);


%}

%Recording how many times each edge is traversed:
% down_edges = zeros(L,1);
%Harder to do than I originally thought;  not so bad, just need to store
%the edges in the downstream paths from each node.  On balance, storing an
%L by 1 binary vector for each node might be a bit excessive.  *possibly* a
%more clever way to do it?  I just didn't see it.


%Number of chains of each type:
%1 = chains
%2 = chains endig in scc
%3 = chains containing non-terminal scc
%4 = chains passing thorugh and endig in a scc.
num_path_all = zeros(4,1);


%Indicator for the current path type; in get_nbrs this is added to the
%count each time we hit a terminal node.
path_type = zeros(4,1);

%number of each type of chain when first entering a node, sum of lengths,
%and sum-square of lengths
num_path0 = zeros(N,4);


basal = setdiff(res,con);
tops = setdiff(con,res);

    function [] = get_nbrs(cur_sp)
        
        %Set of (upstream) neighbors:
        nbrs = con(res==cur_sp);
        
        %total number of paths (of each type) before entering this species.
        num_path0(cur_sp,:) = num_chn_all;
        
        %Update the current path length
        in_path(cur_sp) = max(in_path) + 1;
        
        for ii = nbrs'
            if tops(ii) == true
                if sccs(ii) == true
                    path_type(2) = true;
                    path_type(1) = false;
                end
                num_path_all = num_path + path_type;
                
                
            else  %If it's not in the path or a top 'species'
                %Do it all over again for this neighbor.
                if sccs(ii) == true
                    path_type(3) = true;
                    path_type(1) = false;
                end
                get_nbrs(ii);
            end
            
        end
        
        if (sccs(cur_sp) == true)&&(tops(cur_sp)==true)
            
        end
        %Number of downstream paths for the node whose neighbors we just
        %visited:
        down_paths(cur_sp) = num_path_all - num_path0(cur_sp);
        
        
        %Take the current species out of the path so it can be included in
        %other chains; only after we have gone through all of its neighbors.
        in_path(cur_sp) = 0;
    end





for kk = basal'
    get_nbrs(kk);
end

end