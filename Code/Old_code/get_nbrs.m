function [num_chn,num_scc_end,num_scc_in,num_scc_and,in_path,chn_lens,down_paths] = ...
    get_nbrs(num_chn,num_scc_end,num_scc_in,num_scc_and,cur_sp,res,con,in_path,chn_lens,down_paths)

global h tops sccs

%'exhaustive' depth first search.  Finds **all** possible paths from the 
%input node to sink nodes.  Doesn't go through loops.  This is practically 
%guaranteed to have a very bad scaling for large networks.  On the plus 
%side, the recursion level cannot go above the largest possible path length
%of $#species - #basal$.

% concern is the size of  chn_lens getting out of hand.  
% Consider updating sum of chain lengths and sum of squared chain lengths
% for mean and standard deviations of the lengths.

%Set of neighbors:
h.Nodes(cur_sp).color = [1 0 0];
pause()
nbrs = con(res==cur_sp);

num_chn0(cur_sp) = num_chn; %original number of paths
%Update the current path length
in_path(cur_sp) = max(in_path) + 1;
isempty(nbrs)
if isempty(nbrs) %if we reach a top species...
    h.Nodes(cur_sp).color = [0.5 0.5 0.5];
    pause();
    %increment the count of the chains
    num_chn = num_chn + 1;
    %save the path length
    chn_lens = [chn_lens; in_path(cur_sp)];
    down_paths(cur_sp) = 0;
    tops(cur_sp) = 1;
    h.Nodes(cur_sp).ID = sprintf('Top Species, %u',cur_sp);
    pause()
    in_path(cur_sp) = 0;
else %If we aren't at a top species
    for ii = nbrs'
        if in_path(ii) > 0 %Check if the neighbor is in the path (loop!)
            h.Nodes(ii).color = [0 0 0];
            pause()
            loop_num = loop_num + 1;
            in_loop(ii) = 1;
            h.Nodes(ii).color = [1 1 0.7];
            pause()
            continue
        elseif down_paths(ii) >= 0 %Check if we know how many chains from ii
            h.Nodes(ii).color = [1 0 1];
            pause()
            num_chn = num_chn + down_paths(ii);
           h.Nodes(ii).color = [0 0 1];
            pause()
            continue
        elseif tops(ii) == 1
            h.Nodes(ii).color = [1 1 0];
            pause();
            num_chn = num_chn+1;
            h.Nodes(ii).color = [0.5 0.5 0.5];
            pause();
        else  %If it's not in the path...
            %Do it all over again for this neighbor.
            h.Nodes(cur_sp).color = [1,.9,.9];
            pause()
            [num_chn,num_scc_end,num_scc_in,num_scc_and,in_path,edge_visits,num_chn,chn_lens,down_paths,tops,loop_num,in_loop] =...
                get_nbrs(num_chn,num_scc_end,num_scc_in,num_scc_and,ii,res,con,in_path,edge_visits,num_chn,chn_lens,down_paths,tops,loop_num,in_loop);
        end
        
    end
    
down_paths(cur_sp) = num_chn - num_chn0(cur_sp);
h.Nodes(cur_sp).color = [0 0 1];
pause()
h.Nodes(cur_sp).ID = sprintf('Node %u, paths through: %u',cur_sp,down_paths(cur_sp));    
end

%number of new paths since entering cur_node.
    

%Take the current species out of the path so it can be included in
    %other chains; only after we have gone through all of its neighbors.
    %note that this automatically skips self-loops.
    in_path(cur_sp) = 0; 
end