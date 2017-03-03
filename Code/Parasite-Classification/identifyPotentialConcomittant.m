function [concomittant] = identifyPotentialConcomittant(res,con,para)

%This function creates a link list of potential concomittant links given a
%resource list,(res) a consumer list (con) and a parasite indicator vector.
concomittant = [];
S = max([res;con]);
spList = 1:S;

%List of parasites
parasites = spList(para);

for ii = parasites
    %List of hosts of that parasite
    hosts_ii = res(con==ii);
    for jj = hosts_ii'

        %skip jj if the 'host' is a parasite
        if sum(jj==parasites)
            continue
        else
            %otherwise, find consumers of that host and loop
            con_jj = con(res==jj);
            
            for kk = con_jj'
                %if there is already a link from this consumer to the first
                %parasite, remove from the list.
                if sum((con==kk)&(res==ii))
                    con_jj(con_jj==kk) = [];
                end
                %IF the consumer is a parasite, delete the consumer fromt
                %he list
                if sum(kk==parasites)
                    con_jj(con_jj==kk) = [];
                end
            end
            concomittant = [concomittant;...
            [ii*ones(size(con_jj)) con_jj ]];
        end
        
        
    end
end

mx = sparse(concomittant(:,1),concomittant(:,2),1,S,S);
[res1, con1] = find(mx);
concomittant = [res1 con1];

end







