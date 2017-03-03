function [connected] = walk1(source,connected,res,cons) 

%This function recursively finds all nodes connected to the source using 
%the link list

neighbors = cons(res==source)';

connected = [connected;source];
for kk = neighbors
    
    if sum(kk == connected)
        continue
    else
        
        connected = walk1(kk,connected,res(res~=source),cons(res~=source));
    end
    
    
end