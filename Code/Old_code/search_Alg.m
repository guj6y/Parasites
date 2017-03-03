function [conn] = search_Alg(Sources,Targets,N,source,ord)

Nodes = 1:N;


Pred = zeros(N,1);
Order = zeros(N,1);
Marked = zeros(N,1);
Marked(1) = source;
next = source;
Pred(1) = 0;
Order(1) = 1;
Active_node = 1;
List = 1;

Admissable_arc_targets = Targets(Targets~=source);
Admissable_arc_sources = Sources(Targets~=source);

Admissable_arcs =[Admissable_arc_sources,Admissable_arc_targets];

%First in First out
while ~(isempty(List))
    %Select a node i in list
    if ord
        Active_node = List(end); %First in First out
    else
        Active_node = List(1); %Last in First out
    end
    Active_admissable_targets...
        = Admissable_arc_targets(Admissable_arc_sources==Active_node);
    
    if ~isempty(Active_admissable_targets)
        Target = Active_admissable_targets(1);
        Admissable_arcs;
        Marked(Target) = 1;
        Admissable_arc_sources =...
            Admissable_arc_sources(Admissable_arc_targets~=Target);
        Admissable_arc_targets =...
            Admissable_arc_targets(Admissable_arc_targets~=Target);
        
        Admissable_arcs =[Admissable_arc_sources,Admissable_arc_targets];
        Pred(Target) = Active_node;
        next = next+1;
        Order(Target) = next;
        List = [List;Target]; %Last in First out
    
    else
    List(List==Active_node) = [];

    end
    
end

if length(Pred(Pred>0)) == (N-1)
    conn = 1;
else
    conn = 0;
end
    

end