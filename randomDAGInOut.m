%{ 

Create a random weakly-connected Directed Acyclic Graph (DAG) on n nodes
with m edges with the additional constraint of x start-nodes and y exit-nodes.

Start-nodes are nodes with no in-edges, and exit nodes are nodes with no
out-edges, this also gaurantees that all nodes which are not start or exit
nodes will have both in and out edges.

Much like how an Erdos-Renyi graph can be created in a G(n, m) form, where n represents the number of 
nodes and m the number of edges, this algorithm assembles a simple weakly-connected DAG in much the same way.

This algorithm starts by creating a complete DAG on n nodes, that is the
upper triangle of the adjacency matrix.

Next it removes all in edges from the first x nodes, and all out edges of
the last y nodes.

Next the algorithm starts to randomly removes edges from that semi-complete
DAG until it arrives at the final configuration.

Because it is a destructive process which is randomly removing edges it can sometimes trap itself in
topologies that cannot be reduced to the desired number of edges, so it
also randomly adds edges if there are 2 generations in a row where an edge
cannot be removed.
%}


function dag = randomDAGInOut(nodes, edges, inCount, outCount)

    % get the number of middle nodes
    centerCount = nodes-inCount-outCount;

    % check to make sure that this graph is possible to construct
    if edges > (centerCount*(2*inCount + 2*outCount + centerCount - 1) + 2*inCount*outCount)/2
        % this formula describes the max number of edges possible on a
        % di-graph with this structure
        fprintf("Too many edges to create a DAG \n")
        dag = NaN;
        return
    elseif edges < nodes - 1
        fprintf("Not enough edges to make a basic weakly-connected graph on the given number of vertices. \n")
        dag = NaN;
        return
    elseif nodes < 2
        fprintf("At least two nodes are required to make a graph. \n")
        dag = NaN;
        return
    elseif(inCount+outCount)>nodes
        fprintf("There cannot be more in and out nodes than total nodes. \n")
        dag = NaN;
        return
    end
    
    % Generate all possible edges for a fully connected graph
    [s, t] = find(triu(ones(nodes), 1));
    
    % Create the graph
    G = digraph(s, t);
    
    % Delete al in edges to in nodes and out edges for out nodes
    for i = 1:nodes
        if i <= inCount
    
            % get in edges
            inNeighbors = inedges(G, i);
            G = rmedge(G,inNeighbors);
            
        end
    
        if i > nodes-outCount
    
            % get in edges
            outNeighbors = outedges(G, i);
            G = rmedge(G,outNeighbors);
        end
    end
    
    % Get the adjacency matrix of this full version
    fullAdjMatrix = full(adjacency(G));
    
    while height(G.Edges) > edges
        
        % pick a random edge to check
        randomIndex = randi(numedges(G));
        
        % Count the number of connected components in the original graph
        numComponentsOriginal = conncomp(G,'Type','weak');
        
        % Remove the specified edge
        G_modified = rmedge(G, randomIndex);
        
        % Count the number of connected components after removing the edge
        numComponentsModified = conncomp(G_modified,'Type','weak');
        
        lastTime = true;

        % Check if the edge is a bridge
        if numComponentsModified == numComponentsOriginal
    
            % start with assumption you will delete this edge
            deleteEdge = true;
    
            % Check if the first node can have this edge deleted
            G.Edges(randomIndex,:).EndNodes(1);
            outNeighbors = outedges(G, G.Edges(randomIndex,:).EndNodes(1));
            if length(outNeighbors) <= 1
                deleteEdge = false;
            end
    
            % Check if the second node can have this edge deleted
            G.Edges(randomIndex,:).EndNodes(2);
            inNeighbors = inedges(G, G.Edges(randomIndex,:).EndNodes(2));
            if length(inNeighbors) <= 1
                    deleteEdge = false;
            end 
    
            if deleteEdge == true
                G = rmedge(G,randomIndex);
            elseif deleteEdge == false
                if lastTime == false
                    newAdjMatrix = full(adjacency(G));
                    if ~isequal(fullAdjMatrix, newAdjMatrix)
                        differenceLocations = (fullAdjMatrix ~= newAdjMatrix);
                        [row, col] = find(differenceLocations);
                    
                        % pick a random edge to reconnect
                        randomRevert = randi(length(row));
                        row(randomRevert);
                        col(randomRevert);
                        newAdjMatrix(row(randomRevert),col(randomRevert)) = 1;
                        
                        % Set graph to new matrix
                        G = digraph(newAdjMatrix);
                    
                    end
                end
            end
            lastTime = deleteEdge;
    
        end
    
    end
    
    dag = full(adjacency(G));

end


% Function from DeltaCon
function similarityScore = deltaCon(graphIn1, graphIn2)
    
    A1 = full(adjacency(graphIn1,"weighted"));
    A2 = full(adjacency(graphIn2,"weighted"));
    nodeCnt = size(A2,1);

    p = 0.51;
    inv1 = inverse_LBP( A1, nodeCnt ) .* (p-0.5);
    
    inv2 = inverse_LBP( A2, nodeCnt ) .* (p-0.5);
    
    % Computing the DeltaCon similarity
    similarityScore = 1 / (1 + sqrt( sum(sum( (sqrt(inv1) - sqrt(inv2)).^2 ) )) ); 

end

% Function from DeltaCon
function inv_ = inverse_LBP(A, no_nodes)
    
    max_power = 7;

    % Create sparse
    I = speye(no_nodes);
    
    % Create the sparse degree-diagonal matrix
    D = sparse( diag(sum(A,2)) );
    
    % Compute the about-half homophily factor to guarantee convergence
    c1 = trace(D)+2;
    c2 = trace(D^2) - 1;
    h_h = sqrt((-c1+sqrt(c1^2+4*c2))/(8*c2));
    
    % Compute the constants ah and ch involved in the linear system
    ah = 4*h_h^2 /(1-4*h_h^2);
    ch = 2*h_h / (1-4*h_h^2);
    
    % Invert the matrices M1 and M2
    M = ch*A - ah*D;
    
    % Calculate the inverse of matrix M
    inv_ = I;
    mat_ = M;
    pow = 1;

    while max(max(mat_)) > 10^(-9) && pow < max_power
        inv_ = inv_ + mat_;
        mat_ = mat_*M;
        pow = pow +1;
    end
end