#include "paraGraph.h"
#include "graph.h"


#define INVALID_ID -1

void printDecomp(int *decomp, int width) {
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < width; j++) {
			printf("%d ", decomp[i * width + j]);
		}
		printf("\n");
	}
	printf("\n");
}

class Decomposition
{
  public:
    Decomposition(Graph g, bool *visited, int *decomp)
    {
      	//init pointers
		visited_ = visited;
		num_nodes_ = g->num_nodes;
		decomp_ = decomp;
		//TODO parallelize this part 
		for (int i = 0; i < num_nodes_; i++) {
			decomp_[i] = INVALID_ID;
		}
    }

    bool update(Vertex src, Vertex dst) {
		int old_decomp_id = decomp_[dst];
		int new_decomp_id = decomp_[src] == INVALID_ID ? src : decomp_[src];
		// assign decomp_id to the smaller one
		while (old_decomp_id == INVALID_ID || old_decomp_id > new_decomp_id) {	
			bool status = __sync_bool_compare_and_swap(&decomp_[dst], old_decomp_id, new_decomp_id);
			// successful update new decomp_id
			if (status) {
				break;
			}
			old_decomp_id = decomp_[dst];
		}
		return true;
    }

    bool cond(Vertex v) {
    	// return true if not visited yet
      	return !visited_[v];
    }

    bool operator()(Vertex v) {
		// if (visited_[i] == true) {
  //   		continue;
  //   	}
  //     	if (iter > maxVal - dus[i]) {
  //       	//add vertex
  //       	//TODO do this in parallel
  //       	addVertex(frontier, i);
  //     	}
      	return false;
    }
    
    void claimBalls() {
    	//TODO parallelize this part
    	for (int i = 0; i < num_nodes_; i++) {
    		if (!visited_[i] && decomp_[i] != INVALID_ID) {
    			visited_[i] = true;
    		}
    	}
    }

  private:
    int num_nodes_;
    bool *visited_;
    Vertex *decomp_;
};

/**
	Given a graph, a deltamu per node, the max deltamu value, and the id
	of the node with the max deltamu, decompose the graph into clusters. 
        Returns for each vertex the cluster id that it belongs to inside decomp.
	NOTE: deltamus are given as integers, floating point differences
	are resolved by node id order

**/
void decompose(graph *g, int *decomp, int* dus, int maxVal, int maxId) {
	// both the size of decomp and that of dus are numNodes
	int num_nodes = g->num_nodes;

	//init decomp instance
	bool *visited = (bool *) malloc(sizeof(bool) * num_nodes);
	//TODO use pmemset
	memset(visited, 0, sizeof(bool) * num_nodes);
	Decomposition decomposition(g, visited, decomp);

	//init frontier. vertex with maxDu grows first
	VertexSet *frontier = newVertexSet(SPARSE, num_nodes, num_nodes);
	visited[maxId] = true;
	decomposition.update(maxId, maxId);
	addVertex(frontier, maxId);
	int iter = 0;

	while (frontier->size > 0) {
		VertexSet *new_frontier = edgeMap<Decomposition>(g, frontier, decomposition);
		
		decomposition.claimBalls();
		
		free(frontier);
		frontier = new_frontier;
		
		iter++;

		// start growing all balls i at the next iter with 
    	// unvisited center i and with maxDu - dus[i] < iter
    	// TODO use vertexMap to parallelize this part 
    	// TODO use UnionSet to get the new frontier 
	    for (int i = 0; i < num_nodes; i++) {
	    	//center is already visited
	    	if (visited[i] == true) {
	    		continue;
	    	}
	      	if (iter > maxVal - dus[i]) {
	        	//add vertex
	        	addVertex(frontier, i);
	        	visited[i] = true;
				decomposition.update(i, i);
	      	}
	    }
	}
}
