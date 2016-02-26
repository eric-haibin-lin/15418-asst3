#include "paraGraph.h"
#include "graph.h"
#include <omp.h>


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
    Decomposition(Graph g, bool *visited, int *decomp, int* dus, int maxVal) 
    :num_nodes_(g->num_nodes), visited_(visited), decomp_(decomp), dus_(dus), iter_(0), maxVal_(maxVal)
    {
		pmemset(decomp_, INVALID_ID, num_nodes_);
    }

    bool update(Vertex src, Vertex dst) {
		int old_decomp_id = decomp_[dst];
		int new_decomp_id = decomp_[src] == INVALID_ID ? src : decomp_[src];
		// assign decomp_id to the smaller one
		while (old_decomp_id == INVALID_ID || old_decomp_id > new_decomp_id) {	
			// successful update new decomp_id
			if(__sync_bool_compare_and_swap(&decomp_[dst], old_decomp_id, new_decomp_id)) {
				break;
			}
			// failure, retry
			old_decomp_id = decomp_[dst];
		}
		return true;
    }

    // assign the center of vertex v to itself
    void assignCenter(Vertex v) {
    	decomp_[v] = v;
    }

    // return true if not visited yet
    bool cond(Vertex v) {
      	return !visited_[v];
    }

    bool operator()(Vertex v) {
    	//center is already visited
    	if (visited_[v] == true) {
    		return false;
    	}
      	if (iter_ > maxVal_ - dus_[v]) {
        	visited_[v] = true;
			assignCenter(v);
			return true;
      	}
      	return false;
    }
    
    void claimBalls() {
#pragma omp parallel for schedule(static, 512)
    	for (int i = 0; i < num_nodes_; i++) {
    		if (!visited_[i] && decomp_[i] != INVALID_ID) {
    			visited_[i] = true;
    		}
    	}
    }

    void debug() {
    	printf("\n============ visited array ========================= \n");
	    for (int i = 0; i < num_nodes_; i++) {
	    	//center is already visited
    		printf("%d ", visited_[i]);
	    }
		printf("\n============ dus array ========================= \n");
	    for (int i = 0; i < num_nodes_; i++) {
	    	//center is already visited
    		printf("%d ", dus_[i]);
	    }
    	printf("\nfinding qualified new center with iter %d maxVal %d... \n", iter_, maxVal_);
	    for (int i = 0; i < num_nodes_; i++) {
	    	//center is already visited
	    	if (visited_[i] == true) {
	    		continue;
	    	}
	      	if (iter_ > maxVal_ - dus_[i]) {			        	
	        	printf("found qualified new center %d \n", i);
	      	}
	    }
    }

    // increment to next iteration 
    void nextIterate() {
    	iter_ += 1;
    }

  private:
    int num_nodes_;
    bool *visited_;
    Vertex *decomp_;
    int* dus_;
    int iter_;
    int maxVal_;
};


inline bool *initVisitBitMap(int size) {
	bool *visited = (bool *) malloc(sizeof(bool) * size);
	// init values
#pragma omp parallel for schedule(static, 512)
    for(int i = 0 ; i < size; i++){
        visited[i] = false;
    }
    return visited;
}

// return a full vertex set 
inline VertexSet *initFullVertexSet(int num_nodes) {
	VertexSet *full_vertex_set = (VertexSet *)malloc(sizeof(VertexSet));
    Vertex *vertices = (int *)malloc(sizeof(int) * num_nodes);
#pragma omp parallel for schedule(static, 512)
    for(int i = 0 ; i < num_nodes; i++){
        vertices[i] = 1;
    }
    full_vertex_set->vertices = vertices;
    full_vertex_set->size = num_nodes;
    full_vertex_set->numNodes = num_nodes;
    full_vertex_set->type = DENSE;
    full_vertex_set->capacity = num_nodes;
    return full_vertex_set;
}

// return an initial frontier
inline VertexSet *initFrontier(int num_nodes, int maxId) {
	VertexSet *frontier = (VertexSet *)malloc(sizeof(VertexSet));
    Vertex *vertices = (int *)malloc(sizeof(int) * num_nodes);
#pragma omp parallel for schedule(static, 512)
    for(int i = 0 ; i < num_nodes; i++){
        vertices[i] = 0;
    }
    vertices[maxId] = 1;
    frontier->vertices = vertices;
    frontier->size = 1;
    frontier->numNodes = num_nodes;
    frontier->type = DENSE;
    frontier->capacity = num_nodes;
    return frontier;
}

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
	bool *visited = initVisitBitMap(num_nodes);

	// init full vertex set 
	VertexSet *full_vertex_set = initFullVertexSet(num_nodes);

	Decomposition decomposition(g, visited, decomp, dus, maxVal);

	//init frontier. vertex with maxDu grows first
	VertexSet *frontier = initFrontier(num_nodes, maxId);
	visited[maxId] = true;
	decomposition.assignCenter(maxId);

	while (frontier->size > 0) {
		// do BFS
		VertexSet *new_frontier = edgeMap<Decomposition>(g, frontier, decomposition);

		// claim visited balls
		decomposition.claimBalls();
		freeVertexSet(frontier);
		decomposition.nextIterate();

		// start growing all balls i at the next iter with 
    	// unvisited center i and with maxDu - dus[i] < iter
    	// so that addVertex can be performed in parallel
    	VertexSet *grown_frontier = vertexMap<Decomposition>(full_vertex_set, decomposition, true);
	    frontier = vertexUnion(grown_frontier, new_frontier);
	}
	freeVertexSet(frontier);
}