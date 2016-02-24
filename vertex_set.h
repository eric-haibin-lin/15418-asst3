#ifndef __VERTEX_SET__
#define __VERTEX_SET__

#include "graph.h"
#include <omp.h>
#include <list>

using namespace std;


#define POLICY schedule(dynamic, 512)

typedef enum {
    SPARSE,
    DENSE,
} VertexSetType;

//TODO add bitmap to this structure
typedef struct {
    int size;     // Number of nodes in the set
    int numNodes; // Number of nodes in the graph
    int capacity;
    VertexSetType type; 
    Vertex* vertices;
} VertexSet;

class MemoryManager{
    public:
        int *bufferArray;
        list<VertexSet*> setList;

        MemoryManager(): bufferArray(NULL){};
        ~MemoryManager();
};


VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes);
void freeVertexSet(VertexSet *set);

void addVertex(VertexSet *set, Vertex v);
void removeVertex(VertexSet *set, Vertex v);
void printBitMap(VertexSet *set);
void printVertices(VertexSet *set);

void transform(VertexSet * set);
void pmemset(int * start, int val, int size);

VertexSet*  vertexUnion(VertexSet *u, VertexSet* v);

#endif // __VERTEX_SET__
