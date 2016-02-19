#ifndef __VERTEX_SET__
#define __VERTEX_SET__

#include "graph.h"

typedef enum {
  SPARSE,
} VertexSetType;

//TODO add bitmap to this structure
typedef struct {
  int size;     // Number of nodes in the set
  int numNodes; // Number of nodes in the graph
  int capacity;
  VertexSetType type; 
  Vertex* vertices;
  int* vertices_bitMap;
  bool ifarray;
} VertexSet;

VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes);
void freeVertexSet(VertexSet *set);

void addVertex(VertexSet *set, Vertex v);
void removeVertex(VertexSet *set, Vertex v);
void printBitMap(VertexSet *set);
void printVertices(VertexSet *set);


VertexSet*  vertexUnion(VertexSet *u, VertexSet* v);

#endif // __VERTEX_SET__
