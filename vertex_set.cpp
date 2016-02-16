#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include "mic.h"

/**
 * Creates an empty VertexSet with the given type and capacity.
 * numNodes is the total number of nodes in the graph.
 * 
 * Student may interpret type however they wish.  It may be helpful to
 * use different representations of VertexSets under different
 * conditions, and they different conditions can be indicated by 'type'
 */
VertexSet *newVertexSet(VertexSetType type, int capacity, int numNodes)
{
  // TODO: Implement
  VertexSet *new_vertex_set = (VertexSet *)malloc(sizeof(VertexSet));
  new_vertex_set->size = 0;
  new_vertex_set->numNodes = numNodes;
  new_vertex_set->type = type;
  new_vertex_set->vertices = (Vertex *) malloc(sizeof(Vertex) * capacity);
  return new_vertex_set;
}

void freeVertexSet(VertexSet *set)
{
  // TODO: Implement
  // free the vertices before the set is freed
  free(set->vertices);
  free(set);
}

/**
 * Returns the union of sets u and v. Destroys u and v.
 */
VertexSet* vertexUnion(VertexSet *u, VertexSet* v)
{
  // TODO: Implement

  // STUDENTS WILL ONLY NEED TO IMPLEMENT THIS FUNCTION IN PART 3 OF
  // THE ASSIGNMENT

  return NULL;
}

