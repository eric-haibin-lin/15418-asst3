#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include "mic.h"
#include <ctime>

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
    // One change here is we donot malloc record part when we initialize a new VertexSet;
    // One draw back is we have to allocate array ourside of the class(because I donot want function...)
    VertexSet *new_vertex_set = (VertexSet *)malloc(sizeof(VertexSet));
    new_vertex_set->size = 0;
    new_vertex_set->numNodes = numNodes;
    new_vertex_set->type = type;
    new_vertex_set->capacity = capacity;
    new_vertex_set->ifarray = true;
    new_vertex_set->vertices = NULL;
    new_vertex_set->vertices_bitMap = NULL;
    //new_vertex_set->vertices = (Vertex *) malloc(sizeof(Vertex) * capacity);

    return new_vertex_set;
}

void freeVertexSet(VertexSet *set)
{
    // free the vertices before the set is freed
    if(set->vertices != NULL)
        free(set->vertices);
    if(set->vertices_bitMap != NULL)
        free(set->vertices_bitMap);
    free(set);
}

void addVertex(VertexSet *set, Vertex v)
{
    // TODO: Implement
    // Assume this function is only called by outsider;
    if(set->ifarray){
        if(set->vertices == NULL){
            set->vertices = (Vertex *)malloc(sizeof(Vertex)*set->capacity);
        }
        set->vertices[set->size++] = v;
    } else {
        if(set->vertices_bitMap == NULL){
            set->vertices_bitMap = (int *)malloc(sizeof(int)*(set->numNodes+1));
            memset(set->vertices_bitMap, 0, sizeof(int) * (set->numNodes+1));
        }
        set->vertices_bitMap[v]++;
        set->size++;
    }
}

void removeVertex(VertexSet *set, Vertex v)
{
    // O(n) version;
    if(set->ifarray){
        if(set->vertices == NULL){
            printf("Remove vertex from NULL array\n");
            return;
        }
        Vertex* vs = set->vertices;
        int pos = -1;
        for(int i = 0; i < set->size; ++i){
            if(vs[i] == v){
                pos = i;
                break;
            }
        }
        if(pos == -1){
            // vertice not found;
            printf("Remove vertex not found\n");
            return;
        }
        memcpy(vs+pos, vs+pos+1, (set->size - pos - 1) * sizeof(int));
        set->size--;
    } else {
        if(set->vertices_bitMap == NULL){
            printf("Remove vertex from NULL bitmap array\n");
            return;
        }
        if(set->vertices_bitMap[v] == 0){
            printf("Remove vertex not found\n");
            return;
        }
        set->vertices_bitMap[v]--;
        set->size--;
    }
    return;
}

void printBitMap(VertexSet *set){
    printf("==========================================================\n");
    for(int i = 0; i < set->numNodes; ++i){
        printf("%d\t", set->vertices_bitMap[i]);
    }
    printf("\n");
    printf("==========================================================\n");
    nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
}

void printVertices(VertexSet *set){
    printf("==========================================================\n");
    for(int i = 0; i < set->size; ++i){
        printf("%d\t", set->vertices[i]);
    }
    printf("\n");
    printf("==========================================================\n");
    nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
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

