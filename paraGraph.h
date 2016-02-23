#ifndef __PARAGRAPH_H__
#define __PARAGRAPH_H__

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <omp.h>

#include "vertex_set.h"
#include "graph.h"
#include "mic.h"

#define V_RATIO 50
#define E_RATIO 25

inline int pCountOutEdgeNum(Graph g, VertexSet *u)
{
    int count = 0;
    if(u->type == SPARSE){
        Vertex* vs = u->vertices;
        int size = u->size;
#pragma omp parallel for schedule(dynamic, 512) reduction(+: count)
        for(int i = 0 ; i < size; i++){
            Vertex s = vs[i];    
            if(s < 0) continue;
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
            count += end - start;
        }
    } else {
        int * u_bitmap = u->vertices_bitMap;
        int numNodes = u->numNodes;
#pragma omp parallel for schedule(dynamic, 512) reduction(+: count)
        for(int i = 0 ; i < numNodes; i++){
            if(u_bitmap[i] == 0) continue;
            Vertex s = i;    
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
            count += end - start;
        }

    }
    return count;
}


/*
 * edgeMap --
 * 
 * Students will implement this function.
 * 
 * The input argument f is a class with the following methods defined:
 *   bool update(Vertex src, Vertex dst)
 *   bool cond(Vertex v)
 *
 * See apps/bfs.cpp for an example of such a class definition.
 * 
 * When the argument removeDuplicates is false, the implementation of
 * edgeMap need not remove duplicate vertices from the VertexSet it
 * creates when iterating over edges.  This is a performance
 * optimization when the application knows (and can tell ParaGraph)
 * that f.update() guarantees that duplicate vertices cannot appear in
 * the output vertex set.
 * 
 * Further notes: the implementation of edgeMap is templated on the
 * type of this object, which allows for higher performance code
 * generation as these methods will be inlined.
 */

    template <class F>
static VertexSet *edgeMap(Graph g, VertexSet *u, F &f,
        bool removeDuplicates=true)
{
    VertexSet * results = newVertexSet(DENSE, u->numNodes, u->numNodes);
    int * re_bitmap = results->vertices_bitMap;
    int size = 0;
    if(u->type == SPARSE){
        Vertex * vs = u->vertices;
#pragma omp parallel for schedule(dynamic, 512) reduction(+ : size)
        for(int i = 0 ; i < u->size; ++i){
            Vertex s = vs[i];    
            if(s < 0) continue;
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
            for(const Vertex* v=start; v<end; v++){
                Vertex vn = *v;
                if(f.cond(vn) && f.update(s, vn) && re_bitmap[vn] == 0){
                    re_bitmap[vn] = 1;
                    size++;
                }
            }
        }
    } else {
        int * u_bitmap = u->vertices_bitMap;
#pragma omp parallel for schedule(dynamic, 512) reduction(+:size)  
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;
            if(!f.cond(vn)) continue;
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                if(u_bitmap[s] > 0  && f.update(s, vn) && re_bitmap[vn] == 0){
                    re_bitmap[vn] = 1;
                    size++;
                }
            }
        }
    }

    results->size = size;
    // already remove duplicate;
    int outgoingEdge = pCountOutEdgeNum(g, results);
    //int outgoingEdge = 0;
    if( size * V_RATIO < u->numNodes || outgoingEdge * E_RATIO < u->numNodes){
        int * temp_bitMap = results->vertices_bitMap; //Assume it is not NULL
        inclusiveScan_inplace_yiming(temp_bitMap, results->numNodes);
        results->type = SPARSE;
        results->vertices = (Vertex *)malloc(results->size * sizeof(Vertex));
        Vertex *revs = results->vertices;
        if(temp_bitMap[0]!=0){
            revs[0] = 0;
        }
#pragma omp parallel for schedule(dynamic, 512)
        for(int i = 1; i < u->numNodes; ++i){
            if( temp_bitMap[i] != temp_bitMap[i-1]){
                revs[temp_bitMap[i-1]] = i; 
            }
        }

        free(results->vertices_bitMap);
        results->vertices_bitMap = NULL;
    }
    return results;
}



/*
 * vertexMap -- 
 * 
 * Students will implement this function.
 *
 * The input argument f is a class with the following methods defined:
 *   bool operator()(Vertex v)
 *
 * See apps/kBFS.cpp for an example implementation.
 * 
 * Note that you'll call the function on a vertex as follows:
 *    Vertex v;
 *    bool result = f(v)
 *
 * If returnSet is false, then the implementation of vertexMap should
 * return NULL (it need not build and create a vertex set)
 */
    template <class F>
static VertexSet *vertexMap(VertexSet *u, F &f, bool returnSet=true)
{
    int size = 0;
    if(u->type == SPARSE){
        VertexSet * results = NULL;
        Vertex * start = u->vertices; 
        Vertex * re_vertices = NULL;
        if(returnSet){
            results = newVertexSet(SPARSE, u->size, u->numNodes);
            re_vertices = results->vertices;
        }
#pragma omp parallel for schedule(dynamic, 512) reduction(+: size)                                                        
        for (int i = 0; i < u->size; i++) {                                                      
            if(f(start[i]) && returnSet) {
                re_vertices[i] = start[i];
                size++;
                continue;
            }
            if(returnSet) {
                re_vertices[i] = -1;
            }
        }
        if(returnSet){
            results->size = size;
            return results;
        } else {
            return NULL;
        }
    } else {
        VertexSet * results = NULL;
        if(u->vertices_bitMap == NULL){
            return NULL;
        }
        int *u_bitmap = u->vertices_bitMap;
        int *re_bitmap = NULL;
        if(returnSet){
            results = newVertexSet(DENSE, u->size, u->numNodes);
            re_bitmap = results->vertices_bitMap;
        }
#pragma omp parallel for schedule(dynamic, 512) reduction(+: size)
        for(int i = 0 ; i < u->numNodes; ++i){
            if(u_bitmap > 0 && f(i) && returnSet) {
                re_bitmap[i] = 1;
                size++;
            }
        }
        if(returnSet){
            results->size = size;
            return results;
        } else {
            return NULL;
        }
    }
}

#endif /* __PARAGRAPH_H__ */
