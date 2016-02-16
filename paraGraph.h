#ifndef __PARAGRAPH_H__
#define __PARAGRAPH_H__

#include <stdbool.h>
#include <stdlib.h>
#include <string.h>

#include "vertex_set.h"
#include "graph.h"

#include "mic.h"

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
     VertexSet *edgeMap(Graph g, VertexSet *u, F &f, bool removeDuplicates=true)
{
    // TODO: Implement
    VertexSet * results = newVertexSet(u->type, u->size, u->numNodes);
    Vertex * vs = u->vertices;
    int counter = 0;
    if(!removeDuplicates){
#pragma omp parallel for
        for(int i = 0 ; i < u->size; ++i){
           // TODO make sure vertice here is the corresponding vertices in g; 
            Vertex s = vs[i];    
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
#pragma omp parallel for 
            for(Vertex* v=start; v!=end; v++){
                Vertex vn = *v;
                if(f.cond(vn) && f.update(s, vn)){
#pragma omp critical
                    {
                        results->vertices[counter] = vn;
                        counter++;
                    }
                }
            }
        }
    } else {
        bool* visited = (bool*)malloc(sizeof(bool) * u->numNodes);
#pragma omp prallel for
        for(int i = 0; i < u->size; ++i) {
            ifInU[vs[i]] = true;
        }
#pragma omp parallel for
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = vs[i];
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
#pragma omp prallel for
            for(Vertex* v = start; v!=end; v++){
                Vertex 
            }
        }
    }
    return NULL;
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
VertexSet *vertexMap(VertexSet *u, F &f, bool returnSet=true)
{
    // TODO: Implement
    VertexSet * results = NULL;
    if(returnSet){
        results = newVertexSet(u->type, u->size, u->numNodes); 
    }
    Vertex * start = u->vertices; 
    int counter = 0;
#pragma omp parallel for                                                        
    for (int i = 0; i < u->size; i++) {                                                      
        if(f(start[i]) && returnSet) {
#pragma omp critical 
            {
                results->vertices[counter] = start[i];
                counter++;
            }
        }
    }
    if(returnSet){
        results->size = counter;
        //TODO size may excceed capacity;
        return results;
    } else {
        return NULL;
    }

}

#endif /* __PARAGRAPH_H__ */
