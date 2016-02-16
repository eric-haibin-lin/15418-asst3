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

inline int inclusiveScan_inplace_yiming(int * arr, int n)
{
    int  *partial, *temp;
    int num_threads, work;
    int i, mynum, last;
    if(arr == NULL) return -1;
#pragma omp parallel default(none) private(i, mynum, last) shared(arr, partial, temp, num_threads, work, n)
    {
#pragma omp single
        {
            num_threads = omp_get_num_threads();
            if(!(partial = (int *) malloc (sizeof (int) * num_threads))) exit(-1);
            if(!(temp = (int *) malloc (sizeof (int) * num_threads))) exit(-1);
            work = n / num_threads + 1; /*sets length of sub-arrays*/
        }
        mynum = omp_get_thread_num();
        /*calculate prefix-sum for each subarray*/
        for(i = work * mynum + 1; i < work * mynum + work && i < n; i++)
            arr[i] += arr[i - 1];
        partial[mynum] = arr[i - 1];
#pragma omp barrier
        /*calculate prefix sum for the array that was made from last elements of each of the previous sub-arrays*/
        for(i = 1; i < num_threads; i <<= 1) {
            if(mynum >= i)
                temp[mynum] = partial[mynum] + partial[mynum - i];
#pragma omp barrier
#pragma omp single
            memcpy(partial + 1, temp + 1, sizeof(int) * (num_threads - 1));
        }
        /*update original array*/
        for(i = work * mynum; i < (last = work * mynum + work < n ? work * mynum + work : n); i++)
            arr[i] += partial[mynum] - arr[last - 1];
    }
    //for(int i = 0; i < n; ++i){
    //    printf("%d\t", arr[i]);
    //}
    //printf("\n");
    free(partial);
    free(temp);
    return 0;
}



template <class F>
VertexSet *edgeMap_ES(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet *results){
    if(results == NULL) return NULL;
    Vertex * vs = u->vertices;
    // let's assume malloc will always succeed on lateday
    int *temp = (int *)malloc(sizeof(int) * (u->numNodes + 1));
    memset(temp, 0, sizeof(int)*(u->numNodes+1));
    if(removeDuplicates){
#pragma omp parallel for  
        for(int i = 0; i < u->size; ++i){
            Vertex s = vs[i];    
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
#pragma omp parallel for 
            for(const Vertex* v=start; v<end; v++){
                Vertex vn = *v;
                if(f.cond(vn) && f.update(s, vn)){
                    temp[vn+1] = 1;
                }
            }
        }
    } else {
#pragma omp parallel for
        for(int i = 0; i < u->size; ++i){
            Vertex s = vs[i];    
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
#pragma omp parallel for 
            for(const Vertex* v=start; v<end; v++){
                Vertex vn = *v;
                if(f.cond(vn) && f.update(s, vn)){
                    temp[vn+1]++;
                }
            } 
        }
    }
    // inplace inclusive scan;
    inclusiveScan_inplace_yiming(temp, u->numNodes + 1);

    // now we have 0,1,2,3...
    results->size = temp[u->numNodes];
    Vertex *revs = results->vertices;
    int diff;
#pragma omp parallel for
    for(int i = 0; i < u->numNodes; ++i){
        if((diff = temp[i+1] - temp[i]) != 0){
            // we have diff;
            int offset = temp[i];
            for(int j = 0; j < diff; j++){
                revs[offset + j] = i; 
            }
        }
    }
    free(temp);
    return results;
}

template<class F>
VertexSet *edgeMap_BotUp(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet * results){
    if(results == NULL) return NULL;
    Vertex * vs = u->vertices;
    // let's assume malloc will always succeed on lateday
    bool *temp = (bool *)malloc(sizeof(bool) * (u->numNodes));
    memset(temp, 0, sizeof(bool)*(u->numNodes));
    int count = 0;

#pragma omp parallel for
    for(int i = 0; i < u->size; ++i){
        temp[vs[i]] = true;
    }
    if(removeDuplicates){
#pragma omp parallel for  
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;    
            if(!f.cond(vn)) continue;
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
            bool ifvisited = false;
#pragma omp parallel for 
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                if(temp[s] && f.update(s, vn) && !ifvisited){
#pragma omp critical
                    {
                        results->vertices[count++] = vn; 
                    }
                }
            }
        }
    } else {
#pragma omp parallel for  
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;    
            if(!f.cond(vn)) continue;
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
#pragma omp parallel for 
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                if(temp[s] && f.update(s, vn)){
#pragma omp critical
                    {
                        results->vertices[count++] = vn; 
                    }
                }
            }
        }
    
    }
    free(temp);
    results->size = count;
    return results;
}

template<class F>
VertexSet *edgeMap_TopDown(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet *results)
{
    int counter = 0;
    Vertex *vs = u->vertices;
        bool * visited = NULL;
        if(removeDuplicates){
            visited = (bool*)malloc(sizeof(bool) * u->numNodes);
            memset(visited, 0, sizeof(bool) * u->numNodes);
        }
#pragma omp parallel for
        for(int i = 0 ; i < u->size; ++i){
            // TODO make sure vertice here is the corresponding vertices in g; 
            Vertex s = vs[i];    
            const Vertex* start = outgoing_begin(g, s);
            const Vertex* end = outgoing_end(g, s);
#pragma omp parallel for 
            for(const Vertex* v=start; v<end; v++){
                Vertex vn = *v;
                if(removeDuplicates && f.cond(vn) && f.update(s, vn) && !visited[vn]){
#pragma omp critical
                    {
                        results->vertices[counter] = vn;
                        counter++;
                    }
                    visited[vn] = true;
                }
                if(!removeDuplicates && f.cond(vn) && f.update(s, vn)){
#pragma omp critical
                    {
                        results->vertices[counter] = vn;
                        counter++;
                    }
                }
            }
        }
        if(removeDuplicates){
            free(visited);
        }
        results->size = counter;
    return results;
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
     VertexSet *edgeMap(Graph g, VertexSet *u, F &f, bool removeDuplicates=true)
{
    // TODO: Implement
    VertexSet * results = newVertexSet(u->type, u->numNodes, u->numNodes);
    int dynamic_ratio = 5;
    if(g->num_edges / g->num_nodes < dynamic_ratio){
        edgeMap_TopDown(g, u, f, removeDuplicates, results);
    } else {
        // TODO arge scale
        edgeMap_ES(g, u, f, removeDuplicates, results);
        //edgeMap_BotUp(g,u,f,removeDuplicates, results);
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
VertexSet *vertexMap(VertexSet *u, F &f, bool returnSet=true)
{
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
