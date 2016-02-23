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

#define D_RATIO 1000


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
    free(partial);
    free(temp);
    return 0;
}

inline void pmemset(int * start, int val, int size){
#pragma omp parallel for schedule(static, 512)
    for(int i = 0 ; i < size; i++){
        start[i] = val;
    }
}

template<class F>
VertexSet *edgeMap_BotUp_MKII(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet * results){
    // The dynamic arrangement is done outside
    // BotUp = bitmap now;
    //First we have to read information from u
    int * temp_bitMap;
    temp_bitMap = u->vertices_bitMap;
    int size = 0;
    if(removeDuplicates){
#pragma omp parallel for schedule(dynamic, 512) reduction(+:size) 
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;    
            if(!f.cond(vn)) continue;
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                if(temp_bitMap[s]>0 && f.update(s, vn)){
                    results->vertices_bitMap[vn] = 1;
                    size++;
                }
            }
        }
    } else {
#pragma omp parallel for schedule(dynamic, 512) reduction(+:size)  
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;
            if(!f.cond(vn)) continue;
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                if(temp_bitMap[s] > 0  && f.update(s, vn)){
                    results->vertices_bitMap[vn]++;
                    size++;
                }
            }
        }
    }
    
    results->size = size;
    return results;
}


    template<class F>
VertexSet *edgeMap_TopDown_MKII(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet *results)
{
    // assume array->array creation;
    Vertex *vs = u->vertices;
    int * visited = (int*)malloc(sizeof(int) * (u->numNodes+1));
    pmemset(visited, 0, (u->numNodes+1));
#pragma omp parallel for schedule(dynamic, 512)
    for(int i = 0 ; i < u->size; ++i){
        Vertex s = vs[i];    
        const Vertex* start = outgoing_begin(g, s);
        const Vertex* end = outgoing_end(g, s);
        for(const Vertex* v=start; v<end; v++){
            Vertex vn = *v;
            if(removeDuplicates && f.cond(vn) && f.update(s, vn)){
                visited[vn] = 1;
            }
            if(!removeDuplicates && f.cond(vn) && f.update(s, vn)){
                visited[vn]++;
            }
        }
    }
    inclusiveScan_inplace_yiming(visited, u->numNodes+1);
    results->size = visited[u->numNodes];
    Vertex * revs = results->vertices;
    int diff;
    if(visited[0]!=0){
        revs[0] = 0;
    }
#pragma omp parallel for schedule(dynamic, 512)
    for(int i = 1; i < u->numNodes; ++i){
        if((diff = visited[i] - visited[i-1]) != 0){
            int offset = visited[i-1];
                revs[offset] = i; 
        }
    }
    free(visited);
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
static VertexSet *edgeMap(Graph g, VertexSet *u, F &f,
        bool removeDuplicates=true)
{
    VertexSet * results = newVertexSet(u->type, u->numNodes, u->numNodes);
    if( u->size * D_RATIO < u->numNodes){
        // make sure it is ifarray type;
        if(!u->ifarray){
            u->ifarray = true;
            u->vertices = (Vertex *)malloc(u->size * sizeof(Vertex));
            Vertex *revs = u->vertices;
            int * temp_bitMap = u->vertices_bitMap; //Assume it is not NULL
            inclusiveScan_inplace_yiming(temp_bitMap, u->numNodes+1);
            int diff;
            if(temp_bitMap[0]!=0){
                diff = temp_bitMap[0];
                for(int j = 0 ; j < diff; ++j){
                    revs[j] = 0;
                }
            }
#pragma omp parallel for schedule(dynamic, 512)
            for(int i = 1; i < u->numNodes; ++i){
                if((diff = temp_bitMap[i] - temp_bitMap[i-1]) != 0){
                    revs[temp_bitMap[i-1]] = i; 
                }
            }

            free(u->vertices_bitMap);
            u->vertices_bitMap = NULL;
        }
        results->ifarray = true;
        results->vertices = (Vertex *)malloc(u->numNodes * sizeof(Vertex));
        edgeMap_TopDown_MKII(g,u,f,removeDuplicates, results);
    } else {
        if(u->ifarray){
            u->ifarray = false;
            Vertex * vs = u->vertices;
            u->vertices_bitMap = (int *)malloc(sizeof(int) * (u->numNodes+1));
            pmemset(u->vertices_bitMap, 0, (u->numNodes+1));
#pragma omp parallel for schedule(dynamic, 512)
            for(int i = 0; i < u->size; ++i){
                u->vertices_bitMap[vs[i]]=1;
            }
            free(u->vertices);
            u->vertices = NULL;
        }
        results->ifarray = false;
        results->vertices_bitMap = (int *)malloc((u->numNodes+1) * sizeof(int));
        pmemset(results->vertices_bitMap, 0, (u->numNodes+1));
        edgeMap_BotUp_MKII(g,u,f,removeDuplicates, results);
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
    if(u->ifarray){
        VertexSet * results = NULL;
        int* re_bitmap = NULL;
        if(returnSet){
            results = newVertexSet(u->type, u->size, u->numNodes); 
            results->ifarray = false;
            //results->vertices = (Vertex *)malloc(u->size * sizeof(Vertex));
            results->vertices_bitMap = (int *)malloc(sizeof(int) * (u->numNodes + 1));
            pmemset(results->vertices_bitMap, 0, (u->numNodes+1));
            re_bitmap = results->vertices_bitMap;
        }
        Vertex * start = u->vertices; 
#pragma omp parallel for schedule(dynamic, 512) reduction(+: size)                                                        
        for (int i = 0; i < u->size; i++) {                                                      
            if(f(start[i]) && returnSet) {
                re_bitmap[start[i]] = 1;
                size++;
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
            results = newVertexSet(u->type, u->size, u->numNodes);
            results->ifarray = false;
            results->vertices_bitMap = (int *)malloc(sizeof(int) * (u->numNodes + 1));
            pmemset(results->vertices_bitMap, 0, (u->numNodes+1));
            re_bitmap = results->vertices_bitMap;
        }
#pragma omp parallel for schedule(dynamic, 512) reduction(+: size)
        for(int i = 0 ; i < u->numNodes; ++i){
            if(u_bitmap > 0 && f(i) && returnSet) {
                re_bitmap[i]++;
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
