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

#define D_RATIO 10


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

template<class F>
VertexSet *edgeMap_BotUp_MKII(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet * results){
    // The dynamic arrangement is done outside
    // BotUp = bitmap now;
    //First we have to read information from u
    int * temp_bitMap;
    if(u->ifarray){
        printf("Array VSET in BOTUP\n");
        return NULL;
    } else {
        temp_bitMap = u->vertices_bitMap;
    }

    if(removeDuplicates){
#pragma omp parallel for  
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;    
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                int bitmap_count = temp_bitMap[s];
                while(bitmap_count>0 && f.cond(vn) && f.update(s, vn)){
                    results->vertices_bitMap[vn] = 1;
                    bitmap_count--;
                }
            }
        }
    } else {
#pragma omp parallel for  
        for(int i = 0; i < u->numNodes; ++i){
            Vertex vn = i;
            const Vertex* start = incoming_begin(g, vn);
            const Vertex* end = incoming_end(g, vn);
            for(const Vertex* v=start; v<end; v++){
                Vertex s = *v;
                int bitmap_count = temp_bitMap[s];
                while(bitmap_count>0 && f.cond(vn) && f.update(s, vn)){
                    results->vertices_bitMap[vn]++;
                    bitmap_count--;
                }
            }
        }
    }
    // calculate count using ES
    int *temp = (int *)malloc(sizeof(int) * (u->numNodes + 1));
    memcpy(temp, results->vertices_bitMap, sizeof(int)*(u->numNodes + 1)); 
    temp[u->numNodes] = 0;

    inclusiveScan_inplace_yiming(temp, u->numNodes+1);
    results->size = temp[u->numNodes];
    //printf("Printing Temp:\n");
    //for(int i = 0 ; i < u->numNodes+1; ++i){
    //    printf("%d\t", temp[i]);
    //}
    //printf("\n");
    free(temp);
    //remember to free(temp)
    return results;
}


    template<class F>
VertexSet *edgeMap_TopDown_MKII(Graph g, VertexSet *u, F &f, bool removeDuplicates, VertexSet *results)
{
    // assume array->array creation;
    int counter = 0;
    Vertex *vs = u->vertices;
    bool * visited = NULL;
    if(removeDuplicates){
        // resutls has vertices not NULL;
        visited = (bool*)malloc(sizeof(bool) * u->numNodes);
        memset(visited, 0, sizeof(bool) * u->numNodes);
    }
#pragma omp parallel for
    for(int i = 0 ; i < u->size; ++i){
        Vertex s = vs[i];    
        const Vertex* start = outgoing_begin(g, s);
        const Vertex* end = outgoing_end(g, s);
        //#pragma omp parallel for 
        for(const Vertex* v=start; v<end; v++){
            Vertex vn = *v;
            if(removeDuplicates && f.cond(vn) && f.update(s, vn)){
#pragma omp critical
                {
                if(!visited[vn])
                {
                    results->vertices[counter] = vn;
                    counter++;
                    visited[vn] = true;
                }
                }
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
static VertexSet *edgeMap(Graph g, VertexSet *u, F &f,
    bool removeDuplicates=true)
{
    VertexSet * results = newVertexSet(u->type, u->numNodes, u->numNodes);
   if(u->size * D_RATIO < u->numNodes){
    //if(false){
        // make sure it is ifarray type;
        if(!u->ifarray){
            //printf("U is a BITMAP, Changing it to ARRAY...\n");
            u->ifarray = true;
            u->vertices = (Vertex *)malloc(u->size * sizeof(Vertex));
            Vertex *revs = u->vertices;
            int * temp_bitMap = u->vertices_bitMap; //Assume it is not NULL
            if(temp_bitMap == NULL){ printf("Error EdgeMap: NULL BITMAP\n"); return NULL;}

            inclusiveScan_inplace_yiming(temp_bitMap, u->numNodes+1);

            int diff;
            if(temp_bitMap[0]!=0){
                diff = temp_bitMap[0];
                for(int j = 0 ; j < diff; ++j){
                    revs[j] = 0;
                }
            }
#pragma omp parallel for
            for(int i = 1; i < u->numNodes; ++i){
                if((diff = temp_bitMap[i] - temp_bitMap[i-1]) != 0){
                    // we have diff;
                    int offset = temp_bitMap[i-1];
                    for(int j = 0; j < diff; j++){
                        revs[offset + j] = i; 
                    }
                }
            }

            // free 
            free(u->vertices_bitMap);
            u->vertices_bitMap = NULL;
        }
        //printf("Printing u from ARRAY EdgeMap\n");
        results->ifarray = true;
        results->vertices = (Vertex *)malloc(u->numNodes * sizeof(Vertex));
        edgeMap_TopDown_MKII(g,u,f,removeDuplicates, results);
        //printVertices(results);
    } else {
        if(u->ifarray){
            printf("U is an ARRAY, Changing it to BITMAP.\n");
            u->ifarray = false;
            Vertex * vs = u->vertices;
            if(vs == NULL){printf("ERROR EdgeMap: NULL VERTICES\n"); return NULL;}
            // let's assume malloc will always succeed on lateday
            u->vertices_bitMap = (int *)malloc(sizeof(int) * (u->numNodes+1));
            memset(u->vertices_bitMap, 0, sizeof(int)*(u->numNodes+1));
#pragma omp parallel for
            for(int i = 0; i < u->size; ++i){
                u->vertices_bitMap[vs[i]]++;
            }
            free(u->vertices);
            u->vertices = NULL;
        }
        //printf("Printing u from BITMAP EdgeMap\n");
        //printBitMap(u);
        results->ifarray = false;
        results->vertices_bitMap = (int *)malloc((u->numNodes+1) * sizeof(int));
        memset(results->vertices_bitMap, 0, (u->numNodes+1) * sizeof(int));
        edgeMap_BotUp_MKII(g,u,f,removeDuplicates, results);
        //printBitMap(results);
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
    // because the new set will not be larger than u;
    // Thus we donot need to change for array;
    // Maybe we donot neet to change for !array;
    // printf("Calling VMAP!\n");
    if(u->ifarray){
        //printf("Printing u, %d\n", returnSet);
        //printVertices(u);
        VertexSet * results = NULL;
        if(returnSet){
            results = newVertexSet(u->type, u->size, u->numNodes); 
            results->ifarray = true;
            results->vertices = (Vertex *)malloc(u->size * sizeof(Vertex));
        }
        Vertex * start = u->vertices; 
        int counter = 0;
#pragma omp parallel for                                                        
        for (int i = 0; i < u->size; i++) {                                                      
#pragma omp critical 
            {
                if(f(start[i]) && returnSet) {
                    results->vertices[counter] = start[i];
                    counter++;
                }
            }
        }
        if(returnSet){
            results->size = counter;
            return results;
        } else {
            return NULL;
        }
    } else {
        VertexSet * results = NULL;
        if(u->vertices_bitMap == NULL){
            //printf("VMAP error: NULL bitmap.\n");
            return NULL;
        }
        //printf("Printing U, %d:\n", returnSet);
        //printBitMap(u);
        int *u_bitmap = u->vertices_bitMap;
        int *re_bitmap = NULL;
        if(returnSet){
            results = newVertexSet(u->type, u->size, u->numNodes);
            results->ifarray = false;
            results->vertices_bitMap = (int *)malloc(sizeof(int) * (u->numNodes + 1));
            memset(results->vertices_bitMap, 0, sizeof(int)*(u->numNodes+1));
            re_bitmap = results->vertices_bitMap;
        }
#pragma omp parallel for
        for(int i = 0 ; i < u->numNodes; ++i){
            int bitmap_count = u_bitmap[i];
            while(bitmap_count > 0 && f(i)) {
                //BUG here, might need to do while;
                //printf("vmaping...for %d\n", i);
                if(returnSet){
                    re_bitmap[i]++;
                }
                bitmap_count--;
            }
        }
        if(returnSet){
            int *temp = (int *)malloc(sizeof(int) * (u->numNodes + 1));
            memcpy(temp, results->vertices_bitMap, sizeof(int)*(u->numNodes + 1)); 
            temp[u->numNodes] = 0;

            inclusiveScan_inplace_yiming(temp, u->numNodes+1);
            results->size = temp[u->numNodes];
            free(temp);
            return results;
        } else {
            return NULL;
        }
    }
}

#endif /* __PARAGRAPH_H__ */
