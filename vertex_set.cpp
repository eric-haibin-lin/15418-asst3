#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include "mic.h"
#include <ctime>


void pmemset(int * start, int val, int size){
#pragma omp parallel for schedule(static, 512)
    for(int i = 0 ; i < size; i++){
        start[i] = val;
    }
}

int inclusiveScan_inplace_yiming(int * arr, int n)
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
    new_vertex_set->vertices = NULL;
    new_vertex_set->vertices_bitMap = NULL;
    if(type == SPARSE){
        new_vertex_set->vertices = (Vertex *) malloc(sizeof(Vertex) * capacity);
    } else{
        new_vertex_set->vertices_bitMap = (int *)malloc(sizeof(int) * numNodes);
        pmemset(new_vertex_set->vertices_bitMap, 0, numNodes);
    }
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
    if(set->type == SPARSE){
        if(set->vertices == NULL){
            set->vertices = (Vertex *)malloc(sizeof(Vertex)*set->capacity);
        }
        set->vertices[set->size++] = v;
    } else {
        if(set->vertices_bitMap == NULL){
            set->vertices_bitMap = (Vertex *)malloc(sizeof(int)*set->capacity);
        }
        set->vertices_bitMap[v] = 1;
        set->size++;
        return;
    }
}

void removeVertex(VertexSet *set, Vertex v)
{
    // O(n) version;
    if(set->type == SPARSE){
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
    printf("BITMAP:%d==========================================================\n", set->numNodes);
    for(int i = 0; i < set->numNodes; ++i){
        printf("%d\t", set->vertices_bitMap[i]);
    }
    printf("\n");
    printf("==========================================================\n");
    nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
}

void printVertices(VertexSet *set){
    printf("ARRAY:%d==========================================================\n", set->size);
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

