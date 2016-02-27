#include "vertex_set.h"

#include <stdlib.h>
#include <string.h>
#include <cassert>
#include <stdio.h>
#include "mic.h"
#include <ctime>

using namespace std;

#define SCAN_BLOCK_SIZE 512

void pmemset(int * start, int val, int size){
#pragma omp parallel for schedule(static, 512)
    for(int i = 0 ; i < size; i++){
        start[i] = val;
    }
}

int denseToSparse(int* vertices, int* buffer, int n)
{
    int *base_counter;
#pragma omp parallel shared(base_counter)
{
    int total_num_blocks = omp_get_num_threads();
    // buffer to store total block counts 
    int block_size = (n + total_num_blocks - 1) / total_num_blocks;
    int thread_id = omp_get_thread_num();
#pragma omp single 
    {
       base_counter = new int[total_num_blocks];
       base_counter[0] = 0;
    }
    //for each block, do sequential scan
    int base = thread_id * block_size;
    // sequential scan 
    for (int j = 1; j < block_size && base + j < n; j++) {
        vertices[base + j] += vertices[base + j - 1];
    }

    // make sure everyone finishes
#pragma omp barrier
#pragma omp single
    // collect the base, linearly
    for (int i = 1; i < total_num_blocks; i++) {
        base_counter[i] = base_counter[i - 1] + vertices[i * block_size - 1];
    }

    // update vertices with the correct prefix sum
    const int base_count = base_counter[thread_id];
    for (int j = 0; j < block_size && base + j < n; j++) {
        vertices[base + j] += base_count;
    }
    //prefix sum done
#pragma omp barrier
#pragma omp for schedule(static, 512)
    for (int i = 1; i < n; i++) {
        if (vertices[i] != vertices[i - 1]) {
            buffer[vertices[i - 1]] = i;
        }
    }
}
    if (vertices[0] == 1) {
        buffer[0] = 0;
    }
    delete base_counter;
    return vertices[n - 1];
}


MemoryManager::~MemoryManager()
{
    if(bufferArray)
        delete[] bufferArray;

    while (!setList.empty()) {
        delete setList.front();
        setList.pop_front();
    }
}

MemoryManager memManager;

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
    if (memManager.bufferArray == NULL) {
        memManager.bufferArray = new int[numNodes];
    }
    VertexSet *new_vertex_set = NULL;
    if(memManager.setList.empty()){
        new_vertex_set = (VertexSet *)malloc(sizeof(VertexSet));
        new_vertex_set->vertices = (int *)malloc(sizeof(int) * numNodes);
    } else {
        new_vertex_set = memManager.setList.front();
        memManager.setList.pop_front();
    }
    pmemset(new_vertex_set->vertices, 0, numNodes);
    new_vertex_set->size = 0;
    new_vertex_set->numNodes = numNodes;
    new_vertex_set->type = type;
    new_vertex_set->capacity = capacity;
    return new_vertex_set;
}

void freeVertexSet(VertexSet *set)
{
    // free the vertices before the set is freed
    memManager.setList.push_front(set);
}

void addVertex(VertexSet *set, Vertex v)
{
    // Assume this function is only called by outsider;
    if(set->type == SPARSE){
        set->vertices[set->size++] = v;
    } else {
        set->vertices[v] = 1;
        set->size++;
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
        if(set->vertices == NULL){
            printf("Remove vertex from NULL bitmap array\n");
            return;
        }
        if(set->vertices[v] == 0){
            printf("Remove vertex not found\n");
            return;
        }
        set->vertices[v]--;
        set->size--;
    }
    return;
}

void printVertexSet(VertexSet *set){
    if (set->type == DENSE) {
        printf("=====================BITMAP %d==================================\n", set->numNodes);
        for(int i = 0; i < set->numNodes; ++i){
            printf("%d ", set->vertices[i]);
        }
    } else {
        printf("======================ARRAY %d==================================\n", set->size);
        for(int i = 0; i < set->size; ++i){
            printf("%d ", set->vertices[i]);
        }
    }
    printf("\n");
    printf("================================================================\n");
    nanosleep((const struct timespec[]){{0, 500000000L}}, NULL);
}

void transform(VertexSet * set)
{
    if (set->type == DENSE) {
        set->size = denseToSparse(set->vertices, memManager.bufferArray, set->numNodes);
        swap(set->vertices, memManager.bufferArray);
        set->type = SPARSE;
    }
}

// helper function for union(input set, output set)
inline void vertexUnionHelper(VertexSet *input, VertexSet *output) {
    Vertex* input_vertices = input->vertices;
    Vertex* output_vertices = output->vertices;
    int size = 0;
    int count;
    if (input->type == DENSE) {
        count = input->numNodes;
#pragma omp parallel for schedule(static, 512) reduction(+: size)
        for (int i = 0; i < count; i++) {
            if (output_vertices[i] == 0 && input_vertices[i] == 1) {
                output_vertices[i] = 1;
                size += 1;
            }
        }
    } else {
        count = input->size;
#pragma omp parallel for schedule(static, 512) reduction(+: size)
        for (int i = 0; i < count; i++) {
            int vertex = input_vertices[i];
            if (output_vertices[vertex] == 0) {
                output_vertices[vertex] = 1;
                size += 1;
            }
        }
    }
    output->size += size;
}

/**
 * Returns the union of sets u and v. Destroys u and v. 
 */
VertexSet* vertexUnion(VertexSet *u, VertexSet* v)
{
    // assume u is in type DENSE
    // DANGER: change u in place! 
    assert(u->type == DENSE);
    vertexUnionHelper(v, u);
    freeVertexSet(v);
    return u;
}

