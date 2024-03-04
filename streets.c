#include "streets.h"

#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

struct node {
    int id;
    double lat;
    double lon;
    int num_ways;
    int *way_ids;
};

struct way {
    int id;
    char *name;
    float maxspeed;
    bool oneway;
    int num_nodes;
    int *node_ids;
};

struct ssmap {
    struct node **nodes;
    struct way **ways;
    int num_nodes;
    int num_ways;
};

struct ssmap *ssmap_create(int nr_nodes, int nr_ways) {
    if (nr_nodes == 0 || nr_ways == 0) {
        return NULL;
    }

    struct ssmap *m = malloc(sizeof(struct ssmap));
    if (m == NULL) {
        return NULL;
    }

    m->nodes = malloc(sizeof(struct node *) * nr_nodes);
    if (m->nodes == NULL) {
        free(m);
        return NULL;
    }
    m->ways = malloc(sizeof(struct way *) * nr_ways);
    if (m->ways == NULL) {
        free(m->nodes);
        free(m);
        return NULL;
    }

    m->num_nodes = nr_nodes;
    m->num_ways = nr_ways;

    return m;
}

bool ssmap_initialize(struct ssmap *m) {
    /* TODO: task 2
     * additional initialization code can be added here */
    return true;
}

void ssmap_destroy(struct ssmap *m) {
    // loop through all node ids in m
    for (int i = 0; i < m->num_nodes; i++) {
        free(m->nodes[i]->way_ids);
        free(m->nodes[i]);
    }
    free(m->nodes);

    // loop through all way ids in m
    for (int i = 0; i < m->num_ways; i++) {
        free(m->ways[i]->name);
        free(m->ways[i]->node_ids);
        free(m->ways[i]);
    }
    free(m->ways);

    free(m);
}

struct way *ssmap_add_way(struct ssmap *m, int id, const char *name,
                          float maxspeed, bool oneway, int num_nodes,
                          const int node_ids[num_nodes]) {
    struct way *w = malloc(sizeof(struct way));
    if (w == NULL) {
        return NULL;
    }
    w->id = id;
    w->maxspeed = maxspeed;
    w->oneway = oneway;
    w->num_nodes = num_nodes;

    // strlen(name) + 1 is the space needed to store name (+ 1 to include the
    // extra byte required for the null terminator)
    w->name = malloc(sizeof(char) * (strlen(name) + 1));
    if (w->name == NULL) {
        free(w);
        return NULL;
    }
    strcpy(w->name, name);

    w->node_ids = malloc(sizeof(int) * num_nodes);
    if (w->node_ids == NULL) {
        free(w->name);
        free(w);
        return NULL;
    }
    // copy node_ids into w.node_ids
    for (int i = 0; i < num_nodes; i++) {
        w->node_ids[i] = node_ids[i];
    }

    m->ways[id] = w;

    return w;
}

struct node *ssmap_add_node(struct ssmap *m, int id, double lat, double lon,
                            int num_ways, const int way_ids[num_ways]) {
    struct node *n = malloc(sizeof(struct node));
    if (n == NULL) {
        return NULL;
    }

    n->id = id;
    n->lat = lat;
    n->lon = lon;
    n->num_ways = num_ways;

    n->way_ids = malloc(sizeof(int) * num_ways);
    if (n->way_ids == NULL) {
        free(n);
        return NULL;
    }
    // copy way_ids into node_ptr->way_ids
    for (int i = 0; i < num_ways; i++) {
        n->way_ids[i] = way_ids[i];
    }

    m->nodes[id] = n;

    return n;
}

void ssmap_print_way(const struct ssmap *m, int id) {
    if (id < 0 || id >= m->num_ways) {
        printf("error: way %d does not exist.\n", id);
        return;
    }
    struct way *w = m->ways[id];
    printf("Way %d: %s\n", id, w->name);
}

void ssmap_print_node(const struct ssmap *m, int id) {
    if (id < 0 || id >= m->num_nodes) {
        printf("error: node %d does not exist.\n", id);
        return;
    }
    struct node *n = m->nodes[id];
    printf("Node %d: (%.7lf, %.7lf)\n", id, n->lat, n->lon);
}

void ssmap_find_way_by_name(const struct ssmap *m, const char *name) {
    struct way **ways = m->ways;
    for (int i = 0; i < m->num_ways; i++) {
        if (strstr(ways[i]->name, name) != NULL) {
            printf("%d ", i);
        }
    }
    printf("\n");
}

void ssmap_find_node_by_names(const struct ssmap *m, const char *name1,
                              const char *name2) {
    struct node **nodes = m->nodes;
    struct way **ways = m->ways;

    // criterion = 1 or 2 depending on if it contains ways whose names contain
    // n1 (1), both (2), or none of the names (0)

    // designed this way to easily take care of "edge cases" when using this
    // function with two names (name2 != NULL)
    int criterion[m->num_nodes];
    for (int i = 0; i < m->num_nodes; i++) {
        criterion[i] = 0;
    }
    for (int i = 0; i < m->num_nodes; i++) {
        // i is a node_id
        for (int j = 0; j < nodes[i]->num_ways; j++) {
            // j is an index for this node's array of way_ids
            int way_id = nodes[i]->way_ids[j];
            // check if name1 is a substring in ways[way_id]->name
            if (strstr(ways[way_id]->name, name1) != NULL) {
                criterion[i] = 1;
                break;
            }
        }
    }

    if (name2 == NULL) {
        for (int i = 0; i < m->num_nodes; i++) {
            if (criterion[i] == 1) {
                // node i contains ways whose names contain n1
                printf("%d ", i);
            }
        }
    } else {
        for (int i = 0; i < m->num_nodes; i++) {
            for (int j = 0; j < nodes[i]->num_ways; j++) {
                int way_id = nodes[i]->way_ids[j];
                // check if name1 is a substring in ways[way_id]->name
                if (criterion[i] == 1 &&
                    strstr(ways[way_id]->name, name2) != NULL) {
                    criterion[i] = 2;
                    // this node (with id=i) has potentially different ways
                    // that have each given name appearing in the way's
                    // name; we now proceed to check for uniqueness
                    // seperately AND find the right combination (which way
                    // should be used for a certain name, this affects validity
                    // of a node)
                    break;
                }
            }
        }
        for (int i = 0; i < m->num_nodes; i++) {
            // we can immediately skip nodes that do not contain (in terms
            // of way names) both n1 and n2
            if (criterion[i] != 2) continue;

            // brute force by checking every combination of possible way_ids for
            // this node corresponding to each name; using a double for loop
            // over node i's ways
            for (int j = 0; j < nodes[i]->num_ways; j++) {
                bool found1 = false;
                bool found2 = false;

                int way_id_j = nodes[i]->way_ids[j];

                if (strstr(ways[way_id_j]->name, name1) != NULL) {
                    found1 = true;
                } else if (strstr(ways[way_id_j]->name, name2) != NULL) {
                    found2 = true;
                }

                for (int k = 0; k < nodes[i]->num_ways; k++) {
                    // we need the way for n2 to be different from n1
                    if (j == k) continue;

                    int way_id_k = nodes[i]->way_ids[k];

                    if (found1 && strstr(ways[way_id_k]->name, name2) != NULL) {
                        found2 = true;
                        break;
                    } else if (found2 &&
                               strstr(ways[way_id_k]->name, name1) != NULL) {
                        found1 = true;
                        break;
                    }
                }
                if (found1 && found2) {
                    // found a valid combination so this node is valid, we can
                    // move onto the next node
                    printf("%d ", i);
                    break;
                }
            }
        }
    }
    printf("\n");
}

/**
 * Converts from degree to radian
 *
 * @param deg The angle in degrees.
 * @return the equivalent value in radian
 */
#define d2r(deg) ((deg)*M_PI / 180.)

/**
 * Calculates the distance between two nodes using the Haversine formula.
 *
 * @param x The first node.
 * @param y the second node.
 * @return the distance between two nodes, in kilometre.
 */
static double distance_between_nodes(const struct node *x,
                                     const struct node *y) {
    double R = 6371.;
    double lat1 = x->lat;
    double lon1 = x->lon;
    double lat2 = y->lat;
    double lon2 = y->lon;
    double dlat = d2r(lat2 - lat1);
    double dlon = d2r(lon2 - lon1);
    double a = pow(sin(dlat / 2), 2) +
               cos(d2r(lat1)) * cos(d2r(lat2)) * pow(sin(dlon / 2), 2);
    double c = 2 * atan2(sqrt(a), sqrt(1 - a));
    return R * c;
}

bool is_forward_adj(const struct node *n1, const struct node *n2,
                    const int *node_ids, int size) {
    for (int i = 0; i < size - 1; i++) {
        // if ith node id is equal to node1 id, then the (i+1)th node id
        // must be equal to node2 id for forward adjacency
        if (node_ids[i] == n1->id) {
            return (node_ids[i + 1] == n2->id);
        }
    }
    // never reached since node1 id and node2 id must be in node ids under
    // the context of this program
    return false;
}

/**
 * Helper designed this way to assist in the next part; this helper checks
 * for errors 2 to 4.
 *
 * @param m
 * @param n1
 * @param n2
 * @return travel time between n1 and n2, or an appropriate value
 * representing an error.
 */

double travel_time_helper(const struct ssmap *m, struct node *n1,
                          struct node *n2) {
    double travel_time = 0.0;

    int shared_w_ids[n1->num_ways];
    int num_shared_w_ids = 0;

    for (int i = 0; i < n1->num_ways; i++) {
        for (int j = 0; j < n2->num_ways; j++) {
            if (n1->way_ids[i] == n2->way_ids[j]) {
                shared_w_ids[num_shared_w_ids] = n1->way_ids[i];
                num_shared_w_ids++;
            }
        }
    }

    // no shared way
    if (num_shared_w_ids == 0) {
        return -2.0;
    }

    for (int i = 0; i < num_shared_w_ids; i++) {
        int shared_w_id = shared_w_ids[i];
        struct way *shared_w = m->ways[shared_w_id];

        bool forward_adj =
            is_forward_adj(n1, n2, shared_w->node_ids, shared_w->num_nodes);
        bool backward_adj =
            is_forward_adj(n2, n1, shared_w->node_ids, shared_w->num_nodes);

        if (!forward_adj && !backward_adj) {
            // no adjacency in either direction
            travel_time = -3.0;
        } else if (shared_w->oneway && !forward_adj) {
            // one way but not forward adjacent
            travel_time = -4.0;
        } else {
            return 60 * distance_between_nodes(n1, n2) /
                   m->ways[shared_w_id]->maxspeed;
        }
    }
    return travel_time;
}

double ssmap_path_travel_time(const struct ssmap *m, int size,
                              int node_ids[size]) {
    // error 1 check
    for (int i = 0; i < size; i++) {
        if (node_ids[i] < 0 || node_ids[i] >= m->num_nodes) {
            printf("error: node %d does not exist.\n", node_ids[i]);
            return -1.0;
        }
    }

    double travel_time = 0.0;

    // used for checking error 5
    int visited_node_ids[m->num_nodes];
    for (int i = 0; i < m->num_nodes; i++) {
        visited_node_ids[i] = 0;
    }

    for (int i = 0; i < size - 1; i++) {
        struct node *n1 = m->nodes[node_ids[i]];
        struct node *n2 = m->nodes[node_ids[i + 1]];

        // this function handles errors 2-4
        double pair_travel_time = travel_time_helper(m, n1, n2);

        if (pair_travel_time == -2.0) {
            printf("error: there are no roads between node %d and node %d.\n",
                   n1->id, n2->id);
            return -1.0;
        }
        if (pair_travel_time == -3.0) {
            printf("error: cannot go directly from node %d to node %d.\n",
                   n1->id, n2->id);
            return -1.0;
        }
        if (pair_travel_time == -4.0) {
            printf("error: cannot go in reverse from node %d to node %d.\n",
                   n1->id, n2->id);
            return -1.0;
        }
        // error 5 check
        if (visited_node_ids[n1->id]) {
            printf("error: node %d appeared more than once.\n", n1->id);
            return -1.0;
        }
        // consider n1 visited
        visited_node_ids[n1->id] = 1;

        // n2 has been visited
        if (visited_node_ids[n2->id]) {
            printf("error: node %d appeared more than once.\n", n2->id);
            return -1.0;
        }

        travel_time += pair_travel_time;
    }
    return travel_time;
}

/**
 * Min heap implementation of a priority queue based on item's path_time.
 *
 * @param ids Array of ids of nodes in the heap.
 * @param size Number of items/ids in the heap;
 */
typedef struct heap {
    int *ids;
    int size;
} Heap;

typedef struct heap Heap;

/**
 * @param h
 * @return true if h is empty, false otherwise.
 */
bool is_empty(const Heap *h) {
    return (h->size == 0);
    ;  // for formatting
}

/**
 * Swap two items in the heap
 *
 * @param h Heap in which we would like to perform a swap in.
 * @param index1 Index of one of the items in the swap.
 * @param index2 Index of the other item in the swap.
 */
void heap_swap(const Heap *h, int index1, int index2) {
    int item1 = h->ids[index1];
    h->ids[index1] = h->ids[index2];
    h->ids[index2] = item1;
}

/**
 * Min-heapify h rooted at root_i, given its children are min heaps. The
 * priority of k is given by times[k].
 * @param h
 * @param root_i
 * @param times
 */
void min_heapify(const Heap *h, int root_i, const double *times) {
    int root = h->ids[root_i];
    int left_i = root_i * 2 + 1;
    int right_i = root_i * 2 + 2;

    if (right_i < h->size) {
        // has both left and right child
        int left = h->ids[left_i];
        int right = h->ids[right_i];
        if (times[left] < times[root] && times[left] < times[right]) {
            // left is smallest
            heap_swap(h, root_i, left_i);
            min_heapify(h, left_i, times);
        } else if (times[right] < times[root] && times[right] < times[left]) {
            // right is smallest
            heap_swap(h, root_i, right_i);
            min_heapify(h, right_i, times);
        } else {
            // root is smallest of the three, so already a min heap
            return;
        }
    } else if (left_i < h->size) {
        // no right child
        int left = h->ids[left_i];
        if (times[left] < times[root]) {
            // left is smallest
            heap_swap(h, root_i, left_i);
            min_heapify(h, left_i, times);
        } else {
            // root is smallest of the two, so already min heap
            return;
        }
    } else {
        // no left and right, so already a min heap
        return;
    }
}

/**
 * Extract the lowest priority item in the given heap implementation of a
 * priority queue.
 *
 * @param h A non-empty heap.
 * @param times
 * @param is_in_heap
 * @return the item in h with the lowest priority.
 */
int heap_extract_min(Heap *h, const double *times, bool *is_in_heap) {
    int min_item = h->ids[0];
    h->ids[0] = h->ids[h->size - 1];
    h->size--;
    min_heapify(h, 0, times);
    is_in_heap[min_item] = false;
    return min_item;
}

/**
 * Bubble up in the given heap starting from given index, to min heapify.
 *
 * @param h The heap.
 * @param root_i The index of the starting item. Assume children are min
 * heaps.
 * @param times The priority array.
 */
void heap_insert_helper(Heap *h, int root_i, double *times) {
    int parent_i = (root_i - 1) / 2;  // implicit floor

    // bubble up
    if (times[h->ids[parent_i]] > times[h->ids[root_i]]) {
        int temp = h->ids[parent_i];

        h->ids[parent_i] = h->ids[root_i];
        h->ids[root_i] = temp;

        heap_insert_helper(h, parent_i, times);
    }
}

/**
 * Insert k into min heap h with priority times[k]
 *
 * @param h
 * @param k Key to be inserted.
 * @param times Array consisting of priorities; times[k] is the priority of k.
 * @param is_in_heap
 */
void heap_insert(Heap *h, int k, double *times, bool *is_in_heap) {
    h->ids[h->size] = k;
    h->size++;
    heap_insert_helper(h, h->size - 1, times);
    is_in_heap[k] = true;
}

/**
 * Find the neighbouring node's ids of the node with id src_id and store it
 * in a given array (neighbours)
 *
 * @param m
 * @param src_id
 * @param neighbours
 */
void update_neighbour_ids(const struct ssmap *m, int src_id, int *neighbours) {
    int num_neighbours = 0;
    struct node **nodes = m->nodes;
    struct node *src = nodes[src_id];
    // iterate over ways
    for (int i = 0; i < src->num_ways; i++) {
        struct way *w = m->ways[src->way_ids[i]];
        for (int j = 0; j < w->num_nodes; j++) {
            // for each way, iterate over the nodes
            if (nodes[w->node_ids[j]]->id == src_id) {
                if (0 <= j && j < w->num_nodes - 1) {
                    int id_right = w->node_ids[j + 1];
                    neighbours[num_neighbours++] = id_right;
                }
                if (!w->oneway && 0 < j && j <= w->num_nodes - 1) {
                    int id_left = w->node_ids[j - 1];
                    neighbours[num_neighbours++] = id_left;
                }
            }
        }
    }
    neighbours[num_neighbours] = -1;
}

void ssmap_path_create(const struct ssmap *m, int start_id, int end_id) {
    if (start_id < 0 || start_id >= m->num_nodes || end_id < 0 ||
        end_id >= m->num_nodes) {
        printf("error: node %d does not exist.\n",
               start_id < 0 || start_id >= m->num_nodes ? start_id : end_id);
        return;
    }
    // is_in_heap[node_id] = 1 if node with id node_id
    // in heap (to be defined), 0 otherwise
    bool is_in_heap[m->num_nodes];
    // times[node_id] stores the current best path time we found from
    // start_id to node_id
    double times[m->num_nodes];
    // prev[node_id] stores the previous node's id in the current best path
    // we found from start_id to node_id
    int prev_id[m->num_nodes];

    for (int i = 0; i < m->num_nodes; i++) {
        // nothing is in the heap at the start
        is_in_heap[i] = false;
        // unknown times at start so we assume arbitrarily large
        times[i] = INFINITY;
        // no paths at the start, so all prev nodes are undefined(-1)
        prev_id[i] = -1;
    }
    // we know only this time at the start, since we start here
    times[start_id] = 0.0;

    Heap *h = malloc(sizeof(Heap));
    if (h == NULL) {
        printf("error: out of memory.");
        return;
    }
    h->size = 0;
    h->ids = malloc(sizeof(int) * m->num_nodes);
    if (h->ids == NULL) {
        free(h);
        printf("error: out of memory.");
        return;
    }

    // initialize heap with starting node, since we know time[start_id]
    heap_insert(h, start_id, times, is_in_heap);

    // node has at most m->ways, and each way provides at most 2
    // neighbours on the left and right of the node
    int *neighbours = malloc(sizeof(int) * 2 * m->num_ways);
    if (neighbours == NULL) {
        free(h->ids);
        free(h);
        printf("error: out of memory.");
        return;
    }

    while (!is_empty(h)) {
        int u = heap_extract_min(h, times, is_in_heap);
        update_neighbour_ids(m, u, neighbours);

        for (int i = 0;
             neighbours[i] != -1 && i < sizeof(int) * 2 * m->num_ways; i++) {
            int v = neighbours[i];
            double u_to_v = travel_time_helper(m, m->nodes[u], m->nodes[v]);

            // negative from travel_time_helper means error, so skip this
            // neighbour, as the current path up to the neighbour is invalid

            // note that error 5 and 1 do not occur in the final path given a
            // correct
            // implementation of Dijkstra's algorithm
            if (u_to_v < 0) continue;

            double alt_time = times[u] + u_to_v;
            if (alt_time < times[v]) {
                times[v] = alt_time;
                prev_id[v] = u;
                if (!is_in_heap[v]) {
                    heap_insert(h, v, times, is_in_heap);
                }
            }
        }
    }
    free(h->ids);
    free(h);
    free(neighbours);

    if (prev_id[end_id] == -1) {
        printf("error: could not find a path from node %d to node %d.\n",
               start_id, end_id);
        return;
    }

    // reverse traversal with prev arr to get the path in the correct order
    int path[m->num_nodes];
    int path_index = m->num_nodes;
    int curr_id = end_id;

    while (curr_id != -1) {
        path[--path_index] = curr_id;
        curr_id = prev_id[curr_id];
    }
    for (; path_index < m->num_nodes; path_index++) {
        printf("%d ", path[path_index]);
    }
    printf("\n");
}