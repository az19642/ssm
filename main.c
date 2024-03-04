#include <stdio.h>
#include <ctype.h>
#include <string.h>
#include <stdlib.h>
#include "streets.h"

// use for reading from file and stdin
#define BUFSIZE 32768
char buffer[BUFSIZE];

#define RET_OK(expr, expected, label) do { \
    if ((expr) != (expected)) goto label; \
} while(0)

static bool
load_int_array(int size, int arr[size], FILE * f)
{
    int ret;
    for (int i = 0; i < size; i++) {
        if (fscanf(f, "%d", arr+i) != 1) {
            return false;
        }
    }
    ret = fscanf(f, "\n");
    return ret == 0;
}

static void
remove_newline(char * string)
{
    char * newline = strchr(string, '\n');
    if (newline) {
        *newline = '\0';
    }
}

static struct ssmap * 
load_map(const char * filename)
{
    FILE * f = fopen(filename, "r");
    struct ssmap * map = NULL;
    int nr_nodes, nr_ways;

    if (f == NULL) {
        fprintf(stderr, "error: could not open %s\n", filename);
        return NULL;
    }

    if (fgets(buffer, BUFSIZE, f) == NULL) {
        goto done;
    }

    RET_OK(strcmp(buffer, "Simple Street Map\n"), 0, invalid);
    RET_OK(fscanf(f, "%d ways\n", &nr_ways), 1, invalid);
    RET_OK(fscanf(f, "%d nodes\n", &nr_nodes), 1, invalid);

    map = ssmap_create(nr_nodes, nr_ways);
    if (map == NULL) {
        fprintf(stderr, "error: ssmap_create(%d ,%d) failed\n", nr_nodes, nr_ways);
        goto done;
    }

    for (int i = 0; i < nr_ways; i++) {
        int id, num_nodes;
        float maxspeed;
        char which_way[8];

        /* note: we are intentionally not loading the OSM id */
        RET_OK(fscanf(f, "way %d %*d ", &id), 1, invalid);
        RET_OK(fgets(buffer, BUFSIZE, f), buffer, invalid);
        RET_OK(fscanf(f, " %f %7s %d\n", &maxspeed, which_way, &num_nodes), 3, invalid);

        remove_newline(buffer);
        bool oneway = strcmp(which_way, "oneway") == 0;
        struct way * way = NULL;

        if (num_nodes > 0) {
            int node_ids[num_nodes];
            RET_OK(load_int_array(num_nodes, node_ids, f), true, invalid);
            way = ssmap_add_way(map, id, buffer, maxspeed, oneway, num_nodes, node_ids);
        }
        else {
            fprintf(stderr, "error: way %d has no associated node(s)\n", id);
            goto cleanup;
        }

        if (way == NULL) {
            fprintf(stderr, "error: ssmap_add_way(%d) failed\n", id);
            goto cleanup;
        }
    }

    for (int i = 0; i < nr_nodes; i++) {
        int id, num_ways;
        double lat, lon;

        /* note: we are intentionally not loading the OSM id */
        RET_OK(fscanf(f, "node %d %*d %lf %lf %d\n", &id, &lat, &lon, &num_ways), 4, invalid);
        struct node * node = NULL;
        
        if (num_ways > 0) {
            int way_ids[num_ways];
            RET_OK(load_int_array(num_ways, way_ids, f), true, invalid);
            node = ssmap_add_node(map, id, lat, lon, num_ways, way_ids);
        }
        else {
            fprintf(stderr, "error: node %d has no associated way(s)\n", id);
            goto cleanup;
        }
         
        if (node == NULL) {
            fprintf(stderr, "error: ssmap_add_node(%d) failed\n", id);
            goto cleanup;
        }
    }

    // custom initialization after all nodes and ways have been added
    if (!ssmap_initialize(map)) {
        fprintf(stderr, "error: ssmap_initialize() failed\n");
        goto cleanup;
    }

    printf("%s successfully loaded. %d nodes, %d ways.\n", filename, nr_nodes, nr_ways);
    goto done;
invalid:
    fprintf(stderr, "error: %s has invalid file format\n", filename);
cleanup:
    if (map != NULL) {
        ssmap_destroy(map);
        map = NULL;
    }
done:
    fclose(f);
    return map;
}

static bool
get_integer_argument(char * line, int * iptr)
{
    char * endptr;
    remove_newline(line);

    *iptr = strtol(line, &endptr, 10);
    if (endptr && *endptr != '\0') {
        printf("error: '%s' is not an integer.\n", line);
        return false;
    }

    return true;
}

static void
handle_print(char * line, struct ssmap * map, void (* print)(const struct ssmap *, int))
{
    char * first = strtok_r(line, " \t\r\n\v\f", &line);
    char * second = strtok_r(line, " \t\r\n\v\f", &line);
    int id;

    if (first == NULL || second != NULL) {
        printf("error: invalid number of arguments.\n");
    }
    else if (get_integer_argument(first, &id)) {
        print(map, id);
        return;
    }

    printf("usage: node|way id\n");
}

static void
handle_find(char * line, struct ssmap * map)
{
    char * command = strtok_r(line, " \t\r\n\v\f", &line);
    char * first = strtok_r(line, " \t\r\n\v\f", &line);
    char * second = strtok_r(line, " \t\r\n\v\f", &line);
    char * third = strtok_r(line, " \t\r\n\v\f", &line);

    if (command == NULL) {
        /* fall through */
    }
    else if (strcmp(command, "node") == 0) {
        if (first == NULL || third != NULL) {
            printf("error: invalid number of arguments.\n");
        }
        else {
            ssmap_find_node_by_names(map, first, second);
            return;
        }
    }
    else if (strcmp(command, "way") == 0) {
        if (first == NULL || second != NULL) {
            printf("error: invalid number of arguments.\n");
        }
        else {
            ssmap_find_way_by_name(map, first);
            return;
        }    
    }
    else {
        printf("error: first argument must be either node or way.\n");
    }

    printf("usage: find way keyword | find node keyword [keyword]\n");
}

static bool
handle_path_time(char * line, struct ssmap * map, 
    double (* time_function)(const struct ssmap *, int, int *))
{
    int capacity = 1;
    int n = 0;

    // use number of space characters to determine approximate array size
    for (int i = 0; line[i] != '\0'; i++) {
        if (isspace((int)line[i])) {
            capacity++;
        }
    }

    int node_ids[capacity];
    while(true) {
        char * token = strtok_r(line, " \t\r\n\v\f", &line);

        if (token == NULL)
            break;

        if (!get_integer_argument(token, &node_ids[n++])) {
            return false;
        }
    }

    if (n < 2) {
        printf("error: must specify at least two nodes.\n");
        return false;
    }

    double result = time_function(map, n, node_ids);
    if (result >= 0.) {
        printf("%.4f minutes\n", result);
    }
    
    return true;
}

static bool
handle_path_create(char * line, struct ssmap * map)
{
    char * start = strtok_r(line, " \t\r\n\v\f", &line);
    char * finish = strtok_r(line, " \t\r\n\v\f", &line);
    char * third = strtok_r(line, " \t\r\n\v\f", &line);
    int start_id, end_id;

    if (start == NULL || finish == NULL) {
        printf("error: must specify start node and finish node.\n");
        return false;
    }
    else if (third != NULL) {
        printf("error: too many arguments.\n");
        return false;
    }

    if (!get_integer_argument(start, &start_id) || 
        !get_integer_argument(finish, &end_id)) {
        return false;
    }

    ssmap_path_create(map, start_id, end_id);
    return true;
}

static void
handle_path(char * line, struct ssmap * map)
{
    char * command = strtok_r(line, " \t\r\n\v\f", &line);

    if (command == NULL) {
        /* fall through */
    }
    else if (strcmp(command, "time") == 0) {
        if (handle_path_time(line, map, ssmap_path_travel_time))
            return;
    }
    else if (strcmp(command, "create") == 0) {
        if (handle_path_create(line, map))
            return;
    }
    else {
        printf("error: first argument must be either time or create.\n");
    }

    printf("usage: path create start finish | path time node1 node2 [nodes...]\n");
}

int 
main(int argc, const char * argv[])
{
    if (argc != 2) {
        fprintf(stderr, "usage: %s FILE\n", argv[0]);
        return 0;
    }

    struct ssmap * map = load_map(argv[1]);
    if (map == NULL) {     
        return 1;
    }

    while(true) {
        printf(">> ");
        fflush(stdout);
        char * ptr = fgets(buffer, BUFSIZE, stdin);
        if (ptr == NULL) {
            break;
        }

        char * command = strtok_r(buffer, " \t\r\n\v\f", &ptr);

        if (command == NULL) {
            /* fall through */
        }
        else if (strcmp(command, "quit") == 0) {
            break;
        }
        else if (strcmp(command, "node") == 0) {
            handle_print(ptr, map, ssmap_print_node);
        }
        else if (strcmp(command, "way") == 0) {
            handle_print(ptr, map, ssmap_print_way);
        }
        else if (strcmp(command, "find") == 0) {
            handle_find(ptr, map);
        }
        else if (strcmp(command, "path") == 0) {
            handle_path(ptr, map);
        }
        else {
            printf("error: unknown command %s. Available commands are:\n"
                   "\tnode, way, find, path, quit\n", command);
        }
    }
    
    ssmap_destroy(map);
    return 0;
}