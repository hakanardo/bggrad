#include <stdio.h>
#include <math.h>
#include "graph.h"

typedef Graph<int,int,int> GraphType;
        
extern "C" {
    #include "capi.h"

    void * new_graph(int nodes, int edges) {
        GraphType *g = new GraphType(nodes, edges);
        g->add_node(nodes);
        g->maxflow(); // To initialize structures used by mark_node
        return g;
    }
    void add_capacity(void *p, int nodes, int *source, int *sink) {
        GraphType *g = (GraphType *) p;
        for (int i=0; i<nodes; i++) {
            int old = g->get_trcap(i);
            g->add_tweights(i, source[i], sink[i]);
            int cap = g->get_trcap(i);
            if (!( (old<0 && cap<0) || (old>0 && cap>0) )) {
                g->mark_node(i);
            }
        }
    }
    void add_edges(void *p, int edges, int *head, int *tail, int *cap) {
        GraphType *g = (GraphType *) p;
        for (int i=0; i<edges; i++) {
            g->add_edge(head[i], tail[i], cap[i], cap[i]);
        }
    }
    void add_edge(void *p, int head, int tail, int cap) {
        GraphType *g = (GraphType *) p; 
        g->add_edge(head, tail, cap, cap);
    }
    int maxflow(void *p, int reuse) {
        GraphType *g = (GraphType *) p; 
        return g->maxflow(reuse);
    }
    void get_segments(void *p, int nodes, int *segs) {
        GraphType *g = (GraphType *) p; 
        for (int i=0; i<nodes; i++) {
            segs[i] = (g->what_segment(i) == GraphType::SOURCE);
        }
    }
    void free_graph(void *p) {
        GraphType *g = (GraphType *) p;
        delete g;
    }
    
    void set_capacity_from_image(void *ptr, int nodes, int *old_source, int *old_sink,
                                 unsigned char *image, float pri, float scale) {
        GraphType *g = (GraphType *) ptr; 
        for (int i=0; i<nodes; i++) {
            float p = image[i];
            float pp =  pri * p / (pri * p + (1.0 - pri) * (255.0 - p));
            int source = (int) (-logf(1.0 - pp + 0.000001) * scale);
            int sink = (int) (-logf(pp + 0.000001) * scale);
            int old = g->get_trcap(i);
            g->add_tweights(i, source - old_source[i], sink - old_sink[i]);
            int cap = g->get_trcap(i);
            if (!( (old<0 && cap<0) || (old>0 && cap>0) )) {
                g->mark_node(i);
            }
            old_source[i] = source;
            old_sink[i] = sink;
        }
    }
}
