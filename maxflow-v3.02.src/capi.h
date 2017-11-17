
void * new_graph(int nodes, int edges);
void add_capacity(void *p, int nodes, int *source, int *sink);
void add_edges(void *p, int edges, int *head, int *tail, int *cap);
void add_edge(void *p, int head, int tail, int cap);
int maxflow(void *p, int reuse);
void get_segments(void *p, int nodes, int *segs);
void free_graph(void *p);
void set_capacity_from_image(void *p, int nodes, int *old_source, int *old_sink, unsigned char *image, float pri, float scale);

