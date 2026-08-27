#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

typedef struct {
    int n, m;
    int *u, *v;
    int *offset, *incident;
    unsigned char *active;
} Graph;

typedef struct {
    int *edge, *key;
    uint64_t *tie;
    int size, capacity;
} Heap;

static uint64_t rng_state;

static uint64_t rng_u64(void)
{
    uint64_t z = (rng_state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30)) * UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27)) * UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}

static double rng_open01(void)
{
    return ((rng_u64() >> 11) + 0.5) * (1.0 / 9007199254740992.0);
}

static uint64_t rng_below(uint64_t bound)
{
    uint64_t x, limit = UINT64_MAX - UINT64_MAX % bound;
    do x = rng_u64(); while (x >= limit);
    return x % bound;
}

static void *xcalloc(size_t count, size_t size)
{
    void *p = calloc(count, size);
    if (!p && count != 0) {
        fprintf(stderr, "out of memory\n");
        exit(EXIT_FAILURE);
    }
    return p;
}

static int other(const Graph *g, int e, int x)
{
    return g->u[e] == x ? g->v[e] : g->u[e];
}

static int cmp_u64(const void *a, const void *b)
{
    uint64_t x = *(const uint64_t *)a, y = *(const uint64_t *)b;
    return (x > y) - (x < y);
}

static void build_csr(Graph *g)
{
    int *degree = xcalloc((size_t)g->n, sizeof(*degree));
    g->active = xcalloc((size_t)g->m, sizeof(*g->active));
    for (int e = 0; e < g->m; e++) {
        degree[g->u[e]]++;
        degree[g->v[e]]++;
        g->active[e] = 1;
    }
    g->offset = xcalloc((size_t)g->n + 1, sizeof(*g->offset));
    for (int i = 0; i < g->n; i++) g->offset[i + 1] = g->offset[i] + degree[i];
    g->incident = xcalloc((size_t)2 * g->m, sizeof(*g->incident));
    int *next = xcalloc((size_t)g->n, sizeof(*next));
    memcpy(next, g->offset, (size_t)g->n * sizeof(*next));
    for (int e = 0; e < g->m; e++) {
        g->incident[next[g->u[e]]++] = e;
        g->incident[next[g->v[e]]++] = e;
    }
    free(next);
    free(degree);
}

static Graph read_graph(const char *path)
{
    Graph g = {0};
    FILE *in = fopen(path, "r");
    if (!in) {
        perror(path);
        exit(EXIT_FAILURE);
    }
    if (fscanf(in, "%d %d", &g.n, &g.m) != 2 || g.n < 1 || g.m < 0) {
        fprintf(stderr, "%s: expected header 'N M'\n", path);
        exit(EXIT_FAILURE);
    }
    g.u = xcalloc((size_t)g.m, sizeof(*g.u));
    g.v = xcalloc((size_t)g.m, sizeof(*g.v));
    uint64_t *keys = xcalloc((size_t)g.m, sizeof(*keys));
    for (int e = 0; e < g.m; e++) {
        if (fscanf(in, "%d %d", &g.u[e], &g.v[e]) != 2) {
            fprintf(stderr, "%s: expected %d edges, stopped at edge %d\n",
                    path, g.m, e);
            exit(EXIT_FAILURE);
        }
        if (g.u[e] < 0 || g.u[e] >= g.n || g.v[e] < 0 || g.v[e] >= g.n ||
            g.u[e] == g.v[e]) {
            fprintf(stderr, "%s: invalid edge %d: %d %d\n",
                    path, e, g.u[e], g.v[e]);
            exit(EXIT_FAILURE);
        }
        int lo = g.u[e] < g.v[e] ? g.u[e] : g.v[e];
        int hi = g.u[e] < g.v[e] ? g.v[e] : g.u[e];
        keys[e] = (uint64_t)(unsigned)lo * (uint64_t)(unsigned)g.n +
                  (uint64_t)(unsigned)hi;
    }
    char extra[2];
    if (fscanf(in, "%1s", extra) == 1) {
        fprintf(stderr, "%s: extra data after the declared edges\n", path);
        exit(EXIT_FAILURE);
    }
    fclose(in);
    qsort(keys, (size_t)g.m, sizeof(*keys), cmp_u64);
    for (int e = 1; e < g.m; e++) {
        if (keys[e] == keys[e - 1]) {
            fprintf(stderr, "%s: duplicate edges are not supported\n", path);
            exit(EXIT_FAILURE);
        }
    }
    free(keys);
    build_csr(&g);
    return g;
}

static Graph generate_er_graph(int n, double mean_degree, uint64_t seed)
{
    if (n < 2) {
        fprintf(stderr, "N must be at least 2\n");
        exit(EXIT_FAILURE);
    }
    if (mean_degree < 0.0 || mean_degree > n - 1.0) {
        fprintf(stderr, "mean degree must be in [0, N-1]\n");
        exit(EXIT_FAILURE);
    }
    double p = mean_degree / (n - 1.0);
    double log_q = p < 1.0 ? log1p(-p) : 0.0;
    long m = 0;
    rng_state = seed;
    if (p > 0.0) {
        int v = 1;
        long w = -1;
        while (v < n) {
            w += p == 1.0 ? 1 : 1 + (long)floor(log(rng_open01()) / log_q);
            while (w >= v && v < n) { w -= v; v++; }
            if (v < n && ++m > INT_MAX) {
                fprintf(stderr, "generated graph has too many edges\n");
                exit(EXIT_FAILURE);
            }
        }
    }

    Graph g = {.n = n, .m = (int)m};
    g.u = xcalloc((size_t)g.m, sizeof(*g.u));
    g.v = xcalloc((size_t)g.m, sizeof(*g.v));
    rng_state = seed;
    int e = 0;
    if (p > 0.0) {
        int v = 1;
        long w = -1;
        while (v < n) {
            w += p == 1.0 ? 1 : 1 + (long)floor(log(rng_open01()) / log_q);
            while (w >= v && v < n) { w -= v; v++; }
            if (v < n) {
                g.u[e] = (int)w;
                g.v[e] = v;
                e++;
            }
        }
    }
    build_csr(&g);
    fprintf(stderr, "generated ER G(N,p): N=%d M=%d p=%.8g seed=%llu\n",
            n, g.m, p, (unsigned long long)seed);
    return g;
}

static void free_graph(Graph *g)
{
    free(g->u); free(g->v); free(g->offset); free(g->incident); free(g->active);
}

/* Peel vertices of degree < 2. Returns the number of vertices in the 2-core. */
static int find_two_core(const Graph *g, unsigned char *in_core, int *core_degree,
                         int *queue)
{
    int head = 0, tail = 0, count = g->n;
    memset(in_core, 1, (size_t)g->n);
    memset(core_degree, 0, (size_t)g->n * sizeof(*core_degree));
    for (int e = 0; e < g->m; e++) if (g->active[e]) {
        core_degree[g->u[e]]++;
        core_degree[g->v[e]]++;
    }
    for (int i = 0; i < g->n; i++) if (core_degree[i] < 2) queue[tail++] = i;
    while (head < tail) {
        int x = queue[head++];
        if (!in_core[x]) continue;
        in_core[x] = 0;
        count--;
        for (int k = g->offset[x]; k < g->offset[x + 1]; k++) {
            int e = g->incident[k];
            if (!g->active[e]) continue;
            int y = other(g, e, x);
            if (in_core[y] && --core_degree[y] == 1) queue[tail++] = y;
        }
    }
    return count;
}

static int heap_better(const Heap *h, int a, int b)
{
    if (h->key[a] != h->key[b]) return h->key[a] > h->key[b];
    return h->tie[a] < h->tie[b];
}

static void heap_push(Heap *h, int key, uint64_t tie, int edge)
{
    if (h->size == h->capacity) {
        h->capacity = h->capacity ? 2 * h->capacity : 1024;
        h->key = realloc(h->key, (size_t)h->capacity * sizeof(*h->key));
        h->edge = realloc(h->edge, (size_t)h->capacity * sizeof(*h->edge));
        h->tie = realloc(h->tie, (size_t)h->capacity * sizeof(*h->tie));
        if (!h->key || !h->edge || !h->tie) {
            fprintf(stderr, "out of memory\n");
            exit(EXIT_FAILURE);
        }
    }
    int i = h->size++;
    h->key[i] = key;
    h->tie[i] = tie;
    h->edge[i] = edge;
    while (i > 0) {
        int parent = (i - 1) / 2;
        if (!heap_better(h, i, parent)) break;
        int t = h->key[i]; h->key[i] = h->key[parent]; h->key[parent] = t;
        t = h->edge[i]; h->edge[i] = h->edge[parent]; h->edge[parent] = t;
        uint64_t z = h->tie[i]; h->tie[i] = h->tie[parent]; h->tie[parent] = z;
        i = parent;
    }
}

static int heap_pop(Heap *h, int *key)
{
    if (!h->size) return -1;
    int edge = h->edge[0];
    *key = h->key[0];
    h->size--;
    if (h->size) {
        h->edge[0] = h->edge[h->size];
        h->key[0] = h->key[h->size];
        h->tie[0] = h->tie[h->size];
        for (int i = 0;;) {
            int left = 2 * i + 1, right = left + 1, best = i;
            if (left < h->size && heap_better(h, left, best)) best = left;
            if (right < h->size && heap_better(h, right, best)) best = right;
            if (best == i) break;
            int t = h->key[i]; h->key[i] = h->key[best]; h->key[best] = t;
            t = h->edge[i]; h->edge[i] = h->edge[best]; h->edge[best] = t;
            uint64_t z = h->tie[i]; h->tie[i] = h->tie[best]; h->tie[best] = z;
            i = best;
        }
    }
    return edge;
}

static double l1_normalize(double *x, const unsigned char *edge_ok, int m)
{
    double norm = 0.0;
    for (int e = 0; e < m; e++) if (edge_ok[e])
        norm += fabs(x[2 * e]) + fabs(x[2 * e + 1]);
    if (!(norm > 0.0)) return 0.0;
    for (int e = 0; e < m; e++) {
        if (edge_ok[e]) {
            x[2 * e] /= norm;
            x[2 * e + 1] /= norm;
        } else {
            x[2 * e] = x[2 * e + 1] = 0.0;
        }
    }
    return norm;
}

static int normalize_or_reseed(double *x, const unsigned char *edge_ok, int m)
{
    double norm = l1_normalize(x, edge_ok, m);
    if (norm > 0.0 && isfinite(1.0 / norm)) return 1;
    for (int e = 0; e < m; e++) {
        if (edge_ok[e]) x[2 * e] = x[2 * e + 1] = 1.0;
    }
    return l1_normalize(x, edge_ok, m) > 0.0;
}

static void apply_nb_right(const Graph *g, const unsigned char *edge_ok,
                           const double *x, double *out, double *sum)
{
    memset(sum, 0, (size_t)g->n * sizeof(*sum));
    memset(out, 0, (size_t)2 * g->m * sizeof(*out));
    for (int e = 0; e < g->m; e++) if (edge_ok[e]) {
        sum[g->u[e]] += x[2 * e];
        sum[g->v[e]] += x[2 * e + 1];
    }
    for (int e = 0; e < g->m; e++) if (edge_ok[e]) {
        out[2 * e] = sum[g->v[e]] - x[2 * e + 1] + x[2 * e];
        out[2 * e + 1] = sum[g->u[e]] - x[2 * e] + x[2 * e + 1];
    }
}

static void apply_nb_left(const Graph *g, const unsigned char *edge_ok,
                          const double *x, double *out, double *sum)
{
    memset(sum, 0, (size_t)g->n * sizeof(*sum));
    memset(out, 0, (size_t)2 * g->m * sizeof(*out));
    for (int e = 0; e < g->m; e++) if (edge_ok[e]) {
        sum[g->u[e]] += x[2 * e + 1];
        sum[g->v[e]] += x[2 * e];
    }
    for (int e = 0; e < g->m; e++) if (edge_ok[e]) {
        out[2 * e] = sum[g->u[e]] - x[2 * e + 1] + x[2 * e];
        out[2 * e + 1] = sum[g->v[e]] - x[2 * e] + x[2 * e + 1];
    }
}

/* NBC h(e)=l_e r_e+l_reverse r_reverse, using shifted power iteration. */
enum { NBC_MIN_ITERS = 8, NBC_MAX_ITERS = 20000 };

static int nbc_scores(const Graph *g, const unsigned char *edge_ok,
                      double *score, double *r, double *l, double *next_r,
                      double *next_l, double *sum)
{
    double *r_store = r, *l_store = l;
    if (!normalize_or_reseed(r, edge_ok, g->m) ||
        !normalize_or_reseed(l, edge_ok, g->m)) return 0;
    int used_iters = 0;
    for (int iter = 1; iter <= NBC_MAX_ITERS; iter++) {
        used_iters = iter;
        apply_nb_right(g, edge_ok, r, next_r, sum);
        apply_nb_left(g, edge_ok, l, next_l, sum);
        if (!l1_normalize(next_r, edge_ok, g->m) ||
            !l1_normalize(next_l, edge_ok, g->m)) break;
        double diff_r = 0.0, diff_l = 0.0;
        for (int e = 0; e < g->m; e++) if (edge_ok[e]) {
            diff_r += fabs(next_r[2 * e] - r[2 * e]);
            diff_r += fabs(next_r[2 * e + 1] - r[2 * e + 1]);
            diff_l += fabs(next_l[2 * e] - l[2 * e]);
            diff_l += fabs(next_l[2 * e + 1] - l[2 * e + 1]);
        }
        double *tmp = r; r = next_r; next_r = tmp;
        tmp = l; l = next_l; next_l = tmp;
        if (iter >= NBC_MIN_ITERS && fmax(diff_r, diff_l) < 1e-10) break;
    }
    for (int e = 0; e < g->m; e++)
        score[e] = edge_ok[e] ? l[2 * e] * r[2 * e] +
                                l[2 * e + 1] * r[2 * e + 1] : 0.0;
    if (r != r_store) memcpy(r_store, r, (size_t)2 * g->m * sizeof(*r));
    if (l != l_store) memcpy(l_store, l, (size_t)2 * g->m * sizeof(*l));
    return used_iters;
}

/* Production-compatible deterministic tie-breaking for a=+/-infinity. */
static int choose_extreme(const int *candidate, int count, const double *score,
                          int maximize)
{
    double best = maximize ? -DBL_MAX : DBL_MAX;
    int chosen = -1;
    for (int i = 0; i < count; i++) {
        int e = candidate[i];
        int better = maximize ? score[e] > best : score[e] < best;
        if (better || (score[e] == best && (chosen < 0 || e < chosen))) {
            best = score[e];
            chosen = e;
        }
    }
    return chosen;
}

/* Stable evaluation of weights score^a via logarithms. */
static int choose_finite(const int *candidate, int count, const double *score,
                         double a)
{
    if (a == 0.0) return candidate[rng_below((uint64_t)count)];
    double max_log_weight = -DBL_MAX;
    int positive = 0;
    for (int i = 0; i < count; i++) if (score[candidate[i]] > 0.0) {
        double z = a * log(score[candidate[i]]);
        if (z > max_log_weight) max_log_weight = z;
        positive++;
    }
    if (!positive) return candidate[rng_below((uint64_t)count)];
    double total = 0.0;
    for (int i = 0; i < count; i++) if (score[candidate[i]] > 0.0)
        total += exp(a * log(score[candidate[i]]) - max_log_weight);
    double target = rng_open01() * total;
    int last = -1;
    for (int i = 0; i < count; i++) if (score[candidate[i]] > 0.0) {
        int e = candidate[i];
        last = e;
        target -= exp(a * log(score[e]) - max_log_weight);
        if (target <= 0.0) return e;
    }
    return last;
}

static int sum_key(const Graph *g, const int *core_degree, int e, int minimize)
{
    int key = core_degree[g->u[e]] + core_degree[g->v[e]];
    return minimize ? -key : key;
}

static uint64_t endpoint_tie(const Graph *g, int e)
{
    int u = g->u[e], v = g->v[e];
    if (u > v) { int t = u; u = v; v = t; }
    return (uint64_t)(unsigned)u * (uint64_t)(unsigned)g->n +
           (uint64_t)(unsigned)v;
}

static void rebuild_sum_heap(Heap *heap, const Graph *g,
                             const unsigned char *in_core,
                             const int *core_degree, int minimize)
{
    heap->size = 0;
    for (int e = 0; e < g->m; e++)
        if (g->active[e] && in_core[g->u[e]] && in_core[g->v[e]])
            heap_push(heap, sum_key(g, core_degree, e, minimize),
                      endpoint_tie(g, e), e);
}

static int choose_sum_heap(Heap *heap, const Graph *g,
                           const unsigned char *in_core,
                           const int *core_degree, int minimize)
{
    for (;;) {
        int stored;
        int e = heap_pop(heap, &stored);
        if (e < 0) return -1;
        if (!g->active[e] || !in_core[g->u[e]] || !in_core[g->v[e]]) continue;
        int current = sum_key(g, core_degree, e, minimize);
        if (stored != current) {
            heap_push(heap, current, endpoint_tie(g, e), e);
            continue;
        }
        return e;
    }
}

/* Delete one core edge and perform only the resulting local peeling cascade. */
static int remove_core_edge(Graph *g, int chosen, unsigned char *in_core,
                            int *core_degree, int *queue, int core_size,
                            int *touched, int *touch_mark, int epoch,
                            int *ntouched_out)
{
    int head = 0, tail = 0, ntouched = 0;
    int u = g->u[chosen], v = g->v[chosen];
    g->active[chosen] = 0;
    if (!in_core[u] || !in_core[v]) {
        *ntouched_out = 0;
        return core_size;
    }

    int endpoints[2] = {u, v};
    for (int i = 0; i < 2; i++) {
        int x = endpoints[i];
        core_degree[x]--;
        if (touch_mark[x] != epoch) {
            touch_mark[x] = epoch;
            touched[ntouched++] = x;
        }
        if (core_degree[x] == 1) queue[tail++] = x;
    }
    while (head < tail) {
        int x = queue[head++];
        if (!in_core[x] || core_degree[x] >= 2) continue;
        in_core[x] = 0;
        core_size--;
        for (int k = g->offset[x]; k < g->offset[x + 1]; k++) {
            int e = g->incident[k];
            if (!g->active[e]) continue;
            int y = other(g, e, x);
            if (!in_core[y]) continue;
            core_degree[y]--;
            if (touch_mark[y] != epoch) {
                touch_mark[y] = epoch;
                touched[ntouched++] = y;
            }
            if (core_degree[y] == 1) queue[tail++] = y;
        }
    }
    *ntouched_out = ntouched;
    return core_size;
}

static void refresh_sum_heap(Heap *heap, const Graph *g,
                             const unsigned char *in_core,
                             const int *core_degree, int minimize,
                             const int *touched, int ntouched,
                             int *edge_mark, int epoch)
{
    for (int i = 0; i < ntouched; i++) {
        int x = touched[i];
        if (!in_core[x]) continue;
        for (int k = g->offset[x]; k < g->offset[x + 1]; k++) {
            int e = g->incident[k];
            if (edge_mark[e] == epoch || !g->active[e] ||
                !in_core[g->u[e]] || !in_core[g->v[e]]) continue;
            edge_mark[e] = epoch;
            heap_push(heap, sum_key(g, core_degree, e, minimize),
                      endpoint_tie(g, e), e);
        }
    }
}

static int uf_find(int *parent, int x)
{
    while (parent[x] != x) {
        parent[x] = parent[parent[x]];
        x = parent[x];
    }
    return x;
}

static void giant_component_curve(const Graph *g, const int *order, int *curve,
                                  int *parent, int *size)
{
    for (int i = 0; i < g->n; i++) {
        parent[i] = i;
        size[i] = 1;
    }
    int largest = g->n ? 1 : 0;
    curve[g->m] = largest;
    for (int k = g->m - 1; k >= 0; k--) {
        int e = order[k];
        int a = uf_find(parent, g->u[e]);
        int b = uf_find(parent, g->v[e]);
        if (a != b) {
            if (size[a] < size[b]) { int t = a; a = b; b = t; }
            parent[b] = a;
            size[a] += size[b];
            if (size[a] > largest) largest = size[a];
        }
        curve[k] = largest;
    }
}

static double parse_bias(const char *text)
{
    if (!strcmp(text, "inf") || !strcmp(text, "+inf")) return INFINITY;
    if (!strcmp(text, "-inf")) return -INFINITY;
    char *end;
    errno = 0;
    double a = strtod(text, &end);
    if (errno || *text == '\0' || *end != '\0' || !isfinite(a)) {
        fprintf(stderr, "invalid bias a: %s (use a number, inf, or -inf)\n", text);
        exit(EXIT_FAILURE);
    }
    return a;
}

static uint64_t parse_seed(const char *text)
{
    char *end;
    errno = 0;
    unsigned long long seed = strtoull(text, &end, 10);
    if (errno || *text == '\0' || *end != '\0') {
        fprintf(stderr, "invalid seed: %s\n", text);
        exit(EXIT_FAILURE);
    }
    return (uint64_t)seed;
}

static int parse_n(const char *text)
{
    char *end;
    errno = 0;
    long value = strtol(text, &end, 10);
    if (errno || *text == '\0' || *end != '\0' || value < 2 || value > INT_MAX) {
        fprintf(stderr, "invalid N: %s\n", text);
        exit(EXIT_FAILURE);
    }
    return (int)value;
}

static double parse_mean_degree(const char *text)
{
    char *end;
    errno = 0;
    double value = strtod(text, &end);
    if (errno || *text == '\0' || *end != '\0' || !isfinite(value)) {
        fprintf(stderr, "invalid mean degree: %s\n", text);
        exit(EXIT_FAILURE);
    }
    return value;
}

int main(int argc, char **argv)
{
    Graph g;
    const char *score_text, *bias_text, *removal_seed_text;
    const char *results_path, *sequence_path;
    if (argc == 7 && strcmp(argv[1], "--er")) {
        g = read_graph(argv[1]);
        score_text = argv[2];
        bias_text = argv[3];
        removal_seed_text = argv[4];
        results_path = argv[5];
        sequence_path = argv[6];
    } else if (argc == 10 && !strcmp(argv[1], "--er")) {
        int n = parse_n(argv[2]);
        double mean_degree = parse_mean_degree(argv[3]);
        uint64_t graph_seed = parse_seed(argv[4]);
        g = generate_er_graph(n, mean_degree, graph_seed);
        score_text = argv[5];
        bias_text = argv[6];
        removal_seed_text = argv[7];
        results_path = argv[8];
        sequence_path = argv[9];
    } else {
        fprintf(stderr,
                "usage:\n"
                "  %s network.txt {nbc|2CDC} a seed results.csv sequence.txt\n"
                "  %s --er N mean_degree graph_seed {nbc|2CDC} a seed "
                "results.csv sequence.txt\n", argv[0], argv[0]);
        return EXIT_FAILURE;
    }
    int use_nbc;
    if (!strcmp(score_text, "nbc")) use_nbc = 1;
    else if (!strcmp(score_text, "2CDC")) use_nbc = 0;
    else {
        fprintf(stderr, "score must be 'nbc' or '2CDC'\n");
        return EXIT_FAILURE;
    }
    double a = parse_bias(bias_text);
    rng_state = parse_seed(removal_seed_text);
    FILE *csv = fopen(results_path, "w");
    FILE *seq = fopen(sequence_path, "w");
    if (!csv || !seq) {
        perror(!csv ? results_path : sequence_path);
        return EXIT_FAILURE;
    }

    unsigned char *in_core = xcalloc((size_t)g.n, 1);
    unsigned char *edge_ok = xcalloc((size_t)g.m, 1);
    int *core_degree = xcalloc((size_t)g.n, sizeof(*core_degree));
    int *queue = xcalloc((size_t)g.n, sizeof(*queue));
    int *candidate = xcalloc((size_t)g.m, sizeof(*candidate));
    int *order = xcalloc((size_t)g.m, sizeof(*order));
    int *core_curve = xcalloc((size_t)g.m + 1, sizeof(*core_curve));
    int *gc_curve = xcalloc((size_t)g.m + 1, sizeof(*gc_curve));
    int *touched = xcalloc((size_t)g.n, sizeof(*touched));
    int *touch_mark = xcalloc((size_t)g.n, sizeof(*touch_mark));
    int *edge_mark = xcalloc((size_t)g.m, sizeof(*edge_mark));
    int *uf_parent = xcalloc((size_t)g.n, sizeof(*uf_parent));
    int *uf_size = xcalloc((size_t)g.n, sizeof(*uf_size));
    double *score = xcalloc((size_t)g.m, sizeof(*score));
    double *r = xcalloc((size_t)2 * g.m, sizeof(*r));
    double *l = xcalloc((size_t)2 * g.m, sizeof(*l));
    double *next_r = xcalloc((size_t)2 * g.m, sizeof(*next_r));
    double *next_l = xcalloc((size_t)2 * g.m, sizeof(*next_l));
    double *sum = xcalloc((size_t)g.n, sizeof(*sum));
    for (int e = 0; e < g.m; e++) {
        r[2 * e] = r[2 * e + 1] = 1.0;
        l[2 * e] = l[2 * e + 1] = 1.0;
    }

    fprintf(seq, "# N %d M %d\n", g.n, g.m);
    fprintf(seq, "step edge_id u v\n");
    int core_size = find_two_core(&g, in_core, core_degree, queue);
    core_curve[0] = core_size;
    int step = 0;
    int nbc_it_max = 0, nbc_capped = 0;
    int pool_count = 0;
    Heap heap = {0};
    int heap_mode = !use_nbc && isinf(a);
    int minimize = a < 0.0;
    if (!use_nbc && a == 0.0) {
        for (int e = 0; e < g.m; e++)
            if (g.active[e] && in_core[g.u[e]] && in_core[g.v[e]])
                candidate[pool_count++] = e;
    } else if (heap_mode) {
        rebuild_sum_heap(&heap, &g, in_core, core_degree, minimize);
    }

    while (core_size > 0) {
        int count = 0;
        int chosen = -1;
        if (use_nbc) {
            for (int e = 0; e < g.m; e++) {
                edge_ok[e] = g.active[e] && in_core[g.u[e]] && in_core[g.v[e]];
                if (edge_ok[e]) candidate[count++] = e;
            }
            int used = nbc_scores(&g, edge_ok, score, r, l,
                                  next_r, next_l, sum);
            if (used > nbc_it_max) nbc_it_max = used;
            if (used == NBC_MAX_ITERS) nbc_capped++;
            if (isinf(a)) chosen = choose_extreme(candidate, count, score, a > 0.0);
            else chosen = choose_finite(candidate, count, score, a);
        } else if (heap_mode) {
            chosen = choose_sum_heap(&heap, &g, in_core, core_degree, minimize);
        } else if (a == 0.0) {
            while (pool_count > 0) {
                int i = (int)rng_below((uint64_t)pool_count);
                int e = candidate[i];
                candidate[i] = candidate[--pool_count];
                if (g.active[e] && in_core[g.u[e]] && in_core[g.v[e]]) {
                    chosen = e;
                    break;
                }
            }
        } else {
            for (int e = 0; e < g.m; e++) {
                if (!g.active[e] || !in_core[g.u[e]] || !in_core[g.v[e]]) continue;
                candidate[count++] = e;
                score[e] = (double)core_degree[g.u[e]] + core_degree[g.v[e]];
            }
            chosen = choose_finite(candidate, count, score, a);
        }
        if (chosen < 0) {
            fprintf(stderr, "internal error: no removable 2-core edge at step %d\n", step);
            return EXIT_FAILURE;
        }

        order[step] = chosen;
        fprintf(seq, "%d %d %d %d\n", step + 1, chosen,
                g.u[chosen], g.v[chosen]);
        int epoch = step + 1;
        int ntouched;
        core_size = remove_core_edge(&g, chosen, in_core, core_degree, queue,
                                     core_size, touched, touch_mark, epoch,
                                     &ntouched);
        step++;
        core_curve[step] = core_size;
        if (heap_mode) {
            refresh_sum_heap(&heap, &g, in_core, core_degree, minimize,
                             touched, ntouched, edge_mark, epoch);
            if (heap.size > 3 * g.m)
                rebuild_sum_heap(&heap, &g, in_core, core_degree, minimize);
        }
    }
    int core_steps = step;

    /* Production-style forest tail: shuffle once, then evaluate GC in reverse. */
    int tail_count = 0;
    for (int e = 0; e < g.m; e++)
        if (g.active[e]) candidate[tail_count++] = e;
    for (int i = tail_count - 1; i > 0; i--) {
        int j = (int)rng_below((uint64_t)(i + 1));
        int t = candidate[i]; candidate[i] = candidate[j]; candidate[j] = t;
    }
    for (int i = 0; i < tail_count; i++) {
        int e = candidate[i];
        order[step] = e;
        fprintf(seq, "%d %d %d %d\n", step + 1, e, g.u[e], g.v[e]);
        step++;
    }
    if (step != g.m) {
        fprintf(stderr, "internal error: removal order has %d of %d edges\n", step, g.m);
        return EXIT_FAILURE;
    }

    giant_component_curve(&g, order, gc_curve, uf_parent, uf_size);
    fprintf(csv, "removed,p,GC,2C\n");
    for (int k = 0; k <= g.m; k++) {
        double p = g.m ? (double)k / g.m : 0.0;
        fprintf(csv, "%d,%.17g,%.17g,%.17g\n", k, p,
                (double)gc_curve[k] / g.n, (double)core_curve[k] / g.n);
    }

    fclose(csv);
    fclose(seq);
    fprintf(stderr, "wrote %s and %s (%s, a=%s, seed=%s)\n",
            results_path, sequence_path, score_text, bias_text, removal_seed_text);
    if (use_nbc)
        fprintf(stderr, "NBC iterations: maximum=%d, capped_steps=%d/%d\n",
                nbc_it_max, nbc_capped, core_steps);
    free(in_core); free(edge_ok); free(core_degree); free(queue); free(candidate);
    free(order); free(core_curve); free(gc_curve); free(touched); free(touch_mark);
    free(edge_mark); free(uf_parent); free(uf_size); free(score); free(r); free(l);
    free(next_r); free(next_l); free(sum); free(heap.edge); free(heap.key);
    free(heap.tie);
    free_graph(&g);
    return EXIT_SUCCESS;
}
