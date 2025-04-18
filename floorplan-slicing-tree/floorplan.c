#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdbool.h>
#include "floorplan.h"

// Slicing tree node for floorplanning
typedef struct TreeNode {
    int isLeaf;       // 1 if leaf (module), 0 if internal node
    char op;          // 'H' for horizontal cut, 'V' for vertical cut
    struct TreeNode *left, *right;
    Module *module;   // non-NULL for leaf nodes
    int w, h;         // subtree width and height
} TreeNode;

// Create a leaf node for a module
static TreeNode *create_leaf(Module *m) {
    TreeNode *node = (TreeNode *)malloc(sizeof(TreeNode));
    node->isLeaf = 1;
    node->op = 0;
    node->left = node->right = NULL;
    node->module = m;
    node->w = m->w;
    node->h = m->h;
    return node;
}

// Create an internal node with given cut operator
static TreeNode *create_internal(char op, TreeNode *left, TreeNode *right) {
    TreeNode *node = (TreeNode *)malloc(sizeof(TreeNode));
    node->isLeaf = 0;
    node->op = op;
    node->left = left;
    node->right = right;
    node->module = NULL;
    node->w = node->h = 0;
    return node;
}

// Post-order evaluation of subtree dimensions
static void evaluate_node(TreeNode *node) {
    if (node->isLeaf) return;
    evaluate_node(node->left);
    evaluate_node(node->right);
    if (node->op == 'H') {
        node->w = node->left->w > node->right->w ? node->left->w : node->right->w;
        node->h = node->left->h + node->right->h;
    } else {
        node->w = node->left->w + node->right->w;
        node->h = node->left->h > node->right->h ? node->left->h : node->right->h;
    }
}

// Assign absolute positions to modules based on evaluated tree
static void assign_positions(TreeNode *node, int x, int y) {
    if (node->isLeaf) {
        node->module->x = x;
        node->module->y = y;
        return;
    }
    if (node->op == 'H') {
        // horizontal cut: left subtree below, right above
        assign_positions(node->left, x, y);
        assign_positions(node->right, x, y + node->left->h);
    } else {
        // vertical cut: left subtree left, right subtree right
        assign_positions(node->left, x, y);
        assign_positions(node->right, x + node->left->w, y);
    }
}

// Recursively free slicing tree nodes
static void free_tree(TreeNode *node) {
    if (!node) return;
    if (!node->isLeaf) {
        free_tree(node->left);
        free_tree(node->right);
    }
    free(node);
}
// Build a slicing tree for soft modules using greedy divide-and-conquer
static TreeNode *build_slicing_tree(Module *mods[], int start, int end) {
    if (start + 1 == end) {
        return create_leaf(mods[start]);
    }
    int mid = (start + end) / 2;
    TreeNode *left = build_slicing_tree(mods, start, mid);
    TreeNode *right = build_slicing_tree(mods, mid, end);
    int w_h = left->w > right->w ? left->w : right->w;
    int h_h = left->h + right->h;
    long area_h = (long)w_h * h_h;
    int w_v = left->w + right->w;
    int h_v = left->h > right->h ? left->h : right->h;
    long area_v = (long)w_v * h_v;
    char op = (area_h < area_v ? 'H' : 'V');
    TreeNode *node = create_internal(op, left, right);
    if (op == 'H') {
        node->w = w_h;
        node->h = h_h;
    } else {
        node->w = w_v;
        node->h = h_v;
    }
    return node;
}

int parse_input(const char *filename, Floorplan *fp) {
    FILE *in = fopen(filename, "r");
    if (!in) {
        fprintf(stderr, "Error: cannot open input file %s\n", filename);
        return -1;
    }
    char token[256];
    if (fscanf(in, "%255s %d %d", token, &fp->chipWidth, &fp->chipHeight) != 3 || strcmp(token, "ChipSize") != 0) {
        fprintf(stderr, "Error: expected 'ChipSize <width> <height>'\n");
        fclose(in);
        return -1;
    }
    int nsoft = 0;
    if (fscanf(in, "%255s %d", token, &nsoft) != 2 || strcmp(token, "NumSoftModules") != 0) {
        fprintf(stderr, "Error: expected 'NumSoftModules <count>'\n");
        fclose(in);
        return -1;
    }
    fp->nsoft = nsoft;
    fp->modules = malloc(sizeof(Module) * (nsoft));
    if (!fp->modules) { fclose(in); return -1; }
    for (int i = 0; i < nsoft; i++) {
        if (fscanf(in, "%255s", token) != 1 || strcmp(token, "SoftModule") != 0) {
            fprintf(stderr, "Error: expected 'SoftModule'\n"); fclose(in); return -1;
        }
        char name[256]; int area;
        if (fscanf(in, "%255s %d", name, &area) != 2) {
            fprintf(stderr, "Error: invalid 'SoftModule <name> <area>'\n"); fclose(in); return -1;
        }
        fp->modules[i].name = strdup(name);
        fp->modules[i].min_area = area;
        fp->modules[i].isFixed = 0;
        fp->modules[i].x = fp->modules[i].y = fp->modules[i].w = fp->modules[i].h = 0;
    }
    int nf = 0;
    if (fscanf(in, "%255s %d", token, &nf) != 2 || strcmp(token, "NumFixedModules") != 0) {
        fprintf(stderr, "Error: expected 'NumFixedModules <count>'\n"); fclose(in); return -1;
    }
    fp->nf = nf;
    fp->modules = realloc(fp->modules, sizeof(Module) * (nsoft + nf));
    if (!fp->modules) { fclose(in); return -1; }
    for (int i = 0; i < nf; i++) {
        if (fscanf(in, "%255s", token) != 1 || strcmp(token, "FixedModule") != 0) {
            fprintf(stderr, "Error: expected 'FixedModule'\n"); fclose(in); return -1;
        }
        char name[256]; int x, y, w, h;
        if (fscanf(in, "%255s %d %d %d %d", name, &x, &y, &w, &h) != 5) {
            fprintf(stderr, "Error: invalid 'FixedModule <name> <x> <y> <w> <h>'\n"); fclose(in); return -1;
        }
        int idx = nsoft + i;
        fp->modules[idx].name = strdup(name);
        fp->modules[idx].min_area = 0;
        fp->modules[idx].isFixed = 1;
        fp->modules[idx].x = x;
        fp->modules[idx].y = y;
        fp->modules[idx].w = w;
        fp->modules[idx].h = h;
    }
    fp->nModules = nsoft + nf;
    int nn = 0;
    if (fscanf(in, "%255s %d", token, &nn) != 2 || strcmp(token, "NumNets") != 0) {
        fprintf(stderr, "Error: expected 'NumNets <count>'\n"); fclose(in); return -1;
    }
    fp->nNets = nn;
    fp->nets = malloc(sizeof(Net) * nn);
    if (!fp->nets) { fclose(in); return -1; }
    for (int i = 0; i < nn; i++) {
        if (fscanf(in, "%255s", token) != 1 || strcmp(token, "Net") != 0) {
            fprintf(stderr, "Error: expected 'Net'\n"); fclose(in); return -1;
        }
        char n1[256], n2[256]; int wt;
        if (fscanf(in, "%255s %255s %d", n1, n2, &wt) != 3) {
            fprintf(stderr, "Error: invalid 'Net <mod1> <mod2> <weight>'\n"); fclose(in); return -1;
        }
        int u = -1, v = -1;
        for (int j = 0; j < fp->nModules; j++) {
            if (strcmp(fp->modules[j].name, n1) == 0) u = j;
            if (strcmp(fp->modules[j].name, n2) == 0) v = j;
        }
        if (u < 0 || v < 0) {
            fprintf(stderr, "Error: Net refers to unknown module '%s' or '%s'\n", n1, n2);
            fclose(in);
            return -1;
        }
        fp->nets[i].u = u;
        fp->nets[i].v = v;
        fp->nets[i].weight = wt;
    }
    fclose(in);
    return 0;
}

void compute_floorplan(Floorplan *fp) {
    // determine dimensions for soft modules (square approximation within aspect ratio limits)
    for (int i = 0; i < fp->nsoft; i++) {
        Module *m = &fp->modules[i];
        int side = (int)ceil(sqrt((double)m->min_area));
        m->w = side;
        m->h = side;
    }

    int n = fp->nsoft;
    if (n <= 0) return;
    // build slicing tree for soft modules (greedy divide-and-conquer)
    Module *softMods[n];
    for (int i = 0; i < n; i++) {
        softMods[i] = &fp->modules[i];
    }
    TreeNode *root = build_slicing_tree(softMods, 0, n);
    // compute base y-coordinate to avoid fixed modules
    int base_y = 0;
    for (int j = fp->nsoft; j < fp->nModules; j++) {
        Module *f = &fp->modules[j];
        int top = f->y + f->h;
        if (top > base_y) base_y = top;
    }
    // assign absolute positions starting at (0, base_y)
    assign_positions(root, 0, base_y);
    // free slicing tree nodes
    free_tree(root);
}

int compute_wirelength(const Floorplan *fp) {
    int total = 0;
    for (int i = 0; i < fp->nNets; i++) {
        const Net *net = &fp->nets[i];
        const Module *u = &fp->modules[net->u];
        const Module *v = &fp->modules[net->v];
        int cx_u = u->x + u->w / 2;
        int cy_u = u->y + u->h / 2;
        int cx_v = v->x + v->w / 2;
        int cy_v = v->y + v->h / 2;
        int dx = abs(cx_u - cx_v);
        int dy = abs(cy_u - cy_v);
        total += (dx + dy) * net->weight;
    }
    return total;
}

int write_output(const char *filename, const Floorplan *fp) {
    FILE *out = stdout;
    if (filename) {
        out = fopen(filename, "w");
        if (!out) {
            fprintf(stderr, "Error: cannot open output file %s\n", filename);
            return -1;
        }
    }
    int wl = compute_wirelength(fp);
    fprintf(out, "Wirelength %d\n", wl);
    fprintf(out, "NumSoftModules %d\n", fp->nsoft);
    for (int i = 0; i < fp->nsoft; i++) {
        const Module *m = &fp->modules[i];
        fprintf(out, "%s %d %d %d %d\n", m->name, m->x, m->y, m->w, m->h);
    }
    if (out != stdout) fclose(out);
    return 0;
}

void free_floorplan(Floorplan *fp) {
    if (!fp) return;
    for (int i = 0; i < fp->nModules; i++) {
        free(fp->modules[i].name);
    }
    free(fp->modules);
    free(fp->nets);
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        fprintf(stderr, "Usage: %s <input.txt> <output.floorplan>\n", argv[0]);
        return 1;
    }
    struct timeval tv_start, tv_end;
    gettimeofday(&tv_start, NULL);

    Floorplan fpData = {0};
    if (parse_input(argv[1], &fpData) != 0) return 1;
    compute_floorplan(&fpData);

    int hpwl = compute_wirelength(&fpData);
    if (write_output(argv[2], &fpData) != 0) {
        free_floorplan(&fpData);
        return 1;
    }

    printf("HPWL: %d\n", hpwl);
    // Check for any overlapping modules (pairwise)
    bool has_overlap = false;
    for (int i = 0; i < fpData.nModules && !has_overlap; i++) {
        for (int j = i + 1; j < fpData.nModules; j++) {
            Module *mi = &fpData.modules[i];
            Module *mj = &fpData.modules[j];
            if (mi->x < mj->x + mj->w && mj->x < mi->x + mi->w &&
                mi->y < mj->y + mj->h && mj->y < mi->y + mi->h) {
                has_overlap = true;
                break;
            }
        }
    }
    printf("Overlap: %s\n", has_overlap ? "YES" : "NO");

    gettimeofday(&tv_end, NULL);
    double elapsed = (tv_end.tv_sec - tv_start.tv_sec) +
                     (tv_end.tv_usec - tv_start.tv_usec) / 1e6;
    printf("Total execution time: %.6f seconds\n", elapsed);

    free_floorplan(&fpData);
    return 0;
}
