 #ifndef FLOORPLAN_H
 #define FLOORPLAN_H

 #include <stdio.h>

 typedef struct {
     char *name;
     int min_area;
     int isFixed;
     int x, y;
     int w, h;
 } Module;

 typedef struct {
     int u, v;
     int weight;
 } Net;

 typedef struct {
     int chipWidth;
     int chipHeight;
     int nModules;
     Module *modules;
     int nsoft;
     int nf;
     int nNets;
     Net *nets;
 } Floorplan;

/**
 * Parse input file into Floorplan struct.
 * Returns 0 on success, non-zero on error.
 */
int parse_input(const char *filename, Floorplan *fp);

/**
 * Compute a simple floorplan for soft modules.
 */
void compute_floorplan(Floorplan *fp);

/**
 * Compute total weighted half-perimeter wirelength (HPWL).
 */
int compute_wirelength(const Floorplan *fp);

/**
 * Write floorplan to output file (or stdout if filename is NULL).
 * Returns 0 on success.
 */
int write_output(const char *filename, const Floorplan *fp);

/**
 * Free memory allocated in Floorplan.
 */
void free_floorplan(Floorplan *fp);

#endif /* FLOORPLAN_H */