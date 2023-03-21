#ifndef CONFIG_H
#define CONFIG_H

/**
Set the compile configuration. Don't forget to set the according constants in the file "constants".
**/

// debug mode on/off
#define NDEBUG

#undef CONSTANTS_DEBUG
#undef SWEEPLINE_DEBUG
#undef RED_BLACK_DEBUG
#undef CORE_DEBUG
#undef GRAPH_NUMBERS_DEBUG
#undef BINARY_INPUT_OUTPUT_DEBUG
#undef BOX_DEBUG
#undef NEIGHBOUR_DEBUG
#undef POINT_LOCATION_DEBUG
#undef VARIABLES_DEBUG

#undef FUNCTION_POINTER_DEBUG
#undef RIEMANN_SOLVER_DEBUG

#undef PARTICLE_DEBUG
#undef HALF_EDGE_DEBUG

#undef AREA_DEBUG_PARTICLE
#undef AREA_DEBUG_HALF_EDGE

// debug area
#define AREA_XMIN    -1.25
#define AREA_YMIN    7.85
#define AREA_XMAX    -0.8
#define AREA_YMAX    8.15

// sweepline algorithm
#undef SEARCH_ARC_2

// fix same point degeneracy
#undef SPD_FIX_IN_PROCESS // works only for two points

#endif
