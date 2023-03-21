#ifndef DEBUG_H
#define DEBUG_H

#include "config.h"

#ifdef SWEEPLINE_DEBUG
#define SD_PRINTF(...) printf(__VA_ARGS__)
#else
#define SD_PRINTF(...)
#endif

#ifdef RED_BLACK_DEBUG
#define RB_PRINTF(...) printf(__VA_ARGS__)
#else
#define RB_PRINTF(...)
#endif

#ifdef CORE_DEBUG
#define CO_PRINTF(...) printf(__VA_ARGS__)
#else
#define CO_PRINTF(...)
#endif

#ifdef GRAPH_NUMBERS_DEBUG
#define GN_PRINTF(...) printf(__VA_ARGS__)
#else
#define GN_PRINTF(...)
#endif

#ifdef BINARY_INPUT_OUTPUT_DEBUG
#define IO_PRINTF(...) printf(__VA_ARGS__)
#else
#define IO_PRINTF(...)
#endif

#ifdef BOX_DEBUG
#define BOX_PRINTF(...) printf(__VA_ARGS__)
#else
#define BOX_PRINTF(...)
#endif

#ifdef NEIGHBOUR_DEBUG
#define NE_PRINTF(...) printf(__VA_ARGS__)
#else
#define NE_PRINTF(...)
#endif

#ifdef VARIABLES_DEBUG
#define VAR_PRINTF(...) printf(__VA_ARGS__)
#else
#define VAR_PRINTF(...)
#endif

#ifdef GRADIENTS_DEBUG
#define GRAD_PRINTF(...) printf(__VA_ARGS__)
#else
#define GRAD_PRINTF(...)
#endif

#ifdef FUNCTION_POINTER_DEBUG
#define FP_PRINTF(...) printf(__VA_ARGS__)
#else
#define FP_PRINTF(...)
#endif

#ifdef RIEMANN_SOLVER_DEBUG
#define RS_PRINTF(...) printf(__VA_ARGS__)
#else
#define RS_PRINTF(...)
#endif

#ifdef FLUX_DEBUG
#define FLUX_PRINTF(...) printf(__VA_ARGS__)
#else
#define FLUX_PRINTF(...)
#endif

#ifdef PARTICLE_DEBUG
#define VP_PRINTF(...) printf(__VA_ARGS__)
#else
#define VP_PRINTF(...)
#endif

#ifdef HALF_EDGE_DEBUG
#define HE_PRINTF(...) printf(__VA_ARGS__)
#else
#define HE_PRINTF(...)
#endif

#ifdef POINT_LOCATION_DEBUG
#define PL_PRINTF(...) printf(__VA_ARGS__)
#else
#define PL_PRINTF(...)
#endif

#ifdef AREA_DEBUG_PARTICLE
#undef VP_PRINTF
#define VP_PRINTF(...) if(m_location.x <= AREA_XMAX && m_location.x >= AREA_XMIN && m_location.y <= AREA_YMAX && m_location.y >= AREA_YMIN)printf(__VA_ARGS__)
#endif

#ifdef AREA_DEBUG_HALF_EDGE
#undef HE_PRINTF
#define HE_PRINTF(...) if(m_center.x <= AREA_XMAX && m_center.x >= AREA_XMIN && m_center.y <= AREA_YMAX && m_center.y >= AREA_YMIN)printf(__VA_ARGS__)
#endif

#endif