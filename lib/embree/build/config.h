/* #undef EMBREE_RAY_MASK */
/* #undef EMBREE_STAT_COUNTERS */
/* #undef EMBREE_BACKFACE_CULLING */
#define EMBREE_INTERSECTION_FILTER
#define EMBREE_INTERSECTION_FILTER_RESTORE
/* #undef EMBREE_BUFFER_STRIDE */
/* #undef EMBREE_RETURN_SUBDIV_NORMAL */
/* #undef EMBREE_IGNORE_INVALID_RAYS */
#define EMBREE_GEOMETRY_TRIANGLES
#define EMBREE_GEOMETRY_QUADS
#define EMBREE_GEOMETRY_LINES
#define EMBREE_GEOMETRY_HAIR
#define EMBREE_GEOMETRY_SUBDIV
#define EMBREE_GEOMETRY_USER
#define EMBREE_RAY_PACKETS

#if defined(EMBREE_GEOMETRY_TRIANGLES)
  #define IF_ENABLED_TRIS(x) x
#else
  #define IF_ENABLED_TRIS(x)
#endif

#if defined(EMBREE_GEOMETRY_QUADS)
  #define IF_ENABLED_QUADS(x) x
#else
  #define IF_ENABLED_QUADS(x)
#endif

#if defined(EMBREE_GEOMETRY_LINES)
  #define IF_ENABLED_LINES(x) x
#else
  #define IF_ENABLED_LINES(x)
#endif

#if defined(EMBREE_GEOMETRY_HAIR)
  #define IF_ENABLED_HAIR(x) x
#else
  #define IF_ENABLED_HAIR(x)
#endif

#if defined(EMBREE_GEOMETRY_SUBDIV)
  #define IF_ENABLED_SUBDIV(x) x
#else
  #define IF_ENABLED_SUBDIV(x)
#endif

#if defined(EMBREE_GEOMETRY_USER)
  #define IF_ENABLED_USER(x) x
#else
  #define IF_ENABLED_USER(x)
#endif




