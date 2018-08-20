#ifndef HWLOC_AUX_H_
#define HWLOC_AUX_H_

#include <hwloc.h>

#ifdef __cplusplus
extern "C" {
#endif


int hwloc_aux_get_depth_from_string(hwloc_topology_t topology, const char * desc);

#ifdef __cplusplus
}
#endif

#endif	/* HWLOC_AUX_H_ */
