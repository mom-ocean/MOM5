#ifndef SCATTERDIM_H
#define SCATTERDIM_H

#include "common.h"

typedef enum {NOSCATTER, SCATTERX, SCATTERY} scatter_t;

typedef struct ScatterDimStruct {
  int id;
  size_t len;
  char name[NC_MAX_NAME];
  
  scatter_t scatter_type;
  size_t* scatter_start;
  size_t* scatter_end;
  size_t* scatter_len;
  size_t scatter_ndiv; /* Number of divisions. */
} ScatterDim;

void ScatterDim_free(ScatterDim* p);

ScatterDim* ScatterDim_new(int id, size_t len, const char * name, scatter_t scatter_type, int ndiv);

#endif /* SCATTERDIM_H */