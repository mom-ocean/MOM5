#include "scatterdim.h"

ScatterDim* ScatterDim_new(int id, size_t len, const char * name, scatter_t scatter_type, int ndiv) {
  ScatterDim* p = XMALLOC(ScatterDim,1);
  if (p==NULL) return p;
  
  p->id = id;
  p->len = len;
  strcpy(p->name, name);
  p->scatter_type = scatter_type;
  p->scatter_ndiv = ndiv;

  p->scatter_start = XMALLOC(size_t, ndiv);
  p->scatter_end = XMALLOC(size_t, ndiv);
  p->scatter_len = XMALLOC(size_t, ndiv);
  
  return p;
}
/*-------------------------------------------------------------------*/
void ScatterDim_free(ScatterDim* p) {
  if (p == NULL) return;
  
  if (p->scatter_start) {
    XFREE(p->scatter_start);
  }

  if (p->scatter_end) {
    XFREE(p->scatter_end);
  }

  if (p->scatter_len) {
    XFREE(p->scatter_len);
  }  
}
