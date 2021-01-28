#ifndef CUBOID_C_H
#define CUBOID_C_H

/* C wrapper by CZhao */
#ifdef __cplusplus
extern "C" {
#endif

  extern void* cuboid_init(int *u);

  extern void cuboid_destroy(void *);

  extern void cuboid_transform(void *c, double x1, double x2,
    double x3, double *r1, double *r2, double *r3);

#ifdef __cplusplus
}
#endif

#endif // CUBOID_C_H
