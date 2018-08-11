#ifndef GPF_H_
#define GPF_H_

extern unsigned int *gpf;

void gpf_init(unsigned int);
static inline unsigned int gpf_get(const unsigned long i) {
    return gpf[i];
}

void gpf_clear();

#endif