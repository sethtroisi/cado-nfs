#ifndef RELATION_H_
#define RELATION_H_

#ifdef __cplusplus
extern "C" {
#endif

// Relation I/O
extern void copy_rel(relation_t *Rel, relation_t rel);
extern void clear_relation(relation_t *rel);
extern int read_relation(relation_t *rel, const char *str);
extern int fread_relation(FILE *file, relation_t *rel);
extern unsigned long findroot(long a, unsigned long b, unsigned long p);
extern void computeroots(relation_t * rel);
extern void fprint_relation(FILE *file, relation_t rel);
extern void reduce_exponents_mod2 (relation_t *rel);


#ifdef __cplusplus
}
#endif

#endif	/* RELATION_H_ */
