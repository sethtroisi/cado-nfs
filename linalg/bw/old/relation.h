#ifndef RELATION_H_
#define RELATION_H_

#ifdef	__cplusplus
extern "C" {
#endif

/**/
typedef struct {
    type32  index;
    stype32 expo;
} sp_data_s;

typedef sp_data_s * sp_data;

/**/
typedef struct {
    type32  poly[2];
    stype32 expo;
} lp_data_s;

typedef lp_data_s * lp_data;

/**/
typedef struct {
    type32 info;
    type32 * data;
} relation_s;

typedef relation_s * relation;

/**/
typedef struct {
    type32 info;
    type32 * data;
    short nSPs;
    short nLPs;
    sp_data SP;
    lp_data LP;
    unsigned char degLP[2];
    type32 * checksum;
    size_t data_length;
} ext_relation_s;

typedef ext_relation_s * ext_relation;

/************************************************************************/

void    _relation_alloc         (relation *, type32);
void    _relation_free          (relation *);
int      relation_read          (relation,   FILE *);
int      relation_write         (FILE *,     relation);
#ifdef CHECKSUM_ENABLED
void     relation_checksum      (type32 *,   const type32 *,   relation);
#endif
void    _relation_extend        (ext_relation *, relation *);
unsigned char relation_degLP    (type32 *);

#define relation_alloc(x,l)     _relation_alloc(&(x),l)
#define relation_free(x)        _relation_free(&(x))
#define relation_extend(y,x)    _relation_extend(&(y),&(x))

#define RELATION_W_CHECKSUM      UC32(0x20000000U)
#define RELATION_HAS_CHECKSUM(x) ((((x)->info)&RELATION_W_CHECKSUM)!=0)
#define RELATION_LPS_TYPE(x)     (((x)->info)&UC32(0xC0000000U))
#define RELATION_FF              UC32(0)
#define RELATION_PF_FP           UC32(0x40000000U)
#define RELATION_PP              UC32(0xC0000000U)
#define RELATION_NSPS(x)         (((x)->info)&UC32(0x1FFFFFFFU))
#define RELATION_NLPS(x) ((((x)->info)>>30)^((((x)->info)&UC32(0x80000000U))!=UC32(0)))
#define RELATION_CS(x)  ((type32*) (RELATION_HAS_CHECKSUM(x)?((x)->data):NULL))
#define RELATION_LP(x)  ((lp_data) (RELATION_LPS_TYPE(x)?(((x)->data)+(RELATION_HAS_CHECKSUM(x)?2:0)):NULL))
#define RELATION_SP(x)  ((sp_data) (((x)->data)+(RELATION_HAS_CHECKSUM(x)?2:0)+(RELATION_NLPS(x)*3)))
#define RELATION_DATASIZE(x) (RELATION_HAS_CHECKSUM(x)?2:0)+(RELATION_NLPS(x)*3)+(RELATION_NSPS(x)<<1)
#define RELATION_WITH(x)        ((x)==0?RELATION_FF:((x)==1?RELATION_PF_FP:RELATION_PP))

void    _ext_relation_alloc     (ext_relation *, type32);
void    _ext_relation_free      (ext_relation *);
int      ext_relation_read      (ext_relation,   FILE *);
int      ext_relation_write     (FILE *,     ext_relation);
#ifdef CHECKSUM_ENABLED
void     ext_relation_checksum  (type32 *,   const type32 *, ext_relation);
#endif
void    _ext_relation_dup       (ext_relation *, ext_relation);
void     ext_relation_reap_checksum (ext_relation);
void     ext_relation_adjust    (ext_relation);

#define ext_relation_alloc(x,l) _ext_relation_alloc(&(x),(type32)(l))
#define ext_relation_free(x)    _ext_relation_free(&(x))
#define ext_relation_dup(y,x)   _ext_relation_dup(&(y),x)

int sp_compare(sp_data, sp_data);

#ifdef	__cplusplus
}
#endif

/* vim:set sw=8: */
#endif	/* RELATION_H_ */
