#ifndef SLAVE_CERTIFICATE_H_
#define SLAVE_CERTIFICATE_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct t_filename;
struct t_filename {
	int x;
	char filename[FILENAME_LENGTH];
	struct t_filename * next;
};

#define	CERTIF_DOTPRODUCT	0x1200
#define CERTIF_POLYNOMIAL	0x1201

#define PINFO_HOST_LEN		40
#define PINFO_VERSION_LEN	120
#define PINFO_DATE_LEN		40

struct pinfo {
	char host[PINFO_HOST_LEN];
	char version[PINFO_VERSION_LEN];
	char date[PINFO_DATE_LEN];
};

struct certificate {
	int type;
	int vecnum;
	int r0;
	int r1;
	struct t_filename     * dotproducts;
	struct t_filename     * coeffs;
	struct t_filename     * vectors;
	struct t_filename     * sums;
	struct pinfo	      * producer;
	FILE * f;
};

extern int pinfo_date(struct pinfo *, time_t *);
extern int pinfo_fill(struct pinfo *);
extern struct pinfo * read_pinfo(const char *);
extern int write_pinfo(FILE *, struct pinfo *, char *);

extern int t_filename_head_insert(struct t_filename **, struct t_filename *);
extern int t_filename_tail_insert(struct t_filename **, struct t_filename *);
extern struct t_filename * new_t_filename(int, const char *);
extern struct t_filename * t_filename_head_pop(struct t_filename **pos);
extern struct t_filename * t_filename_lookup_pop(struct t_filename **, int);
extern struct certificate * new_certificate(int, int, int);
extern struct certificate * read_certificate(const char *);
extern int write_certificate(struct certificate *);
extern int sign_certificate(struct certificate *);
extern int close_certificate(struct certificate *);

#ifdef	__cplusplus
}
#endif

#endif	/* SLAVE_CERTIFICATE_H_ */
