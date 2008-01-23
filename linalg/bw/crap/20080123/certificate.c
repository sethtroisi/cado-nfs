#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/utsname.h>
#include <ctype.h>
#include <errno.h>
#include <dirent.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#ifdef ENABLE_PTHREADS
#include <pthread.h>
#endif

#include "slave_params.h"
#include "params.h"
#include <assert.h>
#include "macros.h"
#include "auxfuncs.h"
#include "types.h"
#include "macros.h"
#include "filenames.h"
#include "bw_scalar.h"
#include "bw_lvblock.h"
#include "options.h"
#include "corners.h"
#include "slave.h"
#include "variables.h"
#include "bw_lvblock_steps.h"
#include "threaded.h"
#include "barrier.h"
/* #include "version.h" */
#include "certificate.h"

const char * ctname_polynomial="POLYNOMIAL\n";
const char * ctname_dotproduct="DOTPRODUCT\n";

#define CT_NTAGS	13

#define CT_FILE		0
#define CT_VERSION	1
#define CT_HOST		2
#define CT_DATE		3
#define CT_VECNUM	4
#define CT_BEGIN	5
#define CT_END		6
#define CT_TYPE		7
#define CT_DOTPRODS	8
#define CT_COEFFS	9
#define CT_VECTORS	10
#define CT_SUMS		11
#define CT_SIGNATURE	12

const char ctag_table[CT_NTAGS][12]	= {
	/* CT_FILE */		"FILE: ",
	/* CT_VERSION */	"VERSION: ",
	/* CT_HOST */		"HOST: ",
	/* CT_DATE */		"DATE: ",
	/* CT_VECNUM */		"VECNUM: ",
	/* CT_BEGIN */		"BEGIN: ",
	/* CT_END */		"END: ",
	/* CT_TYPE */		"TYPE: ",
	/* CT_DOTPRODS */	"DOTPRODS: ",
	/* CT_COEFFS */		"COEFFS: ",
	/* CT_VECTORS */	"VECTORS: ",
	/* CT_SUMS */		"SUMS: ",
	/* CT_SIGNATURE */	"SIG: "
};

int t_filename_head_insert(	struct t_filename ** pos,
				struct t_filename * newdata)
{
	newdata->next=*pos;
	*pos=newdata;
	return 0;
}

int t_filename_tail_insert(	struct t_filename ** pos,
					struct t_filename * newdata)
{
	int i;
	for(i=0;*pos;pos=&((*pos)->next),i++);
	return i+t_filename_head_insert(pos,newdata);
}

struct t_filename * t_filename_head_pop(struct t_filename **pos)
{
	struct t_filename * res;
	res=*pos;
	*pos=res->next;
	res->next=NULL;
	return res;
}

struct t_filename * t_filename_lookup_pop(struct t_filename ** pos,
		int x)
{
	for(;*pos && (*pos)->x !=x ;pos=&((*pos)->next));
	if (*pos==NULL) {
		return NULL;
	}
	return t_filename_head_pop(pos);
}

struct t_filename * new_t_filename(int tag, const char * name)
{
	struct t_filename * res;

	res=malloc(sizeof(struct t_filename));
	res->x=tag;
	mkfname(res->filename,"%s",name);
	res->next=NULL;	/* That's not necessary, but it's probably wiser to
		   	 * do so anyway */
	return res;
}

int close_certificate(struct certificate * certif)
{
	for(;certif->dotproducts!=NULL;)
		free(t_filename_head_pop(&(certif->dotproducts)));
	for(;certif->coeffs!=NULL;)
		free(t_filename_head_pop(&(certif->coeffs)));
	for(;certif->vectors!=NULL;)
		free(t_filename_head_pop(&(certif->vectors)));
	for(;certif->sums!=NULL;)
		free(t_filename_head_pop(&(certif->sums)));
	if (certif->producer!=NULL)
		free(certif->producer);
	if (certif->f!=NULL)
		fclose(certif->f);
	free(certif);
	return 0;
}

struct certificate * new_certificate(int type, int j, int r0)
{
	struct certificate * res;

	res=malloc(sizeof(struct certificate));
	res->type=type;
	res->vecnum=j;
	res->r0 = r0;
	res->r1 = -1;	/* because it is not yet specified */
	res->dotproducts=NULL;
	res->coeffs=NULL;
	res->vectors=NULL;
	res->sums=NULL;
	res->f=NULL;
	res->producer=NULL;
	return res;
}

static int chop_string(char * s, int maxlen)
{
	s[maxlen-1]='\0';
	maxlen=strlen(s);
	for(;maxlen>=0 && isspace(s[maxlen-1]);s[--maxlen]='\0');
	return 0;
}

static int pinfo_secure(struct pinfo * dst)
{
	chop_string(dst->host,PINFO_HOST_LEN);
	chop_string(dst->version,PINFO_VERSION_LEN);
	chop_string(dst->date,PINFO_DATE_LEN);

	if (strlen(dst->host)==0) {
		strcpy(dst->host,"<unknown>");
	}
	if (strlen(dst->version)==0) {
		strcpy(dst->version,"<unknown>");
	}
	if (strlen(dst->date)==0) {
		strcpy(dst->date,"<unknown>");
	}

	return 0;
}

int pinfo_fill(struct pinfo * dst)
{
	struct utsname	buf;
	time_t	tt;

	uname(&buf);
	strcpy(dst->host,buf.nodename);
	time(&tt);
#ifdef __SVR4
	ctime_r(&tt,dst->date,PINFO_DATE_LEN);
#else
	ctime_r(&tt,dst->date);
#endif
	/* strcpy(dst->version,version_string); */

	pinfo_secure(dst);

	return 0;
}

int pinfo_date(struct pinfo * dst, time_t * tp)
{
#ifdef __SVR4
	ctime_r(tp,dst->date,PINFO_DATE_LEN);
#else
	ctime_r(tp,dst->date);
#endif
	chop_string(dst->date,PINFO_DATE_LEN);

	return 0;
}


int sign_certificate(struct certificate * certif)
{
	struct pinfo verifier;

	pinfo_fill(&verifier);

	fprintf(certif->f,"%s *** NEW SIGNATURE BLOCK ***\n",
			ctag_table[CT_SIGNATURE]);
	fprintf(certif->f,"%sChecked by %s on %s\n",
			ctag_table[CT_SIGNATURE],
			verifier.host,verifier.date);
	fprintf(certif->f,"%sChecker version was %s",
			ctag_table[CT_SIGNATURE],
			verifier.version);
	return 0;
}

int write_pinfo(FILE *f, struct pinfo * prod, char * hdr)
{
	fprintf(f,"%s%s%s\n",hdr,ctag_table[CT_HOST],prod->host);
	fprintf(f,"%s%s%s\n",hdr,ctag_table[CT_VERSION],prod->version);
	fprintf(f,"%s%s%s\n",hdr,ctag_table[CT_DATE],prod->date);

	return 0;
}

static int parse_pinfo_tagged_line(struct pinfo * dst, int tag, const char * data)
{
	switch(tag) {
		case CT_VERSION:
			sprintf(dst->version,data);
			break;
		case CT_HOST:
			sprintf(dst->host,data);
			break;
		case CT_DATE:
			sprintf(dst->date,data);
			break;
		default:
			break;
	}

	return 0;
}

static int parse_tagged_line(struct certificate * certif, int tag, const char * data)
{
	int n;
	char *t;

	switch(tag) {
		case CT_FILE:	/* do nothing */
		case CT_SIGNATURE:
			break;
		case CT_VERSION:
		case CT_HOST:
		case CT_DATE:
			parse_pinfo_tagged_line(certif->producer,tag,data);
			break;
		case CT_VECNUM:
			certif->vecnum=strtol(data,&t,0);
			if (*t!='\n' && *t !='\0')
				goto parse_error;
			break;
		case CT_BEGIN:
			certif->r0=strtol(data,&t,0);
			if (*t!='\n' && *t !='\0')
				goto parse_error;
			break;
		case CT_END:
			certif->r1=strtol(data,&t,0);
			if (*t!='\n' && *t !='\0')
				goto parse_error;
			break;
		case CT_TYPE:
			if (strcmp(data,ctname_polynomial)==0) {
				certif->type=CERTIF_POLYNOMIAL;
			} else if (strcmp(data,ctname_dotproduct)==0) {
				certif->type=CERTIF_DOTPRODUCT;
			} else {
				if (*data!='\n' && *data !='\0')
					goto parse_error;
			}
			break;
		case CT_DOTPRODS:
			n=strtol(data,&t,0);
			if (!isspace(*t))
				goto parse_error;
			for(;isspace(*t);t++);
			for(;isspace(t[strlen(t)-1]);t[strlen(t)-1]=0);
			t_filename_tail_insert(&(certif->dotproducts),
					new_t_filename(n,t));
			break;
		case CT_COEFFS:
			n=strtol(data,&t,0);
			if (!isspace(*t))
				goto parse_error;
			for(;isspace(*t);t++);
			for(;isspace(t[strlen(t)-1]);t[strlen(t)-1]=0);
			t_filename_tail_insert(&(certif->coeffs),
					new_t_filename(n,t));
			break;
		case CT_VECTORS:
			n=strtol(data,&t,0);
			if (!isspace(*t))
				goto parse_error;
			for(;isspace(*t);t++);
			for(;isspace(t[strlen(t)-1]);t[strlen(t)-1]=0);
			t_filename_tail_insert(&(certif->vectors),
					new_t_filename(n,t));
		case CT_SUMS:
			n=strtol(data,&t,0);
			if (!isspace(*t))
				goto parse_error;
			for(;isspace(*t);t++);
			for(;isspace(t[strlen(t)-1]);t[strlen(t)-1]=0);
			t_filename_tail_insert(&(certif->sums),
					new_t_filename(n,t));
		default:
			break;
	}
	return 0;
parse_error:
	fprintf(stderr,"Parse error while parsing certificate file\n"
			"Failed before: %s",t);
	return -1;
}


static int read_one_line(char **s, char * buf, FILE * f)
{
	char * res;
	int i;

	for(;;) {
		res=fgets(buf,LINEBUF_LENGTH,f);

		if (res==NULL) return -1;
		if (buf[0]!='#' && buf[0]!='\n') break;
	}

	for(i=0;i<CT_NTAGS;i++) {
		if (strncmp(buf,ctag_table[i],strlen(ctag_table[i]))==0)
			break;
	}

	if (i==CT_NTAGS)
		return -1;

	*s=buf+strlen(ctag_table[i]);

	return i;
}

static int check_has(const char * file, const void * p, int field)
{
	if (!p) {
		fprintf(stderr,"File %s lacks a ``%s'' field\n",
				file,ctag_table[field]);
		exit(1);
	}
	return 0;
}

static int check_hasnot(const char * file, const void * p, int field)
{
	if (p) {
		fprintf(stderr,"File %s has an unwanted ``%s'' field\n",
				file,ctag_table[field]);
		exit(1);
	}
	return 0;
}

struct pinfo * read_pinfo(const char * filename)
{
	struct pinfo  * res;
	FILE  * f;
	int	tag;
	char  * s;
	char	linebuf[LINEBUF_LENGTH];

	f=fopen(filename,"r");
	res=malloc(sizeof(struct pinfo));

	for(;;) {
		tag=read_one_line(&s,linebuf,f);
		if (tag==-1)
			break;
		parse_pinfo_tagged_line(res,tag,s);	/* Never fails */
	}

	fclose(f);

	pinfo_secure(res);

	return res;
}

struct certificate * read_certificate(const char * filename)
{
	struct certificate * res;
	
	char	linebuf[LINEBUF_LENGTH];
	FILE  * f;
	int	tag;
	char  * s;
	long	offset;
	
	f=fopen(filename,"r");
	res=malloc(sizeof(struct certificate));
	res->type=-1;
	res->vecnum=-1;
	res->r0 = -1;
	res->r1 = -1;
	res->dotproducts = NULL;
	res->coeffs = NULL;
	res->vectors = NULL;
	res->producer=malloc(sizeof(struct pinfo));

	for(;;) {
		tag=read_one_line(&s,linebuf,f);
		if (tag==-1)
			break;
		if (parse_tagged_line(res,tag,s)<0) {	/* parse error */
			fprintf(stderr,"Failed file was: %s\n",filename);
			exit(1);
		}
	}

	if (!(res->r0 < res->r1)) {
		fprintf(stderr,"File %s : wrong ordering\n",filename);
		exit(1);
	}

	if ((res->vecnum < 0 || res->vecnum > nbys)) {
		fprintf(stderr,"File %s : wrong ordering\n",filename);
		exit(1);
	}

	switch(res->type) {
		case CERTIF_DOTPRODUCT:
			check_has(filename,res->dotproducts,CT_DOTPRODS);
			check_has(filename,res->vectors,CT_VECTORS);
			check_hasnot(filename,res->coeffs,CT_COEFFS);
			break;
		case CERTIF_POLYNOMIAL:
			check_hasnot(filename,res->dotproducts,CT_DOTPRODS);
			check_has(filename,res->vectors,CT_VECTORS);
			check_has(filename,res->coeffs,CT_COEFFS);
			check_has(filename,res->vectors,CT_SUMS);
			break;
		default:
			fprintf(stderr,"File %s : unknown type\n",filename);
			exit(1);
	}
	
	offset=ftell(f);
	fclose(f);
	res->f=fopen(filename,"a");
	if (ftell(res->f)!=offset) {
		fprintf(stderr,"File %s : inconsistent reopen\n",filename);
		exit(1);
	}

	return res;
}


int write_certificate(struct certificate * certif)
{
	struct t_filename     * llp;
	struct pinfo		my_id;

	FILE  * f;
	char	filename[FILENAME_LENGTH];

	assert(certif->f==NULL);	/* Otherwise, we can only append... */
		
	mkfname(filename,
			certif_meta_filename,
			certif->vecnum,
			certif->r0,
			certif->r1);
	if (exist(filename)) {
		fprintf(stderr,"Warning: Certificate %s already exists\n",
				filename);
	}


	f=fopen(filename,"w");

	fprintf(f,"%s%s\n",ctag_table[CT_FILE],filename);
	if (certif->producer==NULL) {
		pinfo_fill(&my_id);
		write_pinfo(f,&my_id,"");
	} else {
		write_pinfo(f,certif->producer,"");
	}
	fprintf(f,"%s%d\n",ctag_table[CT_VECNUM],certif->vecnum);
	fprintf(f,"%s%d\n",ctag_table[CT_BEGIN],certif->r0);
	fprintf(f,"%s%d\n",ctag_table[CT_END],certif->r1);
	fputs(ctag_table[CT_TYPE],f);
	switch(certif->type) {
		case CERTIF_POLYNOMIAL:
			fputs(ctname_polynomial,f);
			break;
		case CERTIF_DOTPRODUCT:
			fputs(ctname_dotproduct,f);
			break;
		default:
			assert(0);
	}
	for(llp=certif->dotproducts;llp;llp=llp->next) {
		fprintf(f,"%s%d %s\n",
				ctag_table[CT_DOTPRODS],
				llp->x,llp->filename);
	}
	for(llp=certif->coeffs;llp;llp=llp->next) {
		fprintf(f,"%s%d %s\n",
				ctag_table[CT_COEFFS],
				llp->x,llp->filename);
	}
	for(llp=certif->vectors;llp;llp=llp->next) {
		fprintf(f,"%s%d %s\n",
				ctag_table[CT_VECTORS],
				llp->x,llp->filename);
	}
	for(llp=certif->sums;llp;llp=llp->next) {
		fprintf(f,"%s%d %s\n",
				ctag_table[CT_SUMS],
				llp->x,llp->filename);
	}
	certif->f=f;
		
	printf("ISSUED CERTIFICATE: %s\n",filename);
	return 0;
}
