#include <sys/stat.h>
#include <sys/time.h>
#include <sys/types.h>
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
#include "recovery.h"

/* Variables globales utilisées:
 *
 * nbxs
 * ncols
 * bw_filesize
 * keep_periodicity
 * vs_checkpoint_barrier
 * task_type
 */


int keep_periodicity = DEFAULT_KEEP;
barrier_t vs_checkpoint_barrier;


/************************************************************************/

/* lookup_a_files
 * 
 * Cherche les fichiers A pour une colonne donnée.
 *
 * Erreurs si:	deux révisions concurrentes existent.
 * 		un fichier est absent.
 * 		un fichier a une taille incohérente avec son indice.
 *
 * en cas d'erreur, arrêt du prg.
 *
 * RETURN: 0 si aucun fichier n'est présent (message) ; sinon, la
 * révision(commune), qui peut être 0, d'ailleurs.
 */

static int lookup_a_files(DIR * dirp, int vecnum)
{
	struct dirent * curr;
	struct stat	sbuf;

	int   * rev;
	char  * pattern;
	int	res;
	int	i,j,stamp;
	char	filename[FILENAME_LENGTH];

	pattern=strrchr(a_meta_filename,'/');
	if (pattern==NULL) pattern=a_meta_filename; else pattern++;

	rev=malloc(nbxs*sizeof(int));
	for(i=0;i<nbxs;rev[i++]=-1);

	rewinddir(dirp);
	for(;(curr=readdir(dirp))!=NULL;) {
		res=sscanf(curr->d_name,pattern,&i,&j,&stamp);
		if (res!=3 || j!=vecnum)	continue;
		mkfname(filename,pattern,i,j,stamp);

		/* The following comparison might fail in case the parsed
		 * file has a suffix. We want to discard it in this case.
		 */
		if (strcmp(filename,curr->d_name)!=0)
			continue;

		printf("Found %s\n",curr->d_name);
		if (rev[i]>=0) {
			fprintf(stderr,a_meta_filename,i,j,rev[i]);
			fprintf(stderr," clashes with ");
			fprintf(stderr,a_meta_filename,i,j,stamp);
			exit(1);
		}
		rev[i]=stamp;
		mkfname(filename,a_meta_filename,i,j,stamp);
		if (stat(filename,&sbuf)<0) {
			perror(filename);
			exit(1);
		}
		if (sbuf.st_size/STAMPSIZE!=stamp) {
			fprintf(stderr,"%s has wrong size\n",filename);
			exit(1);
		}
	}
	
	stamp=rev[0];
	for(i=1;i<nbxs;i++) {
		if (rev[i]==stamp) continue;
		fprintf(stderr,"There is no A<%02d,%02d> file\n",i,vecnum);
		exit(1);
	}

	if (stamp==-1) {
		printf("Found no A file\n");
		stamp=0;
	}

	free(rev);

	return stamp;
}

/* lookup_v_files
 *
 * Cherche les fichiers V pour une colonne donnée.
 *
 * Erreurs si:	deux révisions concurrentes existent.
 * 		le fichier n'a pas la bonne taille.
 *
 * RETURN: 0 si aucun fichier n'est présent (message) ; sinon, la
 * révision (qui ne peut pas être 0)
 *
 * Il est normalement impossible que cete fonction renvoie 0.
 */
static int lookup_v_files(DIR * dirp, int vecnum)
{
	struct dirent * curr;
	struct stat	sbuf;

	int	rev;
	char  * pattern;
	int	res;
	int	j,stamp;
	char	filename[FILENAME_LENGTH];

	pattern=strrchr(v_meta_filename,'/');
	if (pattern==NULL) pattern=v_meta_filename; else pattern++;

	rev=-1;

	rewinddir(dirp);
	for(;(curr=readdir(dirp))!=NULL;) {
		res=sscanf(curr->d_name,pattern,&j,&stamp);
		if (res!=2 || j!=vecnum)	continue;
		mkfname(filename,pattern,j,stamp);

		/* The following comparison might fail in case the parsed
		 * file has a suffix. We want to discard it in this case.
		 */
		if (strcmp(filename,curr->d_name)!=0)
			continue;

		printf("Found %s\n",curr->d_name);
		if (rev>=0) {
			fprintf(stderr,v_meta_filename,j,rev);
			fprintf(stderr," clashes with ");
			fprintf(stderr,v_meta_filename,j,stamp);
			exit(1);
		}
		rev=stamp;
		mkfname(filename,v_meta_filename,j,stamp);
		if (stat(filename,&sbuf)<0) {
			perror(filename);
			exit(1);
		}
		if (sbuf.st_size!=ncols*bw_filesize*sizeof(mp_limb_t)) {
			fprintf(stderr,"%s has wrong size\n",filename);
			exit(1);
		}
	}
	
	if (rev==-1) {
		printf("Found no V file\n");
		rev=0;
	}

	return rev;
}

/* lookup_s_files
 *
 * Cherche les fichiers S pour une colonne donnée.
 *
 * Erreurs si:	deux révisions concurrentes existent.
 * 		le fichier n'a pas la bonne taille.
 *
 * RETURN: 0 si aucun fichier n'est présent (message) ; sinon, la
 * révision (qui ne peut pas être 0)
 *
 * Il est normalement impossible que cete fonction renvoie 0. En effet,
 * S0 est toujours le vecteur nul.
 */
static int lookup_s_files(DIR * dirp, int vecnum)
{
	struct dirent * curr;
	struct stat	sbuf;

	int	rev;
	char  * pattern;
	int	res;
	int	j,stamp;
	char	filename[FILENAME_LENGTH];

	pattern=strrchr(s_meta_filename,'/');
	if (pattern==NULL) pattern=s_meta_filename; else pattern++;

	rev=-1;

	rewinddir(dirp);
	for(;(curr=readdir(dirp))!=NULL;) {
		res=sscanf(curr->d_name,pattern,&j,&stamp);
		if (res!=2 || j!=vecnum)	continue;
		mkfname(filename,pattern,j,stamp);

		/* The following comparison might fail in case the parsed
		 * file has a suffix. We want to discard it in this case.
		 */
		if (strcmp(filename,curr->d_name)!=0)
			continue;

		printf("Found %s\n",curr->d_name);
		if (rev>=0) {
			fprintf(stderr,s_meta_filename,j,rev);
			fprintf(stderr," clashes with ");
			fprintf(stderr,s_meta_filename,j,stamp);
			exit(1);
		}
		rev=stamp;
		mkfname(filename,s_meta_filename,j,stamp);
		if (stat(filename,&sbuf)<0) {
			perror(filename);
			exit(1);
		}
		if (sbuf.st_size!=ncols*bw_filesize*sizeof(mp_limb_t)) {
			fprintf(stderr,"%s has wrong size\n",filename);
			exit(1);
		}
	}
	
	if (rev==-1) {
		printf("Found no S file\n");
		rev=0;
	}

	return rev;
}

/* lookup_trace_file
 *
 * Cherche les restes d'un run précédent.
 *
 * Erreurs si:	La révision atteinte n'est pas correcte.
 * 		Le cstart est > au stamp
 *
 * RETURN: -1 si aucun fichier n'est présent (message) ; 0 sinon
 *
 */
static int lookup_trace_file(struct vs_data * vsd)
{
	struct stat	sbuf;
	struct pinfo  * old_tenant;

	char	filename[FILENAME_LENGTH];
	char	linebuf[LINEBUF_LENGTH];
	FILE  * f;
	int	stamp,cstart;

	mkfname(filename,l_meta_filename,vsd->vecnum);

	if (!exist(filename))
		return -1;	/* no trace file */

	old_tenant=read_pinfo(filename);

	f=fopen(filename,"r");
	if (f==NULL) {
		perror(filename);	/* must be worse than ENOENT */
		exit(1);
	}


	for(;!feof(f);) {
		fgets(linebuf,1024,f);
		if (strncmp(linebuf,"ID_END",6)==0)
			break;
	}
	if (feof(f)) {
		fclose(f);
		fprintf(stderr,"Incomplete trace file\n");
		exit(1);	/* we could return -1 as well */
	}
	fscanf(f,"%*[^:]: %d\n",&stamp);
	fscanf(f,"%*[^:]: %d\n",&cstart);

	if (stamp != vsd->rev_v || cstart > vsd->rev_v) {
		fprintf(stderr,"File %s has incoherent data\n",filename);
		exit(1);	/* User repairs by himself */
	}

	/* Obtain the previous producer info */
	fstat(fileno(f),&sbuf);
	pinfo_date(old_tenant,&sbuf.st_mtime);
	vsd->tenant=old_tenant;
	fclose(f);

	return cstart;
}

static void write_trace_file(struct vs_data *vsd, int n)
{
	FILE  * f;
	char	filename[FILENAME_LENGTH];
	long	pos;

	if (!is_master_thread())
		return;

	mkfname(filename,l_meta_filename,vsd->vecnum);
	f=fopen(filename,"r+");
	fseek(f,vsd->trace_offset,SEEK_SET);
	fprintf(f,"Last stamp: %d\n",vsd->rev_v+n);
	fprintf(f,"Running certificate began at: %d\n",vsd->rc->r0);
	fflush(f);
	pos=ftell(f);
	fclose(f);
	truncate(filename,pos);
}

static void start_trace_file(struct vs_data *vsd)
{
	struct pinfo	my_id;

	char	filename[FILENAME_LENGTH];
	FILE  * f;

	pinfo_fill(&my_id);

	mkfname(filename,l_meta_filename,vsd->vecnum);
	f=fopen(filename,"w");
	write_pinfo(f,&my_id,"");
	fprintf(f,"ID_END\n");
	vsd->trace_offset=ftell(f);
	fclose(f);
}

/* sync_a_files -- do the truncation */
static int	sync_a_files(struct vs_data *vsd)
{
	struct stat	sbuf;

	char	filename[FILENAME_LENGTH];
	int	i;
	long	newsize;

	newsize=vsd->rev_a_s * STAMPSIZE;

	for(i=0;i<nbxs;i++) {
		mkfname(filename,a_meta_filename,i,vsd->vecnum,vsd->rev_a_s);
		if (stat(filename,&sbuf)<0) {
		       if (vsd->rev_a_s != 0 || errno != ENOENT) {
				die("stat(%s) : %s\n",
						1,filename,strerror(errno));
		       }
		       continue;
		}

		if (sbuf.st_size < newsize) {
			die("%s : size should be at least %d\n",
					1,filename,newsize);
		}

		if (sbuf.st_size > newsize) {
			if (truncate(filename,newsize)<0) {
				die("truncate(%s) : %s\n",
					1,filename,strerror(errno));
			}
		}

		stat(filename,&sbuf);

		if (sbuf.st_size!=newsize) {
			die("%s : unable to obtain size %d\n",
					1,filename,newsize);
		}
	}

	return 0;
}

static int	rename_a_files(struct vs_data *vsd)
{
	int i;
	char filename[2][FILENAME_LENGTH];
	
	for(i=tsv()->x0;i<tsv()->x1;i++) {
		mkfname(filename[0],a_meta_filename, i, vsd->vecnum,
				vsd->rev_a_s);
		mkfname(filename[1],a_meta_filename, i, vsd->vecnum,
				vsd->rev_a_s+1);
		if (rename(filename[0],filename[1])<0) {
			fprintf(stderr,"%s->%s: %s\n",
					filename[0],filename[1],
					strerror(errno));
		}
	}
	return 0;
}

static int	open_a_files_appendmode(struct vs_data *vsd)
{
	int i;
	char filename[FILENAME_LENGTH];
	
	for(i=0;i<nbxs;i++) {
		mkfname(filename,a_meta_filename, i, vsd->vecnum,
				vsd->rev_a_s);
		vsd->a_save[i]=fopen(filename,"a");
		if (vsd->a_save[i]==NULL) {
			die("fopen(%s) : %s\n",
					1,filename,strerror(errno));
		}
		if (ftell(vsd->a_save[i])!=vsd->rev_a_s*STAMPSIZE) {
			die("%s is wrongly positioned\n",1,filename);
		}
	}
	return 0;
}

static int	flush_a_files(struct vs_data *vsd)
{
	int i;
	
	for(i=tsv()->x0;i<tsv()->x1;i++) {
		fflush(vsd->a_save[i]);
	}
	return 0;
}

static int	duplicate_v_file(struct vs_data *vsd, int n)
{
	char filename[2][FILENAME_LENGTH];

	if (!is_master_thread())
		return 0;

	mkfname(filename[0],v_meta_filename,vsd->vecnum,vsd->rev_v+n);
	mkfname(filename[1],"%s.saved",filename[0]);

	if (link(filename[0],filename[1])<0) {
		fprintf(stderr,"%s->%s: %s\n",
				filename[0],filename[1],
				strerror(errno));
	}
	t_filename_head_insert(&(vsd->rc->vectors),
			new_t_filename(vsd->rev_v+n,filename[1]));

	return 0;
}

static int	duplicate_s_file(struct vs_data *vsd, int n)
{
	char filename[2][FILENAME_LENGTH];

	if (!is_master_thread())
		return 0;

	if (vsd->msi->differential && n==0) {
		t_filename_head_insert(&(vsd->rc->sums),
				new_t_filename(vsd->rev_a_s+n,"/dev/null"));
		return 0;
	}


	if (vsd->msi->differential) {
		mkfname(filename[0],ds_meta_filename,
				vsd->vecnum,vsd->ds_start,vsd->rev_a_s+n);
	} else {
		mkfname(filename[0],s_meta_filename,vsd->vecnum,vsd->rev_a_s+n);
	}
	mkfname(filename[1],"%s.saved",filename[0]);

	if (link(filename[0],filename[1])<0) {
		fprintf(stderr,"%s->%s: %s\n",
				filename[0],filename[1],
				strerror(errno));
	}
	t_filename_head_insert(&(vsd->rc->sums),
			new_t_filename(vsd->rev_a_s+n,filename[1]));

	if (vsd->msi->differential) {
		bw_lvblock_set_zero_full(vsd->sum);
		vsd->ds_start=vsd->rev_a_s+n;
		return 0;
	}

	return 0;
}

static int	dispose_v_file(struct vs_data *vsd)
{
	char filename[FILENAME_LENGTH];

	if (!is_master_thread())
		return 0;

	mkfname(filename,v_meta_filename,vsd->vecnum,vsd->rev_v);
	unlink(filename);

	return 0;
}

static int	dispose_s_file(struct vs_data *vsd)
{
	char filename[FILENAME_LENGTH];

	if (!is_master_thread())
		return 0;

	if (!vsd->msi->differential) {
		mkfname(filename,s_meta_filename,vsd->vecnum,vsd->rev_a_s);
		unlink(filename);
	} else if (vsd->ds_start < vsd->rev_a_s) {
		mkfname(filename,ds_meta_filename,
				vsd->vecnum,vsd->ds_start,vsd->rev_a_s);
		if (!exist(filename)) {
			printf("%s does not exist...\n",filename);
			return -1;
		}
		unlink(filename);
	}

	return 0;
}

static int	read_v_file(struct vs_data *vsd)
{
	char	filename[FILENAME_LENGTH];
	FILE  * f;

	if (!is_master_thread())
		return 0;

	/* v has been shifted by the appropriate amount on input */
	
	if (vsd->rev_v>0) {
		mkfname(filename,v_meta_filename,vsd->vecnum,vsd->rev_v);
	} else {
		mkfname(filename,y_meta_filename,vsd->vecnum);
	}

	f=fopen(filename,"r");
		
	if (f==NULL) {
		die("fopen(%s) : %s (was OK before !)\n",
				1,filename,strerror(errno));
	}

	bw_lvblock_read(vsd->vec,f);
	fclose(f);

	return 0;
}

static int	read_s_file(struct vs_data *vsd)
{
	char	filename[FILENAME_LENGTH];
	FILE  * f;

	if (!is_master_thread())
		return 0;

	/* v has been shifted by the appropriate amount on input */

	vsd->ds_start=vsd->rev_a_s;
	
	if (vsd->rev_a_s==0 || vsd->msi->differential) {
		bw_lvblock_set_zero_full(vsd->sum);
		return 0;
	}

	mkfname(filename,s_meta_filename,vsd->vecnum,vsd->rev_a_s);

	f=fopen(filename,"r");
		
	if (f==NULL) {
		die("fopen(%s) : %s (was OK before !)\n",
				1,filename,strerror(errno));
	}

	bw_lvblock_read(vsd->sum,f);
	fclose(f);

	return 0;
}

static int	save_v_file(struct vs_data *vsd)
{
	char	filename[FILENAME_LENGTH];
	FILE  * f;

	if (!is_master_thread())
		return 0;

	/* v has been shifted by the appropriate amount on input */
	
	mkfname(filename,v_meta_filename,vsd->vecnum,vsd->rev_v + 1);

	f=fopen(filename,"w");
		
	if (f==NULL)
		die("fopen(%s) : %s\n", 1,filename,strerror(errno));

	bw_lvblock_write(f,vsd->vec);
	fclose(f);

	return 0;
}

static int	save_s_file(struct vs_data *vsd)
{
	char	filename[FILENAME_LENGTH];
	FILE  * f;

	if (!is_master_thread())
		return 0;

	/* s has been shifted by the appropriate amount on input */
	
	if (vsd->msi->differential) {
		mkfname(filename,ds_meta_filename,
				vsd->vecnum,vsd->ds_start,vsd->rev_a_s+1);
	} else {
		mkfname(filename,s_meta_filename,vsd->vecnum,vsd->rev_a_s+1);
	}

	f=fopen(filename,"w");
		
	if (f==NULL)
		die("fopen(%s) : %s\n", 1,filename,strerror(errno));

	bw_lvblock_write(f,vsd->sum);
	fclose(f);

	return 0;
}

static int	rebuild_certificate(struct vs_data * vsd)
{
	char	filename[2][FILENAME_LENGTH];
	int	ctype=-1;


	if (vsd->cstart == -1) {
		return 0;
	}

	if (vsd->cstart == vsd->rev_v) {
		if (vsd->tenant) free(vsd->tenant);
		vsd->tenant=NULL;
		return 0;
	}

	switch(task_type) {
		case TASK_DOTPRODUCTS:
			ctype=CERTIF_DOTPRODUCT;
			break;
		case TASK_POLYNOMIAL:
			ctype=CERTIF_POLYNOMIAL;
			break;
	}
	vsd->rc=new_certificate(ctype,vsd->vecnum,vsd->cstart);

	/* Transfer the producer info to the certificate structure */
	vsd->rc->producer=vsd->tenant;
	vsd->tenant=NULL;

	/* Put the originating v file */
	if (vsd->cstart>0) {
		mkfname(filename[0],v_meta_filename,vsd->vecnum,vsd->cstart);
		mkfname(filename[1],"%s.saved",filename[0]);
	} else {
		mkfname(filename[1],y_meta_filename,vsd->vecnum);
	}
	if (!exist(filename[1])) {
		die("Cannot find %s\n",1,filename[1]);
	}
	t_filename_head_insert(&(vsd->rc->vectors),
			new_t_filename(vsd->cstart,filename[1]));


	return 0;
}

static int	issue_sync_certificate(struct vs_data * vsd, int n)
{
	struct pinfo my_id;
	int	i;
	char	filename[FILENAME_LENGTH];
	FILE  * src,* dst;

	if (vsd->rc==NULL) {
		/* Might happen at the beginning */
		return 0;
	}

	switch(task_type) {
		case TASK_DOTPRODUCTS:
			for(i=0;i<nbxs;i++) {
				mkfname(filename,a_meta_filename,first.i+i,
						vsd->vecnum, vsd->rev_v + n);
				src=fopen(filename,"r");
				if (src==NULL)
					perror(filename);
				mkfname(filename,a_sub_meta_filename,first.i+i,
						vsd->vecnum, vsd->cstart,
						vsd->rev_v + n);
				dst=fopen(filename,"w");
				fseek(src,vsd->cstart*STAMPSIZE,SEEK_SET);
				fcopy(dst,src,
					(vsd->rev_v+n-vsd->cstart)*STAMPSIZE);
				fclose(dst);
				fclose(src);
				t_filename_tail_insert(&(vsd->rc->dotproducts),
						new_t_filename(i,filename));
			}
			break;
		case TASK_POLYNOMIAL:
			sprintf(filename,f_meta_filename,
					vsd->vecnum,vsd->msi->solution_col,
					vsd->msi->t_value);
			src = fopen(filename,"r");
			sprintf(filename,f_sub_meta_filename,
					vsd->vecnum,vsd->msi->solution_col,
					vsd->msi->t_value,vsd->cstart,
					vsd->rev_v + n);
			dst = fopen(filename,"w");
			fseek(src,vsd->msi->valuation*
					bw_filesize*sizeof(mp_limb_t)+
					vsd->cstart*STAMPSIZE,
					SEEK_SET);
			fcopy(dst,src,(vsd->rev_v+n-vsd->cstart)*STAMPSIZE);
			fclose(dst);
			fclose(src);
			/* Put the last s file */
			duplicate_s_file(vsd, n);
			break;
	}

	/* Put the last v file */
	duplicate_v_file(vsd, n);

	vsd->rc->r1=vsd->rev_v + n;

	if (vsd->rc->producer) {
		write_certificate(vsd->rc);
		pinfo_fill(&my_id);
		write_pinfo(vsd->rc->f,&my_id,"# <rebuilder info> ");
	} else {
		write_certificate(vsd->rc);
	}

	close_certificate(vsd->rc);
	vsd->rc=NULL;

	return 0;
}

static int	start_next_certificate(struct vs_data * vsd, int k)
{
	char	filename[2][FILENAME_LENGTH];
	int	ctype=-1;

	vsd->cstart=vsd->rev_v + k;

	switch(task_type) {
		case TASK_DOTPRODUCTS:
			ctype=CERTIF_DOTPRODUCT;
			break;
		case TASK_POLYNOMIAL:
			ctype=CERTIF_POLYNOMIAL;
			break;
	}
	vsd->rc=new_certificate(ctype,vsd->vecnum,vsd->cstart);

	/*
	 * There's no need to duplicate the file, since it has already
	 * been duplicated by the certificate save procedure. However, we
	 * do need to insert the file in the certificate.
	 */

	if (vsd->cstart>0) {
		mkfname(filename[0],v_meta_filename,vsd->vecnum,vsd->cstart);
		mkfname(filename[1],"%s.saved",filename[0]);
	} else {
		mkfname(filename[1],y_meta_filename,vsd->vecnum);
	}

	/* When we start a new job, it is really likely that we get the
	 * wrong number of links. Correct this...
	 */
	if (!exist(filename[1])) {
		fprintf(stderr,"WARNING: %s does not exist (duplicated)\n",
				filename[1]);
		/* We choose to duplicate. */
		link(filename[0],filename[1]);
	}

	t_filename_head_insert(&(vsd->rc->vectors),
			new_t_filename(vsd->cstart,filename[1]));


	if (task_type==TASK_POLYNOMIAL) {
	       	if (vsd->cstart>0 && !vsd->msi->differential) {
			mkfname(filename[0],s_meta_filename,
					vsd->vecnum,vsd->cstart);
			mkfname(filename[1],"%s.saved",filename[0]);
			t_filename_head_insert(&(vsd->rc->sums),
				new_t_filename(vsd->cstart,filename[1]));
		} else {
			t_filename_head_insert(&(vsd->rc->sums),
				new_t_filename(vsd->cstart,"/dev/null"));
		}
	}

	return 0;
}

/***************************** ENTRY POINTS **********************************/

int	vs_try_startup_sync(struct vs_data ** pvsd, int vecnum, int diff)
{
	DIR   * dirp;

#define vsd	(*pvsd)

	if (!is_master_thread()) {
		die("Must only be called from the master thread\n",1);
	}

	assert(nbxs==m_param);
	/* Hélas, j'ai eu la flemme d'implémenter un truc qui marche même
	 * dans un autre cas */

	dirp=opendir(wdir_filename);

	printf("Vector %d ; beginning survey\n",vecnum);

	vsd=malloc(sizeof(struct vs_data));

	vsd->vecnum=vecnum;
	vsd->rev_v=lookup_v_files(dirp,vecnum);
	switch(task_type) {
		case TASK_DOTPRODUCTS:
			vsd->rev_a_s=lookup_a_files(dirp,vecnum);
			break;
		case TASK_POLYNOMIAL:
			if (diff) {
				printf("Enabling differential mode\n");
				vsd->rev_a_s = vsd->rev_v;
			} else {
				vsd->rev_a_s=lookup_s_files(dirp,vecnum);
			}
			break;
	}

	closedir(dirp);

	if (vsd->rev_a_s!=vsd->rev_v) {
		fprintf(stderr,"Inconsistent revisions\n");
		exit(1);
	}

	vsd->rc=NULL;
	vsd->tenant=NULL;
	vsd->vec=NULL;
	vsd->sum=NULL;
	vsd->a_save=NULL;
	vsd->msi=NULL;

	if (task_type==TASK_DOTPRODUCTS) {
		vsd->a_save=malloc(nbxs*sizeof(FILE*));
		memset(vsd->a_save,0,nbxs*sizeof(FILE*));
	}

	vsd->cstart=lookup_trace_file(vsd);

	return vsd->rev_v;
#undef vsd
}

int	vs_do_startup_sync(struct vs_data * vsd)
{
	if (!is_master_thread()) {
		die("Must only be called from the master thread\n",1);
	}

	rebuild_certificate(vsd);

	if (task_type==TASK_DOTPRODUCTS)
		sync_a_files(vsd);

	issue_sync_certificate(vsd,0);

	if (task_type==TASK_DOTPRODUCTS)
		open_a_files_appendmode(vsd);
	else
		read_s_file(vsd);

	read_v_file(vsd);

	start_trace_file(vsd);

	start_next_certificate(vsd,0);

	write_trace_file(vsd,0);

	return 0;
}

int	vs_checkpoint(struct vs_data * vsd, int force_certif)
{
	barrier_wait(&vs_checkpoint_barrier,NULL,NULL);

	switch(task_type) {
		case TASK_DOTPRODUCTS:
			/* Done concurrently */
			flush_a_files(vsd);
			rename_a_files(vsd);
			break;
		case TASK_POLYNOMIAL:
			/* done only by the master thread */
			save_s_file(vsd);	/* save the new one */
			dispose_s_file(vsd);	/* discard the old one */
			break;
		default:
			BUG();
			break;
	}

	/* done only by the master thread */
	save_v_file(vsd);	/* save the new one */
	dispose_v_file(vsd);	/* discard the old one */

	/* This is put here because rename_a_files bus absolutely be
	 * finished by the time we enter issue_sync_certificate */
	barrier_wait(&vs_checkpoint_barrier,NULL,NULL);

	if (is_master_thread() && (
		force_certif || (vsd->rev_v+1) % keep_periodicity == 0))
	{
		issue_sync_certificate(vsd,1);
		start_next_certificate(vsd,1);
	}

	write_trace_file(vsd,1);

	barrier_wait(&vs_checkpoint_barrier,NULL,NULL);

	if (is_master_thread()) {
		/* Dangerous! touches global data */
		vsd->rev_a_s++;
		vsd->rev_v++;
	}

	return 0;
}
