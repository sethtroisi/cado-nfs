#ifndef RECOVERY_H_
#define RECOVERY_H_

#ifdef	__cplusplus
extern "C" {
#endif

#define	STAMPSIZE	(bw_filesize*sizeof(mp_limb_t)*periodicity)

/* Ça, ça va certainement finir dans le .h */

struct vs_data {
	struct certificate    * rc;
	int			rev_a_s;
	int			rev_v;
	int			vecnum;
	FILE		     ** a_save;
	int			cstart;	/* redundant wirh rc->r0 when
					   the latter exists */
	int			ds_start;	/* I guess it's also redundant*/
	long			trace_offset;
	struct pinfo	      * tenant;
	bw_vector_block		vec;
	bw_vector_block		sum;
	struct mksol_info_block * msi;
};

	
/* Points d'entrée. Toutes les fonctions peuvent être appelées en
 * contexte MT, quitte à rendre la main immédiatement pour n-1 threads.
 */

int	vs_try_startup_sync(struct vs_data **, int, int);
/* Fait l'inventaire des fichiers pertinents qui sont présents.
 * Jette à la moindre incohérence ; leur résolution est à la charge de
 * l'utilisateur.
 * Dit quelle révision peut être atteinte, au vu des fichiers présents.
 * Remplit la structure struct vs_data avec ces infos.
 *
 * Révision 0 = on part from scratch.
 *
 * Fait aussi l'inventaire des certificats présents, et cherche un
 * éventuel fichier témoin.
 *
 * RETURN: La révision.
 */

int	vs_do_startup_sync(struct vs_data *);
/* Fait la mise à jour comme demandé.
 * Reconstruit un certificat jusqu'au stamp dernièrement atteint.
 *
 * RETURN: 0, -1 en cas d'erreur, +errno.
 */

int	vs_checkpoint(struct vs_data *, int);
/* Passe d'un stamp à l'autre.
 * Renomme les fichiers.
 * Change de certificat courant s'il est temps d'en mettre un autre.
 */

/* Les clients ont aussi le droit d'écrire dans a_file */

#ifdef	__cplusplus
}
#endif

#endif	/* RECOVERY_H_ */
