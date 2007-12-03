
/* WARNING: This file miscompiles with digital Unix's cc, at least on
 * version 4.0F (read: the default compiler that ships with 4.0F - I
 * don't know this compiler's exact version number), when optimization
 * level -O4 is used.
 *
 * To get something that runs, switch to -O3 at most.
 */

/* #define _BSD_SOURCE */
#include <stdarg.h>
#include <ctype.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#ifdef __SVR4
/* Buggy solaris 2.5 crap */
extern char * strdup(const char *);
extern char * index(const char *, int);
#endif
#include "options.h"
#include <assert.h>

#ifndef MAX
#define MAX(a,b)	(((a)>(b))?(a):(b))
#endif
	
#ifndef MIN
#define MIN(a,b)	(((a)<(b))?(a):(b))
#endif
	
#ifdef __DECC
#if __DECC_VER <= 50990005
#warning "please turn optimization off for this file"
#endif
#pragma nostandard
#endif

#define char_in_string(s,c) ((c) && index((s),(c))!= NULL)

const struct textflag_desc yesno_flag[]={
		{"yes y on",1},
		{"no n off",0},
		{NULL,0}};

struct opt_parser_context;
static void opt_die(const char *, int, struct opt_parser_context *, ...);

#define NO_CTX ((struct opt_parser_context *)NULL)

void new_option(struct opt_desc ** p_optlist, int * p_nopts, ...)
{
	va_list ap;
	struct param_desc * res;
	int i;
	int npar;

	if (*p_optlist==NULL) {
		*p_optlist=malloc(sizeof(struct opt_desc));
	} else {
		*p_optlist=realloc(*p_optlist,
				(1+*p_nopts)*sizeof(struct opt_desc));
	}
	va_start(ap,p_nopts);
	(*p_optlist)[*p_nopts].names	= va_arg(ap,char *);
	npar=va_arg(ap,int);
	res=malloc((npar+1)*sizeof(struct param_desc));
	for(i=0;i<npar;i++) {
		res[i].type=va_arg(ap,unsigned int);
		res[i].default_value=va_arg(ap,char*);
		res[i].additional_data=va_arg(ap,void*);
	}
	res[i].type=0;
	res[i].default_value=NULL;
	res[i].additional_data=NULL;
	(*p_optlist)[*p_nopts].params	= res;
	(*p_optlist)[*p_nopts].process	= va_arg(ap,opt_handler*);
	(*p_optlist)[*p_nopts].argument	= va_arg(ap,void*);
	(*p_optlist)[*p_nopts].flags	= va_arg(ap,int);
	(*p_optlist)[*p_nopts].doc	= va_arg(ap,char*);
	(*p_nopts)++;
	va_end(ap);
}

static char* canonical_name(struct extended_option_desc * opt)
{
	static char buf[12];
	if (opt->long_names) {
		return opt->long_names[0];
	} else if (opt->short_names) {
		sprintf(buf,"-%c",opt->short_names[0]);
	} else {
		sprintf(buf,"(null, #%2d)",opt->number);
	}
	return buf;
}

int 
func_builtin_raise_flag(char **argv, struct extended_option_desc * conf, void * p_val)
{
	*(int *) p_val = 1;
	return 1;
}

int 
func_builtin_process_int(char **argv, struct extended_option_desc * conf, void * p_val)
{
	char           *endptr;

	*(int *) p_val = (int) strtol(argv[0], &endptr, 10);

	return strlen(endptr) == 0;
}

int 
func_builtin_increase_int(char **argv, struct extended_option_desc * conf, void * p_val)
{
	*(int *) p_val = conf->n_hits;

	return 1;
}

int 
func_builtin_process_long(char **argv, struct extended_option_desc * conf, void * p_val)
{
	char           *endptr;

	*(long *) p_val = strtol(argv[0], &endptr, 10);

	return strlen(endptr) == 0;
}

int 
func_builtin_process_mem(char **argv, struct extended_option_desc * conf, void * p_val)
{
	char           *endptr;

	*(int *) p_val = (int) strtol(argv[0], &endptr, 0);

	switch(*endptr) {
		case 'K' :
		case 'k' :	*(int *) p_val *= (1<<10);
				endptr++;
				break;
		case 'M' :
		case 'm' :	*(int *) p_val *= (1<<20);
				endptr++;
		default  :
				break;
	}
	if (char_in_string("BbOo",*endptr))
		endptr++;

	return strlen(endptr) == 0;
}

static int 
word_in_list(const char *haystack, const char *needle)
{
	char           *pos;

	pos = strstr(haystack, needle);
	for(;pos && pos[strlen(needle)]!=' ' && pos[strlen(needle)]!=0;) {
		pos += strlen(needle);
		pos = strstr(pos, needle);
	}

	if (pos == NULL)
		return 0;

	if (pos != haystack && *(pos - 1) != ' ')
		return 0;

	pos += strlen(needle);

	if (*pos != ' ' && *pos != (char) 0)
		return 0;

	return 1;
}

int 
func_builtin_process_textflag(
		char **argv,
		struct extended_option_desc * conf,
		void * p_val)
{
	int i;
	struct textflag_desc * t;

	for(i=0;i<conf->n_params;i++) {
		for(t=conf->params[i].additional_data;t->all_names;t++) {
			if (word_in_list(t->all_names,argv[i])) {
				((int*)p_val)[i]=t->flag_value;
				break;
			}
		}
		if (t->all_names==NULL)
			opt_die("Flag value %s unknown for argument %s\n",
					1,NO_CTX,
					argv[i],canonical_name(conf));
	}
	return 1;
}

int 
func_builtin_process_string_cpy(char **argv, struct extended_option_desc * conf, void * p_val)
{
	strcpy(p_val, argv[0]);
	return 1;
}

int 
func_builtin_process_string_alloc(char **argv, struct extended_option_desc * conf, void * p_val)
{
	if (*(char **) p_val != NULL) {
		printf("\nWarning, changing parameter of option %s."
			" Previous value freed\n",
				canonical_name(conf));
		free(*(char **) p_val);
	}
	*(char **) p_val = strdup(argv[0]);
	return 1;
}

static int deref_strcmp(const void *a, const void *b)
{
	return strcmp(*(char**)a,*(char**)b);
}

static void check_ambiguous(struct extended_option_desc * l, int n)
{
	int total=0;
	int i,j,k;
	char ** blah;

	for(i=0;i<n;i++) total+=l[i].n_long_names+l[i].n_short_names;

	if (total==0)
		return;

	blah=malloc(total*sizeof(char*));
	k=0;
	for(i=0;i<n;i++) {
		for(j=0;j<l[i].n_long_names;j++) {
			blah[k++]=strdup(l[i].long_names[j]);
		}
		for(j=0;j<l[i].n_short_names;j++) {
			blah[k]=malloc(3);
			sprintf(blah[k],"-%c",l[i].short_names[j]);
			k++;
		}
	}
	assert(k==total);
	qsort(blah,k,sizeof(char*),&deref_strcmp);
	for(i=1;i<total;i++) {
		if (strcmp(blah[i],blah[i-1])!=0)
			continue;
		opt_die("Options %s is multiply defined\n",1,NO_CTX,blah[i]);
	}

	for(i=0;i<total;i++) {
		free(blah[i]);
	}
	free(blah);
}

static struct extended_option_desc *
preparse_optlist(struct opt_desc list[], int n)
{
	struct extended_option_desc * res;
	int i;

	res = malloc(n*sizeof(struct extended_option_desc));
	if (res==NULL) return NULL;

	for(i=0;i<n;i++) {
		int			j_long,j_short;
		const char	      * s;
		int			space_approx;
		struct param_desc     * v;

		/* get a trivial upper bound on the number of names */

		s=list[i].names;

		if (s!=NULL) {
			for(space_approx=0;*s!='\0';space_approx++) {
				s+=strcspn(s," ");
				s+=strspn(s," ");
			}
		} else {
			res[i].long_names	= NULL;
			res[i].n_long_names	= 0;
			res[i].short_names	= NULL;
			res[i].n_short_names	= 0;
			res[i].matches_null	= 1;
			goto no_opt_list;
		}
		res[i].long_names	= malloc(space_approx*sizeof(char*));;
		res[i].n_long_names	= 0;
		res[i].short_names	= malloc(space_approx+1);
					/* +1 : safety measure */
		res[i].n_short_names	= 0;
		j_long=j_short=0;
		for(s=list[i].names;*s!='\0';s+=strspn(s," ")) {
			int len;
			len=strcspn(s," ");
			if (strncmp(s,"--",2)==0) { /* long option */
#define t	res[i].long_names[j_long]
				t=malloc(len+1);
				memcpy(t,s,len);
				t[len]='\0';
				if (strlen(t)==2) {
					opt_die(
	"Option %s : double dash alone forbidden as an option name\n",1,NO_CTX,
						list[i].names);
				}
				if (isdigit((int) t[2])) {
					opt_die(
	"Option %s : leading digits forbidden in long options\n", 1, NO_CTX, t);
				}
				if (strcspn(t,",=") != len) {
					opt_die(
	"Option %s : = or , forbidden in option names\n", 1, NO_CTX,  t);
				}
				j_long++;
#undef t
			} else if (s[0]=='-') { /* short option */
				if (len==1) {
					opt_die(
	"Option %s : single dash alone forbidden as an option name\n", 1, NO_CTX, 
						list[i].names);
				}
				if (len!=2) {
					opt_die(
	"Option %s : single dash options must be one-letter long\n", 1, NO_CTX, 
						list[i].names);
				}
				if (char_in_string("0123456789=,",s[1])) {
					opt_die(
	"Option -%c : character %c is forbidden\n", 1, NO_CTX,  s[1], s[1]);
				}
				res[i].short_names[j_short]=s[1];
				j_short++;
			} else {
				opt_die(
	"Option %s : must have a leading - or --\n",1, NO_CTX, list[i].names);
				exit(1);
			}
			s+=len;
		}
		res[i].n_long_names	= j_long;
		res[i].n_short_names	= j_short;
		res[i].long_names	= j_long?realloc(res[i].long_names,
						j_long*sizeof(char*)):NULL;
		res[i].short_names	= realloc(res[i].short_names,j_short+1);
		res[i].short_names[j_short]	= '\0';
					/* useless fool-proof measure */
		res[i].matches_null	= (list[i].flags & OPT_MATCH_NULL)!=0;
no_opt_list:
		res[i].process	= list[i].process;
		res[i].argument	= list[i].argument;
		if (list[i].flags & (OPT_LASTONLY | OPT_NOTLAST) &&
				list[i].names) {
			opt_die("Option %s : OPT_LASTONLY and OPT_NOTLAST only "
				"allowed for barewords\n",1,NO_CTX,
				canonical_name(res+i));
		}
		res[i].flags	= list[i].flags;
		res[i].doc	= list[i].doc;
		res[i].params	= list[i].params;
		res[i].n_hits	= 0;
		res[i].n_params	= 0;
		res[i].minimum_params	= 0;
		res[i].number	= i;
		for(v=res[i].params;v->type!=0;v++) {
			res[i].n_params++;
			if (v->default_value==NULL) {
				res[i].minimum_params++;
				if (res[i].minimum_params < res[i].n_params) {
					opt_die("Option %s : cannot have a default"
					" value for argument %d and not all"
					" the previous ones\n",1, NO_CTX, 
					canonical_name(res+i),
					res[i].n_params-1);
				}
			}
		}
		res[i].found_params = NULL;
	}
	check_ambiguous(res,n);
	return res;
}

struct opt_parser_context {
	int				argc;
	char			     ** argv;
	struct extended_option_desc   * opt;
	int				n_opts;
	int				curr_arg;
	int				pos_in_arg;
};

static void error_nomatch(struct opt_parser_context *ctx, char *s)
{
	opt_die("no match found for argument %s\n", 1, ctx, s);
}

static void error_nomatch_short(struct opt_parser_context *ctx, char c)
{
	opt_die("no match found for argument -%c\n", 1, ctx, c);
}

static struct extended_option_desc *
try_long_matcher(struct opt_parser_context *ctx, char *s, int l)
{
#define i ctx->curr_arg
	int j,k;
	for(j=0;j<ctx->n_opts;j++) {
		if (ctx->opt[j].n_hits && (ctx->opt[j].flags & OPT_ONCE)) {
			if (ctx->opt[j].flags & OPT_IGNFURTHER)
				continue;
		}
		for(k=0;k<ctx->opt[j].n_long_names;k++) {
			if (strncmp(ctx->opt[j].long_names[k],s,l)==0) {
				if (ctx->opt[j].n_hits &&
					(ctx->opt[j].flags & OPT_ONCE)) {
						opt_die(
	"Option %s (#%d) matches a 2nd time against %s !\n",1, ctx, 
					canonical_name(&(ctx->opt[j])),
					j, ctx->argv[i]);
				}
				ctx->opt[j].n_hits++;
				return &(ctx->opt[j]);
			}
		}
	}
	return NULL;
}
static struct extended_option_desc *
try_short_matcher(struct opt_parser_context *ctx, char c)
{
	int j,k;
	for(j=0;j<ctx->n_opts;j++) {
		if (ctx->opt[j].n_hits && (ctx->opt[j].flags & OPT_ONCE)) {
			if (ctx->opt[j].flags & OPT_IGNFURTHER)
				continue;
		}
		for(k=0;k<ctx->opt[j].n_short_names;k++) {
			if (ctx->opt[j].short_names[k]==c) {
				if (ctx->opt[j].n_hits &&
					(ctx->opt[j].flags & OPT_ONCE)) {
						opt_die(
	"Option %s (#%d) matches a 2nd time against %s !\n",1, ctx, 
					canonical_name(&(ctx->opt[j])),
					j, ctx->argv[i]);
				}
				ctx->opt[j].n_hits++;
				return &(ctx->opt[j]);
			}
		}
	}
	return NULL;
#undef i
}

static void
do_call_with_args(struct extended_option_desc * od)
{
	int i;

	if (od->n_params) {
		if (od->long_names) {
			printf("%s=",od->long_names[0]);
		} else if (od->short_names) {
			printf("-%c=",od->short_names[0]);
		}
		for(i=0;i<od->n_params;i++)
			printf("%s%c",od->found_params[i],
					(i==(od->n_params-1))?' ':',');
	} else {
		if (od->long_names) {
			printf("%s ",od->long_names[0]);
		} else if (od->short_names) {
			printf("-%c ",od->short_names[0]);
		} else {
			assert(0);
		}
	}
	(*(od->process))(od->found_params,od,od->argument);
	for(i=0;i<od->n_params;i++) free(od->found_params[i]);
	free(od->found_params);
	od->found_params=NULL;
}

static void
do_call_with_default_args(struct extended_option_desc * od)
{
	if (od->n_params) {
		int i;
		od->found_params = malloc(od->n_params * sizeof(char*));
		for(i=0;i<od->n_params;i++)
			od->found_params[i]=strdup(od->params[i].default_value);
	}

	do_call_with_args(od);
}
	
static int
build_parms_from_string(struct extended_option_desc * od,
	       struct opt_parser_context *ctx,
       	       char *t, char *s)
{
	int i;
	int n;
	int k;

	k=0;

	od->found_params = malloc(od->n_params * sizeof(char*));
	for(i=0;i<od->n_params && s[k] != '\0';i++) {
		if (od->params[i].type == OPTPARAM_NUMERAL) {
			n=strspn(s+k,"0123456789");
			if (n==0) {	/* for integers, NULL=default */
				if (od->params[i].default_value==NULL) {
					opt_die("Argument %s : no default for "
				 		"parameter %d of option %s\n",1,
						ctx, t, i,canonical_name(od));
				}
				od->found_params[i]= strdup(
						od->params[i].default_value);
				continue;
			}
			od->found_params[i]=malloc(n+1);
			strncpy(od->found_params[i],s+k,n);
			od->found_params[i][n]='\0';
			k+=n;
			if (s[k]==',') k++;	/* not necessary */
		} else {
			n=strcspn(s+k,",");
			od->found_params[i]=malloc(n+1);
			strncpy(od->found_params[i],s+k,n);
			od->found_params[i][n]='\0';
			k+=n;
			if (s[k]==',') k++; /* otherwise, s[k]==0 */
		}
	}
	for(;i<od->n_params;i++) {
		if (od->params[i].default_value) {
			od->found_params[i]=strdup(od->params[i].default_value);
		} else {
			opt_die("Argument %s : no default for "
				"parameter %d of option %s\n",1,
				ctx, t, i,canonical_name(od));
		}
	}
	return k;
}

static int
build_parms_from_cmdline(struct extended_option_desc * od,
			struct opt_parser_context *ctx)
{
	int i,k;
	od->found_params = malloc(od->n_params * sizeof(char*));
	for(i=0;i<od->n_params && ((ctx->curr_arg+i) < ctx->argc) ;i++) {
		od->found_params[i]=strdup(ctx->argv[ctx->curr_arg+i]);
	}
	k=i;
	for(;i<od->n_params;i++) {
		if (od->params[i].default_value) {
			od->found_params[i]=strdup(od->params[i].default_value);
		} else {
			opt_die("Argument %s : no default for "
				"parameter %d of option %s\n",1,
				ctx, ctx->argv[ctx->curr_arg-1],
				i,canonical_name(od));
		}
	}
	return k;
}


static int
try_to_match_null(struct opt_parser_context *ctx)
{
#define i ctx->curr_arg
	int j;
	int res=0;
	for (j = 0; j < ctx->n_opts; j++) {
		if (!ctx->opt[j].matches_null)
			continue;
		if ((ctx->opt[j].flags & OPT_LASTONLY) && i != (ctx->argc-1))
			continue;
		if ((ctx->opt[j].flags & OPT_NOTLAST) && i == (ctx->argc-1))
			continue;
		if (ctx->opt[j].n_hits && (ctx->opt[j].flags & OPT_ONCE)) {
			if (ctx->opt[j].flags & OPT_IGNFURTHER)
				continue;
			opt_die(
	"Option %s (#%d) matches a 2nd time against %s !\n",1, ctx, 
					canonical_name(&(ctx->opt[j])),
					j, ctx->argv[i]);
		}
		/* HIT ! */
		ctx->opt[j].n_hits++;
		i+=build_parms_from_cmdline(ctx->opt+j,ctx);
		do_call_with_args(ctx->opt+j);
		res=1;
		break;
	}

	return res;
#undef i
}



#define i ctx->curr_arg
#define s ctx->argv[i]
#define j ctx->pos_in_arg

static void process_long_option(struct opt_parser_context *ctx)
{
	struct extended_option_desc * od;
	/*
	 * Handle the case of the skip directives
	 */
	if (s[2]=='\0') {
		/* skip all remaining args. */
		i++;
		for(;i<ctx->argc && try_to_match_null(ctx););
		if (i<ctx->argc) error_nomatch(ctx, s);
		return;
	} else if (isdigit((int) s[2])) {
		/* skip a given # of args */
		char *c;
		unsigned int n;
		n=strtoul(s+2,&c,10);
		if (*c=='\0') {
			i++;
			n=MIN(n+i,ctx->argc);
			for(;i<n && try_to_match_null(ctx););
			if (i<n) error_nomatch(ctx, s);
			return;
		} else {
			/* huh, there are characters behind... */
			opt_die("Argument %s : leading digits forbidden.\n",
					1,ctx,s);
		}
	}
	/*
	 * Look for the option employed.
	 */
	j=strcspn(s,"= ");
	od=try_long_matcher(ctx,s,j);	/* j : length */
	if (od==NULL) error_nomatch(ctx, s);
	/*
	 * has pending parameters ?
	 */
	if (char_in_string("= ",s[j])) {
		if (od->n_params==0) {
			opt_die("Argument %s : option %s takes no parameters\n", 1,
					ctx, s, canonical_name(od));
		}
		build_parms_from_string(od,ctx,s,s+j+1);
		do_call_with_args(od);
		i++;
	/*
	 * Now, necessarily : s[j]=0
	 */
	} else if (od->minimum_params==0) {
			/* also captures the case of 0 params */
			do_call_with_default_args(od);
			i++;
			return;
	} else {
		if (i == ctx->argc-1) {
			opt_die("Cannot find parameters for arg. %s "
				"(option %s needs at least %d)\n",1,
				ctx, s, canonical_name(od),
				od->minimum_params);
		}
		i++;
		if (char_in_string(s,',') &&
				od->n_params > 1 &&
				!(od->flags & OPT_NOPARSE_COMMA)) {
			build_parms_from_string(od,ctx,ctx->argv[i-1],s);
			do_call_with_args(od);
			i++;
			return;
		}
		i+=build_parms_from_cmdline(od,ctx);
		do_call_with_args(od);
	}
}

static void process_short_option(struct opt_parser_context *ctx)
{
	struct extended_option_desc * od;

	od=try_short_matcher(ctx,s[j]);
	if (od==NULL) error_nomatch_short(ctx, s[j]);
	/*
	 * has pending parameter (explicitly introduced) ?
	 */
	if (char_in_string("= ",s[j+1])) {
		if (od->n_params==0) {
			opt_die("Argument %s : option %s takes no parameters\n", 1,
					ctx, s, canonical_name(od));
		}
		j+=2+build_parms_from_string(od,ctx,s,s+j+2);
		do_call_with_args(od);
	} else if (isdigit((int) s[j+1]) && od->n_params &&
				 od->params[0].type == OPTPARAM_NUMERAL) {
		j+=1+build_parms_from_string(od,ctx,s,s+j+1);
		do_call_with_args(od);
	} else if (od->minimum_params==0) {
		do_call_with_default_args(od);
		j++;
	} else {
		if (i == ctx->argc-1) {
			opt_die("Cannot find parameters for arg. %s "
				"(option %s needs at least %d)\n",1,
				ctx, s, canonical_name(od),
				od->minimum_params);
		}
		i++;
		if (char_in_string(s,',') &&
				od->n_params > 1 &&
			       	!(od->flags & OPT_NOPARSE_COMMA)) {
			build_parms_from_string(od,ctx,ctx->argv[i-1],s);
			do_call_with_args(od);
			i++;
			return;
		}
		i+=build_parms_from_cmdline(od,ctx);
		do_call_with_args(od);
		return;
	}
	if (s[j]==',')  j++;
	if (s[j]=='\0') i++;
	return;
}


static void process_one_option(struct opt_parser_context *ctx)
{
	if (s[0]!='-' || (s[1]=='\0')) {
		/* single dash is treated as a bareword */
		if (try_to_match_null(ctx)) {
			return;
		} else error_nomatch(ctx, s);
	}

	if (s[1]=='-') {	/* long option */
		process_long_option(ctx);
	} else {
		int s_i;
		j=1;
		s_i=i;
		for(;s_i==i;) process_short_option(ctx);
	}
}

static void in_process_options(struct opt_parser_context *ctx)
{
	for (i = 1; i < ctx->argc;)
		process_one_option(ctx);

}

#undef i
#undef j
#undef s

static void showuse(struct opt_parser_context * c)
{
	/* TODO TBC */
}


static void opt_die(const char *fmt, int signal, struct opt_parser_context * c, ...)
{
	va_list         ap;

#ifdef __SUNPRO_C
	va_start(ap, (void) c);
#else
	va_start(ap, c);
#endif
	vfprintf(stderr, fmt, ap);
	va_end(ap);
	if (c!=NULL)
		showuse(c);
	exit(signal);
}


void process_options(int argc, char *argv[], int n_opts, struct opt_desc * tab)
{
	struct opt_parser_context m_ctx;
	struct opt_parser_context * ctx = & m_ctx;
	int i;
	int j;
	
	m_ctx.argc	= argc;
	m_ctx.argv	= argv;
	m_ctx.opt	= preparse_optlist(tab,n_opts);
	m_ctx.n_opts	= n_opts;
	m_ctx.curr_arg	= 1;
	m_ctx.pos_in_arg	= 0;

	printf("Parsing options: ");
	in_process_options(&m_ctx);
	printf("\n");

	for(i=0;i<n_opts;i++) {
		if (ctx->opt[i].flags & OPT_REQUIRED && !ctx->opt[i].n_hits) {
			opt_die("Option %s must be present\n",1,
					ctx,canonical_name(ctx->opt+i));
		}
	}

	for(i=0;i<n_opts;i++) {
		for(j=0;j<ctx->opt[i].n_long_names;j++)
			free(ctx->opt[i].long_names[j]);
		if (ctx->opt[i].long_names)	free(ctx->opt[i].long_names);
		if (ctx->opt[i].short_names)	free(ctx->opt[i].short_names);
		if (tab[i].params)		free(tab[i].params);
	}
	free(ctx->opt);
}
