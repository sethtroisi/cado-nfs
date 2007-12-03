#ifndef OPTIONS_PROCESSING_H_
#define OPTIONS_PROCESSING_H_

#ifdef	__cplusplus
extern "C" {
#endif

struct textflag_desc;
struct param_desc;
struct opt_desc;
struct extended_option_desc;

typedef int (opt_handler) (char **, struct extended_option_desc *, void *);

struct opt_desc {
	const char	      *	names;
	struct param_desc     * params;
	opt_handler	      * process;
	void		      * argument;
	int			flags;
	char		      * doc;
};

struct extended_option_desc {
	struct param_desc     *	params;
	opt_handler	      * process;
	void		      *	argument;
	unsigned int		flags;
	char		      * doc;
	/* These must not be used. Use canonical_name() instead */
	char		     **	long_names;
	int			n_long_names;
	char		      *	short_names;
	int			n_short_names;
	/* The rest is undocumented */
	int			n_params;
	int			minimum_params;
	char		     ** found_params;
	int			n_hits;
	int			matches_null;
	int			number;
};

struct textflag_desc {
	const char	      * all_names;
	unsigned int		flag_value;
};

struct param_desc {
	unsigned int		type;
	char		      * default_value;	/* NULL if no default */
	void		      * additional_data;
};



#define OPTPARAM_NUMERAL	1
#define OPTPARAM_TEXTUAL	2

#define PARAM_NUMERAL(x)	OPTPARAM_NUMERAL,x,NULL
#define PARAM_TEXTUAL(x)	OPTPARAM_TEXTUAL,x,NULL
#define PARAM_TEXTFLAG(x,y)	OPTPARAM_TEXTUAL,x,y


typedef char * arglist[];

void process_options(int, char *[], int, struct opt_desc *);
/* extern struct param_desc      *	compile_param_desc(int, ...); */
/* #define PARAMLIST	compile_param_desc */
extern void new_option(struct opt_desc **, int *, ...);


#define OPT_ONCE		0x0001
#define OPT_REQUIRED		0x0002
#define OPT_IGNFURTHER		0x0004
#define OPT_INDEP		(OPT_ONCE|OPT_IGNFURTHER)
#define OPT_LASTONLY		0x0008
#define OPT_NOTLAST		0x0010
#define OPT_MATCH_NULL		0x0020
#define OPT_NOPARSE_COMMA	0x0040
#define OPT_UNDOCUMENTED	0x0080

/* Built-in flags */
extern const struct textflag_desc yesno_flag[];

/* Built-in handlers */
extern opt_handler func_builtin_raise_flag;
extern opt_handler func_builtin_process_int;
extern opt_handler func_builtin_increase_int;
extern opt_handler func_builtin_process_long;
extern opt_handler func_builtin_process_mem;
extern opt_handler func_builtin_process_textflag;
extern opt_handler func_builtin_process_string_cpy;
extern opt_handler func_builtin_process_string_alloc;

#define builtin_raise_flag		&func_builtin_raise_flag
#define builtin_process_int		&func_builtin_process_int
#define builtin_increase_int		&func_builtin_increase_int
#define builtin_process_long		&func_builtin_process_long
#define builtin_process_mem		&func_builtin_process_mem
#define builtin_process_textflag	&func_builtin_process_textflag
#define builtin_process_string_cpy	&func_builtin_process_string_cpy
#define builtin_process_string_alloc	&func_builtin_process_string_alloc

#define OPT_FLAG(s,v) 	s,					\
			1,PARAM_TEXTFLAG("yes",yesno_flag),	\
			builtin_process_textflag,		\
			v,					\
			OPT_ONCE,				\
			"undocumented flag " # s
#define OPT_MEMPRM(s,v,f) s,				\
			1,PARAM_NUMERAL(NULL),		\
			builtin_process_mem,		\
			v,				\
			OPT_ONCE | (f),			\
			"undocumented flag " # s
#define OPT_INTPRM(s,v,f) s,				\
			1,PARAM_NUMERAL(NULL),		\
			builtin_process_int,		\
			v,				\
			OPT_ONCE | (f),			\
			"undocumented flag " # s
#define OPT_FILENAME(s,v,f) s,				\
			1,PARAM_TEXTUAL(NULL),		\
			builtin_process_string_alloc,	\
			v,				\
			OPT_ONCE | (f),			\
			"undocumented flag " # s

#ifdef	__cplusplus
}
#endif

#endif				/* OPTIONS_PROCESSING_H_ */
