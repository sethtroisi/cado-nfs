/* One more than the highest code number the byte code generator can produce */
#define BC_MAXCODE 32

typedef char literal_t;
typedef char code_t;

/* A dictionary. It contains nr_entries entries, each entry consisting of
   a sequence of literals (a key) and a code. 
   Keys are translated to codes greedily (i.e., at any point we try to get 
   the longest dictionary match before writing a code).
   A code of 0 in the dictionary is special and means: write nothing */
typedef struct {
  int nr_entries;
  size_t *key_len;     /* The lengths of the literal sequences */
  literal_t **key;
  code_t *code;
} bc_dict_t;

typedef struct {
  size_t histlen, nrstored; /* Max window size, nr of literals in history */
  literal_t *history;       /* Current history window */
  code_t *buffer;           /* The codes written so far */
  size_t bufalloc, buffull; /* Allocated number of codes and 
			       current number of codes in buffer */
  const bc_dict_t *dict;    /* Pointer to the dictionary we use */
} bc_state_t;

bc_state_t * bytecoder_init (const bc_dict_t *);
void bytecoder_clear (bc_state_t *);
void bytecoder (const literal_t, bc_state_t *);
void bytecoder_flush (bc_state_t *state);

/* Returns the number of bytes currently in the bytecoder buffer */
size_t bytecoder_size (const bc_state_t *);
/* Writes all the data currently in the bytecoder buffer to the given pointer,
   and clears the buffer */
void bytecoder_read (code_t *, bc_state_t *);

void prac_bytecode (const unsigned long, const double, const double, 
		    const double, const double, bc_state_t *);
