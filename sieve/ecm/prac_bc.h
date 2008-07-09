/* One more than the highest code number the byte code generator can produce */
#define MAXCODE 32
#define PRAC_NR_MULTIPLIERS 10

typedef char literal_t;

void bytecoder_init (int);
void bytecoder_clear ();
void bytecoder (const int c);
void bytecoder_flush ();

/* Returns the number of bytes currently in the bytecoder buffer */
unsigned int bytecoder_size ();
/* Writes all the data currently in the bytecoder buffer to the given pointer,
   and clears the buffer */
void bytecoder_read (char *);

unsigned long prac_best (double *, const unsigned long, const int, 
			 const unsigned int, const unsigned int);
void prac_bytecode (const unsigned long, const unsigned int, 
		    const unsigned int);
 
