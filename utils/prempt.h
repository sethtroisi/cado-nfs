/* Max lenght of the line of the current shell */ 
#define PREMPT_S_CMD (1<<16)

/* Lenght of prempt buffer. Must be a power of 2. */
#define PREMPT_BUF (1<<22)

/* Lenght of one write in prempt buffer. Between 64 and 1024 Ko
   seems the best. */
#define PREMPT_ONE_READ (PREMPT_BUF>>3)

/* #define MIN(A,B) ((A)<(B)?(A):(B)) */

typedef struct {
  char **files, *buf;
  volatile char *pcons, *pprod;
  volatile unsigned int end;
} __prempt_t;
typedef __prempt_t prempt_t[1];

/* Return a unix commands list. Exemple :
   cat file_relation1
   gzip -dc file_relation2.gz file_relation3.gz
   bzip2 -dc file_relation4.gz file_relation5.gz
   [empty string]
*/
extern char ** prempt_open_compressed_rs (char ** ficname);
