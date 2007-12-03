#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <fcntl.h>
#include <unistd.h>
#include <string.h>

char buf1[256];
char buf2[256];
char buf3[256];
char buf4[256];

#define RIEN_A_FOUTRE 0
#define FICHIER_C 1
#define FICHIER_H 2

int filenametype(char name[])
{
    int res=RIEN_A_FOUTRE;
    char *v;
    
    for(v=name;(*v!=0)&&(*v!='.');v++);
    
    if ((*v=='.')&&(*(v+2)==0))
	switch(*(v+1))
	{
	    case 'c' : res=FICHIER_C; *v=0; break;
	    case 'h' : res=FICHIER_H; *v=0; break;
	    default: break;
	}
    return res;
}


int exist(const char * s)
{
    int fd;
    fd=open(s,O_RDONLY,0);
    if (fd>0) close(fd);
	return (fd>0);
}

void print_file_dependencies(
	char *dirname, 
	char *filename, 
	char *context, 
	char filetype)
{
    FILE *f;
    char *n,*l;
    char *cbase;
    int e_context,e_common;
    int context_is_common;

    sprintf(buf2,"%s/%s.%c",dirname,filename,filetype);
    f=fopen(buf2,"r");
    if (f==NULL)
    {
	perror(buf2);
	exit(1);
    }
    
    context_is_common=(strcmp(context,"common")==0);
    
    if (filetype=='h')
	cbase="h/";
    else
    {
	cbase="$(OBJDIR)/";
	filetype='o';
    }

    for(;fgets(buf2,256,f)!=NULL;)
    {
	if (strncmp(buf2,"#include",8)!=0) continue;
	if ((n=strchr(buf2,'\"'))==NULL) continue;
	n++;
	l=strstr(n,".h\"");
	if (l==NULL)
	{
	    fprintf(stderr,
		    "Problème dans le fichier %s/%s.%c",
		    dirname,filename,filetype);
	    continue;
	}
	*l=0;

	sprintf(buf4,"h/common/%s.h",n);
	e_common=exist(buf4);
	    
	if (!context_is_common)
	{
	    sprintf(buf4,"h/%s/%s.h",context,n);
	    e_context=exist(buf4);
	    
	    if (!(e_context^e_common))
		fprintf(stderr,"Problème pour le fichier %s "
			"(inclus dans %s), contexte %s\n",
			n,filename,context);
	}
	if (e_common)
	    l="common";
	else
	    l=context;
	    /*l="$(TARGET)"; */

	printf("%s%s/%s.%c	: h/%s/%s.h\n",
		cbase,context, filename, filetype,
		l, n);
    }
    fclose(f);
}

void traite(char *dir,char *name,char *context, char filetype)
{
    int t;
    
    t=filenametype(name);
    
    if (t==RIEN_A_FOUTRE) return;

    if (((t==FICHIER_C)&&(filetype=='c'))||
	((t==FICHIER_H)&&(filetype=='h')))
	print_file_dependencies(dir,name,context,filetype);
    else
	fprintf(stderr,"Que fait %s.%c dans %s ?\n",
		    name,'c'+'h'-filetype,dir);
}

void get_dependencies(char *context, char filetype)
{
    DIR * mydir;
    struct dirent * fichier;
    
    fprintf(stderr,"Examining [%s] (level '%c')\n",context,filetype);

    printf("######################################"
	   "######################################\n");
    printf("# DEPENDANCES SUR LES .%c, au niveau %s\n\n",filetype,context);
    
    sprintf(buf1,"%c/%s",filetype,context);
    
    mydir = opendir(buf1);
    
    if (mydir==NULL)
    {
	if (errno == ENOENT)
	{
	    printf("# (Aucune, répertoire inexistant)\n");
	    return;
	}
	else
	{
	    perror(buf1);
	    exit(errno);
	}
    }
    
    for(rewinddir(mydir);(fichier=readdir(mydir))!=NULL;)
    {
	strcpy(buf3,fichier->d_name);
	traite(buf1,buf3,context,filetype);
    }

    closedir(mydir);
    
    if ((strcmp(context,"common")==0)&&(filetype=='h'))
    {
	sprintf(buf1,"%c/%s/def",filetype,context);
	mydir = opendir(buf1);
	if (mydir!=NULL)
	{
	    sprintf(buf1,"%c/%s",filetype,context);
	
	    for(rewinddir(mydir);(fichier=readdir(mydir))!=NULL;)
	    {
		sprintf(buf3,"def/%s",fichier->d_name);
		traite(buf1,buf3,context,filetype);
	    }
	    closedir(mydir);
	}
	else if (errno != ENOENT)
	    perror(buf1);	/* c'est pas la peine de virer, cela dit\n" */
    }
    printf("\n\n");
}

int main(int argc, char *argv[])
{
    int i;

    for(i=1;i<argc;i++)
    {
	get_dependencies(argv[i],'c');
	get_dependencies(argv[i],'h');
    }
    
    return 0;
}


