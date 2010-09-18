#include <stdio.h>
#include <unistd.h>

int fileaccess_( char *fname ,int *length) 
{
  fname[*length] = '\0';
  return(access(fname,R_OK));
}
