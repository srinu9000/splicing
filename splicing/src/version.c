
#include "splicing_version.h"
#include <stdio.h>

static const char *splicing_version_string=SPLICING_VERSION;

int splicing_version(const char **version_string,
                   int *major,
                   int *minor,
                   int *subminor) {
  int i1, i2, i3;
  int *p1= major ? major : &i1, 
    *p2= minor ? minor : &i2,
    *p3= subminor ? subminor : &i3;

  if (version_string) { 
    *version_string = splicing_version_string;
  }
  
  *p1 = *p2 = *p3 = 0;
  sscanf(SPLICING_VERSION, "%i.%i.%i", p1, p2, p3);
  
  return 0;
}
