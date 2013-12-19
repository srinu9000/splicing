#ifndef SPLICING_VERSION_H
#define SPLICING_VERSION_H

#define SPLICING_VERSION "@VERSION@"

int splicing_version(const char **version_string,
		     int *major,
		     int *minor,
		     int *subminor);

#endif
