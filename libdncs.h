#ifndef LIBDNCS_H
#define LIBDNCS_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

int sample(const char* seq, const char* folder, const char* ff, size_t sample, const char* method);
char* pdb_to_angle(const char* filename);
void free_string(char* str);

#ifdef __cplusplus
}
#endif

#endif // LIBDNCS_H
