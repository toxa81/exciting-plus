#include <string.h>

void cmemcpy_(void* dest, const void* src, int* length)
{
    memcpy(dest, src, *length);
}
