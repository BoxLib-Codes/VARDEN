#include <stdio.h>

int
main()
{
    printf("FLOAT TYPES\n");
    printf( "sizeof(long double) = %ld\n", sizeof(long double));
    printf( "sizeof(double)      = %ld\n", sizeof(double));
    printf( "sizeof(float)       = %ld\n", sizeof(float));
    printf( "INT   TYPES\n");
    printf( "sizeof(long long)   = %ld\n", sizeof(long long));
    printf( "sizeof(long)        = %ld\n", sizeof(long));
    printf( "sizeof(int)         = %ld\n", sizeof(int));
    printf( "sizeof(short)       = %ld\n", sizeof(short));
    printf( "sizeof(char)        = %ld\n", sizeof(char));
    return(0);
}
