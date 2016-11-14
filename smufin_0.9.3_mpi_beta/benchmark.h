#include <chrono>
#include <iostream>
#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <time.h>

#define START(x) x##start = std::chrono::steady_clock::now()                                         
#define END(x)   x##end   = std::chrono::steady_clock::now()                                         
#define TIME(x)  x##time  = std::chrono::duration_cast<std::chrono::nanoseconds>(x##end-x##start).count()
#define PRINT(x) std::cout << #x << "(), " << x##time << std::endl

#define CSTART(x) \
    struct timespec x##start;\
    clock_gettime(CLOCK_REALTIME, &x##start); \

#define CEND(x)\
    struct timespec x##end;\
    clock_gettime(CLOCK_REALTIME, &x##end); 

#define CTIME(x) double x##time = (x##end.tv_nsec - x##start.tv_nsec);

#define CPRINT(x) printf(#x); printf("(), %f\n", x##time);
