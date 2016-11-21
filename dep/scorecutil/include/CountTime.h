#ifndef COUNTTIME_H_INCLUDED
#define COUNTTIME_H_INCLUDED

#define CPUFREQ 2300.096e6


#define rdtsc(val)\
 do { \
     unsigned int __a,__d; \
     asm volatile("rdtsc" : "=a" (__a), "=d" (__d)); \
     (val) = ((unsigned long)__a) | (((unsigned long)__d)<<32); \
    } while(0)


#endif

