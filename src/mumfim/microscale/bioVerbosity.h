#ifndef MUMFIM_VERBOSITY_H_
#ifdef MUMFIM_VERBOSE_1
#define MUMFIM_V1(out) out
#else
#define MUMFIM_V1(out) while(false);
#endif
#ifdef MUMFIM_VERBOSE_2
#define MUMFIM_V2(out) out
#else
  // need this to compile, but compiler should optimize away
#define MUMFIM_V2(out) while(false);
#endif
#ifdef MUMFIM_VERBOSE_3
#define MUMFIM_V3(out) out
#else
  // need this to compile, but compiler should optimize away
#define MUMFIM_V3(out) while(false);
#endif
#endif
