#ifndef BIO_VERBOSITY_H_
#ifdef BIO_VERBOSE_1
#define BIO_V1(out) out
#else
#define BIO_V1(out) while(false);
#endif
#ifdef BIO_VERBOSE_2
#define BIO_V2(out) out
#else
  // need this to compile, but compiler should optimize away
#define BIO_V2(out) while(false);
#endif
#ifdef BIO_VERBOSE_3
#define BIO_V3(out) out
#else
  // need this to compile, but compiler should optimize away
#define BIO_V3(out) while(false);
#endif
#endif
