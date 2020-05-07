#ifndef BIO_BATCHED_ANALYSIS_H__
#define BIO_BATCHED_ANALYSIS_H__

/*
 * The assumption is that we will not want a list of Batched run objects,
 * so we use template specializations of the internal methods to obtain
 * different behavior for different sizes of networks. This is in contrast
 * to using dynamic polymorphism
 */


class BatchedRun
{
};

#endif
