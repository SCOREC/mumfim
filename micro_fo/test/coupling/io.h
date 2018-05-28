#ifndef BIO_TESTING_IO_H_
#define BIO_TESTING_IO_H_
#include "apfDynamicMatrix.h"
#include <fstream>
std::istream & operator>>(std::istream & in,
                          apf::DynamicMatrix & m);
std::istream & operator>>(std::istream & in,
                           apf::DynamicVector & v);
bool operator==(const apf::DynamicMatrix & a,
                const apf::DynamicMatrix & b);
void rowPermute(const apf::DynamicMatrix & a,
                int * prmt,
                apf::DynamicMatrix & b);
#endif
