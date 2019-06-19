#ifndef BIO_MICRO_FO_PARAMS_H_
#define BIO_MICRO_FO_PARAMS_H_
#include <amsiDetectOscillation.h>
#include <string>
namespace bio
{
  struct MicroCase;
  void printMicroFOCase(const MicroCase & cs);
#ifdef HAS_YAML_CPP
  void loadMicroFOFromYamlFile(const std::string & filename, std::vector<MicroCase> & cs);
#endif
  enum class SolverType
  {
    Implicit,
    Explicit,
    Mixed
  };
  struct DeformationGradient
  {
    double data[9];
    double & operator[](std::size_t idx) { return data[idx]; }
    const double & operator[](std::size_t idx) const { return data[idx]; }
    DeformationGradient operator*(double val)
    {
      DeformationGradient dg;
      for (int i = 0; i < 9; ++i)
      {
        dg.data[i] = data[i] * val;
      }
      return dg;
    }
  };
  struct Axis
  {
    double data[3];
    double & operator[](std::size_t idx) { return data[idx]; }
    const double & operator[](std::size_t idx) const { return data[idx]; }
  };
  struct MicroProblemDefinition
  {
    std::string meshFile;
    DeformationGradient deformationGradient;
  };
  struct DetectOscillationParams
  {
    amsi::DetectOscillationType oscType;
    int maxMicroCutAttempts;
    int microAttemptCutFactor;
    // params for iteration type
    int maxIterations;
    // params for prev norm type
    double prevNormFactor;
  };
  struct MicroSolutionStrategy
  {
    double cnvgTolerance;
    double slvrTolerance;
    SolverType slvrType;
    DetectOscillationParams oscPrms;
  };
  struct MicroOutput
  {
    bool threeDOrntTens;
    bool twoDOrntTens;
    Axis twoDOrntTensAxis;
  };
  struct MicroCase
  {
    std::string name;
    MicroProblemDefinition pd;
    MicroSolutionStrategy ss;
    MicroOutput out;
  };
}  // namespace bio
#endif
