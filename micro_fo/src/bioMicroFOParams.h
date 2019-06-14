#ifndef BIO_MICRO_FO_PARAMS_H_
#define BIO_MICRO_FO_PARAMS_H_
#include <amsiDetectOscillation.h>
#include <string>
#include <memory>
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
  enum class AmplitudeType
  {
    SmoothStep,
    SmoothStepHold
  };
  struct DeformationGradient
  {
    double data[9];
    double & operator[](std::size_t idx) { return data[idx]; }
    const double & operator[](std::size_t idx) const { return data[idx]; }
    DeformationGradient(double d11, double d12, double d13, double d21, double d22, double d23, double d31, double d32, double d33)
    {
      data[0] = d11;
      data[1] = d12;
      data[2] = d13;
      data[3] = d21;
      data[4] = d22;
      data[5] = d23;
      data[6] = d31;
      data[7] = d32;
      data[8] = d33;
    }
    DeformationGradient() {};
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
  // currently we have these parameters here because
  // they are needed for computing derivative terms
  struct MicroSolutionStrategy
  {
    double cnvgTolerance;
    double slvrTolerance;
    SolverType slvrType;
    DetectOscillationParams oscPrms;
  };
  struct MicroSolutionStrategyExplicit : public MicroSolutionStrategy
  {
    unsigned int print_history_frequency;
    unsigned int print_field_frequency;
    bool print_field_by_num_frames;
    double visc_damp_coeff;
    double total_time;
    int serial_gpu_cutoff;
    double crit_time_scale_factor;
    double energy_check_eps;
    AmplitudeType ampType;
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
    std::unique_ptr<MicroSolutionStrategy> ss;
    MicroOutput out;
    MicroCase() : ss(new MicroSolutionStrategy) {}
  };
}  // namespace bio
#endif
