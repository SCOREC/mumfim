#include "yaml-cpp/node/parse.h"
#ifdef HAS_YAML_CPP
#include <yaml-cpp/yaml.h>
#endif
#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include "bioMicroFOParams.h"
namespace bio
{
  static std::ostream & operator<<(std::ostream & os,
                                   const DeformationGradient & dg)
  {
    os << dg[0] << " " << dg[1] << " " << dg[2] << "; ";
    os << dg[3] << " " << dg[4] << " " << dg[5] << "; ";
    os << dg[6] << " " << dg[7] << " " << dg[8];
    return os;
  }
  static std::ostream & operator<<(std::ostream & os, const Axis & ax)
  {
    for (int i = 0; i < 3; ++i)
    {
      os << ax[i] << " ";
    }
    return os;
  }
  static std::ostream & operator<<(std::ostream & os,
                                   const MicroProblemDefinition & pd)
  {
    os << "\t" << pd.meshFile << "\n";
    os << "\tF=[" << pd.deformationGradient << "]\n";
    return os;
  }
  static std::ostream & operator<<(std::ostream & os,
                                   const DetectOscillationParams & dtctOscPrms)
  {
    os << "\tDetect Oscillation Type: ";
    switch (dtctOscPrms.oscType)
    {
      case amsi::DetectOscillationType::IterationOnly:
        os << "Iteration Only\n";
        break;
      case amsi::DetectOscillationType::PrevNorm:
        os << "Previous Norm\n";
        break;
      case amsi::DetectOscillationType::IterationPrevNorm:
        os << "Combined\n";
        break;
    }
    os << "\tMax Cut Attempts: " << dtctOscPrms.maxMicroCutAttempts << "\n";
    os << "\tAttempt Cut Factor: " << dtctOscPrms.microAttemptCutFactor << "\n";
    os << "\tMax Iterations: " << dtctOscPrms.maxIterations << "\n";
    os << "\tPrevious Norm Factor: " << dtctOscPrms.prevNormFactor << "\n";
    return os;
  }
  static std::ostream & operator<<(std::ostream & os,
                                   const MicroSolutionStrategy & ss)
  {
    os << "\tConvergence Tolerance: " << ss.cnvgTolerance << "\n";
    os << "\tSolver Tolerance: " << ss.slvrTolerance << "\n";
    os << "\tSolver Type: ";
    switch (ss.slvrType)
    {
      case SolverType::Implicit:
        os << "Implicit\n";
        break;
      case SolverType::Explicit:
        os << "Explicit\n";
        os << "\t\tLoad Time: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss).total_time
           << "\n";
        os << "\t\tAmplitude Type: ";
        if (static_cast<const MicroSolutionStrategyExplicit &>(ss).ampType ==
            AmplitudeType::SmoothStep)
        {
          os << "Smooth Step\n";
        }
        os << "\t\tViscous Damping Factor: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss)
                  .visc_damp_coeff
           << "\n";
        os << "\t\tPrint History Frequency: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss)
                  .print_history_frequency
           << "\n";
        os << "\t\tPrint Field Frequency: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss)
                  .print_field_frequency
           << "\n";
        os << "\t\tPrint Field By Num Frames: "
           << (static_cast<const MicroSolutionStrategyExplicit &>(ss)
                       .print_field_by_num_frames
                   ? "true"
                   : "false")
           << "\n";
        os << "\t\tCritical Time Scale Factor: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss)
                  .crit_time_scale_factor
           << "\n";
        os << "\t\tEnergy Check Epsilon: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss)
                  .energy_check_eps
           << "\n";
        os << "\t\tSerial GPU Cutoff: "
           << static_cast<const MicroSolutionStrategyExplicit &>(ss)
                  .serial_gpu_cutoff
           << "\n";
        break;
      case SolverType::Mixed:
        os << "Mixed\n";
        break;
    }
    os << ss.oscPrms;
    return os;
  }
  static std::ostream & operator<<(std::ostream & os, const MicroOutput & out)
  {
    os << "\tOutput 3D Orientation Tensor: "
       << (out.threeDOrntTens ? "true\n" : "false\n");
    os << "\tOutput 2D Orientation Tensor: "
       << (out.twoDOrntTens ? "true " : "false");
    if (out.twoDOrntTens)
    {
      os << "[ " << out.twoDOrntTensAxis << "]\n";
    }
    return os;
  }
  static std::ostream & operator<<(std::ostream & os, const MicroCase & cs)
  {
    os << "Case: " << cs.name << "\n";
    os << "Problem Definition\n";
    os << cs.pd;
    os << "Solution Strategy\n";
    os << *cs.ss;
    os << "Output\n";
    os << cs.out;
    return os;
  };
  void printMicroFOCase(const MicroCase & cs) { std::cout << cs << std::endl; }
}  // namespace bio
#ifdef HAS_YAML_CPP
namespace YAML
{
  template <>
  struct convert<bio::DetectOscillationParams>
  {
    static bool decode(const Node & node, bio::DetectOscillationParams & prms)
    {
      auto type = node["Type"].as<std::string>();
      prms.maxMicroCutAttempts = node["Max Micro Cut Attempt"].as<int>();
      prms.microAttemptCutFactor = node["Micro Attempt Cut Factor"].as<int>();
      if (type == "Iteration Only")
      {
        prms.oscType = amsi::DetectOscillationType::IterationOnly;
        prms.maxIterations = node["Max Iterations"].as<int>();
        prms.prevNormFactor = 1.0;
      }
      else if (type == "Previous Norm")
      {
        prms.oscType = amsi::DetectOscillationType::PrevNorm;
        prms.prevNormFactor = node["Previous Norm Factor"].as<double>();
        prms.maxIterations = 1;
      }
      else if (type == "Combined")
      {
        prms.oscType = amsi::DetectOscillationType::IterationPrevNorm;
        prms.maxIterations = node["Max Iterations"].as<int>();
        prms.prevNormFactor = node["Previous Norm Factor"].as<double>();
      }
      else
      {
        std::cerr
            << type
            << " is not a recognized iteration type. Try 'Iteration Only',"
               "'Previous Norm', or 'Combined'."
            << std::endl;
        std::abort();
      }
      return true;
    }
  };
  template <>
  struct convert<std::unique_ptr<bio::MicroSolutionStrategy>>
  {
    static bool decode(const Node & node,
                       std::unique_ptr<bio::MicroSolutionStrategy> & ss)
    {
      auto solver_type = node["Solver Type"].as<std::string>();
      if (solver_type == "Implicit")
      {
        ss = std::make_unique<bio::MicroSolutionStrategy>();
        ss->slvrType = bio::SolverType::Implicit;
      }
      else if (solver_type == "Explicit")
      {
        auto es = std::make_unique<bio::MicroSolutionStrategyExplicit>();
        es->slvrType = bio::SolverType::Explicit;
        es->print_history_frequency = node["Print History Frequency"].as<int>();
        es->print_field_frequency = node["Print Field Frequency"].as<int>();
        es->visc_damp_coeff = node["Viscous Damping Factor"].as<double>();
        es->total_time = node["Load Time"].as<double>();
        es->serial_gpu_cutoff = node["Serial GPU Cutoff"].as<int>();
        es->crit_time_scale_factor =
            node["Critical Time Scale Factor"].as<double>();
        es->energy_check_eps = node["Energy Check Epsilon"].as<double>();
        auto amplitude_type = node["Amplitude Type"].as<std::string>();
        if (amplitude_type == "Smooth Step")
        {
          es->ampType = bio::AmplitudeType::SmoothStep;
        }
        else
        {
          std::cerr << amplitude_type
                    << " is not a valid amplitude type. Try 'Smooth Step'.\n";
          std::abort();
        }
        ss = std::move(es);
      }
      else if (solver_type == "Mixed")
      {
        ss->slvrType = bio::SolverType::Mixed;
        std::cerr << "Mixed Solver Not Currently Implemented" << std::endl;
        std::abort();
      }
      else
      {
        std::cerr << solver_type
                  << " is not a valid solver type. Try 'Implicit', 'Explicit', "
                     "or 'Mixed'"
                  << std::endl;
        std::abort();
      }
      ss->cnvgTolerance = node["Convergence Tolerance"].as<double>();
      // only add get the solver tolerance if it exits
      if (node["Solver Tolerance"])
      {
        ss->slvrTolerance = node["Solver Tolerance"].as<double>();
      }
      else
      {
        std::cout << "Setting the solver tolerance to 1E-6" << std::endl;
        ss->slvrTolerance = 1E-6;
      }
      if (ss->slvrTolerance > ss->cnvgTolerance)
      {
        std::cerr
            << "You should not have the solver tolerance be larger than the "
               "convergence tolerance!"
            << std::endl;
        std::abort();
      }
      ss->oscPrms =
          node["Detect Oscillation"].as<bio::DetectOscillationParams>();
      return true;
    }
  };
  template <>
  struct convert<bio::MicroOutput>
  {
    static bool decode(const Node & node, bio::MicroOutput & out)
    {
      out.threeDOrntTens = node["3D Orientation Tensor"].as<bool>();
      out.twoDOrntTens = node["2D Orientation Tensor"][0].as<bool>();
      if (out.twoDOrntTens)
      {
        auto axis = node["2D Orientation Tensor"][1];
        if (!axis.IsSequence() || axis.size() != 3)
        {
          return false;
        }
        for (int i = 0; i < 3; ++i)
        {
          out.twoDOrntTensAxis[i] = axis[i].as<double>();
        }
      }
      return true;
    }
  };
  template <>
  struct convert<bio::MicroProblemDefinition>
  {
    static bool decode(const Node & node, bio::MicroProblemDefinition & pd)
    {
      pd.meshFile = node["Mesh"].as<std::string>();
      auto dg = node["Deformation Gradient"];
      if (!dg.IsSequence() || dg.size() != 9)
      {
        return false;
      }
      for (size_t i = 0; i < dg.size(); ++i)
      {
        pd.deformationGradient[i] = dg[i].as<double>();
      }
      return true;
    }
  };
  template <>
  struct convert<bio::MicroCase>
  {
    static bool decode(const Node & node, bio::MicroCase & cs)
    {
      cs.name = node["Case"].as<std::string>();
      cs.pd = node["Problem Definition"].as<bio::MicroProblemDefinition>();
      cs.out = node["Output"].as<bio::MicroOutput>();
      cs.ss = node["Solution Strategy"]
                  .as<std::unique_ptr<bio::MicroSolutionStrategy>>();
      return true;
    }
  };
}  // namespace YAML
namespace bio
{
  static void loadMicroFOFromYamlStream(std::istream & fin,
                                        std::vector<MicroCase> & cases)
  {
    try
    {
      auto doc = YAML::Load(fin);
      cases = doc.as<std::vector<MicroCase>>();
    }
    catch (YAML::ParserException & e)
    {
      std::cerr << e.what() << "\n";
    }
  }
  void loadMicroFOFromYamlFile(const std::string & filename,
                               std::vector<MicroCase> & cases)
  {
    std::ifstream fin(filename);
    if (!fin.is_open())
    {
      std::cerr << "Could not open file: " << filename << " for reading"
                << std::endl;
      std::abort();
    }
    loadMicroFOFromYamlStream(fin, cases);
  }
}  // namespace bio
#endif
