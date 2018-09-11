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
    for (int i = 0; i < 9; ++i)
    {
      os << dg[i] << " ";
    }
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
    os << "\tF=" << pd.deformationGradient << "\n";
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
    os << cs.ss;
    os << "Output\n";
    os << cs.out;
    return os;
  };
  void printMicroFOCase(const MicroCase & cs) { std::cout << cs << std::endl; }
#ifdef HAS_YAML_CPP
  static void operator>>(const YAML::Node & node,
                         DetectOscillationParams & prms)
  {
    std::string tmp;
    node["Type"] >> tmp;
    node["Max Micro Cut Attempt"] >> prms.maxMicroCutAttempts;
    node["Micro Attempt Cut Factor"] >> prms.microAttemptCutFactor;
    if (tmp == "Iteration Only")
    {
      prms.oscType = amsi::DetectOscillationType::IterationOnly;
      node["Max Iterations"] >> prms.maxIterations;
      prms.prevNormFactor = 1.0;
    }
    else if (tmp == "Previous Norm")
    {
      prms.oscType = amsi::DetectOscillationType::PrevNorm;
      node["Previous Norm Factor"] >> prms.prevNormFactor;
      prms.maxIterations = 1;
    }
    else if (tmp == "Combined")
    {
      prms.oscType = amsi::DetectOscillationType::IterationPrevNorm;
      node["Max Iterations"] >> prms.maxIterations;
      node["Previous Norm Factor"] >> prms.prevNormFactor;
    }
    else
    {
      std::cerr << tmp
                << " is not a recognized iteration type. Try 'Iteration Only', "
                   "'Previous Norm', or 'Combined'."
                << std::endl;
      std::abort();
    }
  }
  static void operator>>(const YAML::Node & node, MicroSolutionStrategy & ss)
  {
    std::string tmp;
    node["Solver Type"] >> tmp;
    if (tmp == "Implicit")
    {
      ss.slvrType = SolverType::Implicit;
    }
    else if (tmp == "Explicit")
    {
      ss.slvrType = SolverType::Explicit;
    }
    else if (tmp == "Mixed")
    {
      ss.slvrType = SolverType::Mixed;
    }
    else
    {
      std::cerr << tmp
                << " is not a valid solver type. Try 'Implicit', 'Explicit', "
                   "or 'Mixed'"
                << std::endl;
      std::abort();
    }
    node["Convergence Tolerance"] >> ss.cnvgTolerance;
    // only add get the solver tolerance if it exits
    if (const YAML::Node * sTol = node.FindValue("Solver Tolerance"))
    {
      *sTol >> ss.slvrTolerance;
    }
    else
    {
      std::cout << "Setting the solver tolerance to 1E-6" << std::endl;
      ss.slvrTolerance = 1E-6;
    }
    if (ss.slvrTolerance > ss.cnvgTolerance)
    {
      std::cerr
          << "You should not have the solver tolerance be larger than the "
             "convergence tolerance!"
          << std::endl;
      std::abort();
    }
    DetectOscillationParams oscPrms;
    node["Detect Oscillation"] >> oscPrms;
    ss.oscPrms = oscPrms;
  }
  static void operator>>(const YAML::Node & node, Axis & ax)
  {
    for (int i = 0; i < 3; ++i)
    {
      node[i] >> ax[i];
    }
  }
  static void operator>>(const YAML::Node & node, MicroOutput & out)
  {
    node["3D Orientation Tensor"] >> out.threeDOrntTens;
    node["2D Orientation Tensor"][0] >> out.twoDOrntTens;
    if (out.twoDOrntTens)
    {
      node["2D Orientation Tensor"][1] >> out.twoDOrntTensAxis;
    }
  }
  static void operator>>(const YAML::Node & node, DeformationGradient & dg)
  {
    for (int i = 0; i < 9; ++i)
    {
      node[i] >> dg[i];
    }
  }
  static void operator>>(const YAML::Node & node, MicroProblemDefinition & pd)
  {
    node["Mesh"] >> pd.meshFile;
    node["Deformation Gradient"] >> pd.deformationGradient;
  }
  static void operator>>(const YAML::Node & node, MicroCase & cs)
  {
    node["Case"] >> cs.name;
    MicroProblemDefinition pd;
    MicroSolutionStrategy ss;
    MicroOutput out;
    node["Problem Definition"] >> pd;
    node["Solution Strategy"] >> ss;
    node["Output"] >> out;
    cs.pd = pd;
    cs.ss = ss;
    cs.out = out;
  }
  static void loadMicroFOFromYamlStream(std::istream & fin,
                                        std::vector<MicroCase> & cases)
  {
    try
    {
      YAML::Parser parser(fin);
      YAML::Node doc;
      parser.GetNextDocument(doc);
      for (unsigned i = 0; i < doc.size(); ++i)
      {
        MicroCase cs;
        doc[i] >> cs;
        cases.push_back(cs);
      }
    }
    catch (YAML::ParserException & e)
    {
      std::cout << e.what() << "\n";
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
#endif
}  // namespace bio
