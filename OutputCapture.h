#ifndef OUTPUTCAPTURE_H
#define OUTPUTCAPTURE_H

#include <sstream>
#include <string>

// Estrutura para armazenar métricas do solver
struct SolverMetrics {
  // Tempos
  double analysisTime = 0.0;
  double factorizationTime = 0.0;
  double solveTime = 0.0;
  double totalTime = 0.0;

  // Operações
  double estimatedFlops = 0.0;
  double actualFlops = 0.0;

  // Memória
  int realSpaceFactors = 0;
  int integerSpaceFactors = 0;
  int maxFrontalSize = 0;

  // Threads
  int numThreads = 0;

  // Iterações (para PARDISO com CGS)
  int cgsIterations = 0;

  // Ordering based on the output (e.g., "METIS", "AMD", etc.)
  std::string orderingBasedOn;
};

class OutputCapture {
private:
  std::stringstream buffer;
  int saved_stdout;
  int saved_stderr;
  int pipe_fd[2];

public:
  OutputCapture();
  ~OutputCapture();
  std::string getOutput();
  SolverMetrics parseMumpsOutput(const std::string& output);
  SolverMetrics parsePardisoOutput(const std::string& output);
};

#endif // OUTPUTCAPTURE_H
