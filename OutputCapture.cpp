#include "OutputCapture.h"
#include <regex>
#include <iostream>
#include <cmath>
#include <unistd.h>
#include <fcntl.h>

// Implementação da classe OutputCapture
OutputCapture::OutputCapture() : saved_stdout(-1), saved_stderr(-1) {
    // Criar pipe para capturar saída
    if (pipe(pipe_fd) == -1) {
        return;
    }
    
    // Salvar file descriptors originais
    saved_stdout = dup(STDOUT_FILENO);
    saved_stderr = dup(STDERR_FILENO);
    
    // Redirecionar stdout e stderr para o pipe
    dup2(pipe_fd[1], STDOUT_FILENO);
    dup2(pipe_fd[1], STDERR_FILENO);
    
    // Tornar a leitura do pipe não-bloqueante
    int flags = fcntl(pipe_fd[0], F_GETFL, 0);
    fcntl(pipe_fd[0], F_SETFL, flags | O_NONBLOCK);
}

OutputCapture::~OutputCapture() {
    if (saved_stdout != -1) {
        // Restaurar stdout e stderr
        fflush(stdout);
        fflush(stderr);
        dup2(saved_stdout, STDOUT_FILENO);
        dup2(saved_stderr, STDERR_FILENO);
        close(saved_stdout);
        close(saved_stderr);
    }
    
    if (pipe_fd[0] != -1) {
        close(pipe_fd[0]);
    }
    if (pipe_fd[1] != -1) {
        close(pipe_fd[1]);
    }
}

std::string OutputCapture::getOutput() {
    // Fazer flush antes de ler
    fflush(stdout);
    fflush(stderr);
    
    char buf[4096];
    ssize_t n;
    while ((n = read(pipe_fd[0], buf, sizeof(buf))) > 0) {
        buffer.write(buf, n);
    }
    
    return buffer.str();
}

// Parser para saída do MUMPS
SolverMetrics OutputCapture::parseMumpsOutput(const std::string& output) {
    SolverMetrics metrics;
    std::regex regex;
    std::smatch match;
    
    // Extrair número de threads
    regex = std::regex(R"(executing #MPI\s*=\s*\d+\s+and #OMP\s*=\s*(\d+))");
    if (std::regex_search(output, match, regex)) {
        metrics.numThreads = std::stoi(match[1]);
    }
    
    // Extrair tempo de análise
    regex = std::regex(R"(Elapsed time in analysis driver\s*=\s*([\d.]+))");
    if (std::regex_search(output, match, regex)) {
        metrics.analysisTime = std::stod(match[1]);
    }
    
    // Extrair tempo de fatoração
    regex = std::regex(R"(Elapsed time for factorization\s*=\s*([\d.]+))");
    if (std::regex_search(output, match, regex)) {
        metrics.factorizationTime = std::stod(match[1]);
    }
    
    // Extrair tempo total de solve
    regex = std::regex(R"(Elapsed time in solve driver\s*=\s*([\d.]+))");
    if (std::regex_search(output, match, regex)) {
        metrics.solveTime = std::stod(match[1]);
    }
    
    // Extrair operações estimadas (em notação científica)
    regex = std::regex(R"(RINFOG\(1\) Operations during elimination \(estim\)\s*=\s*([\d.]+)D\+(\d+))");
    if (std::regex_search(output, match, regex)) {
        double base = std::stod(match[1]);
        int exp = std::stoi(match[2]);
        metrics.estimatedFlops = base * std::pow(10, exp);
    }
    
    // Extrair operações reais
    regex = std::regex(R"(Operations in node elimination\s*=\s*([\d.]+)D\+(\d+))");
    if (std::regex_search(output, match, regex)) {
        double base = std::stod(match[1]);
        int exp = std::stoi(match[2]);
        metrics.actualFlops = base * std::pow(10, exp);
    }
    
    // Extrair espaço real para fatores
    regex = std::regex(R"(INFOG\s*\(9\) Real space for factors\s*=\s*(\d+))");
    if (std::regex_search(output, match, regex)) {
        metrics.realSpaceFactors = std::stoi(match[1]);
    }
    
    // Extrair espaço inteiro para fatores
    regex = std::regex(R"(INFOG\s*\(10\) Integer space for factors\s*=\s*(\d+))");
    if (std::regex_search(output, match, regex)) {
        metrics.integerSpaceFactors = std::stoi(match[1]);
    }
    
    // Extrair tamanho máximo frontal
    regex = std::regex(R"(INFOG\s*\(11\) Maximum front size\s*=\s*(\d+))");
    if (std::regex_search(output, match, regex)) {
        metrics.maxFrontalSize = std::stoi(match[1]);
    }
     
    // Extrair o Ordering utilizado (MUMPS is "Ordering based on METIS" ou similar)
    regex = std::regex(R"(Ordering based on\s+(\w+))");
    if (std::regex_search(output, match, regex)) {
        metrics.orderingBasedOn = match[1];
    }

    
    metrics.totalTime = metrics.analysisTime + metrics.factorizationTime + metrics.solveTime;
    
    return metrics;
}

// Parser para saída do PARDISO
SolverMetrics OutputCapture::parsePardisoOutput(const std::string& output) {
    SolverMetrics metrics;
    std::regex regex;
    std::smatch match;
    
    // Extrair número de threads
    regex = std::regex(R"(Parallel Direct Factorization is running on (\d+) OpenMP)");
    if (std::regex_search(output, match, regex)) {
        metrics.numThreads = std::stoi(match[1]);
    }
    
    // Extrair tempo de reordenamento
    regex = std::regex(R"(Time spent in reordering.*:\s*([\d.]+)\s*s)");
    double reorderTime = 0.0;
    if (std::regex_search(output, match, regex)) {
        reorderTime = std::stod(match[1]);
    }
    
    // Extrair tempo de fatoração simbólica
    regex = std::regex(R"(Time spent in symbolic factorization.*:\s*([\d.]+)\s*s)");
    double symbolicTime = 0.0;
    if (std::regex_search(output, match, regex)) {
        symbolicTime = std::stod(match[1]);
    }
    
    metrics.analysisTime = reorderTime + symbolicTime;
    
    // Extrair tempo de fatoração numérica
    regex = std::regex(R"(Time spent in factorization step.*:\s*([\d.]+)\s*s)");
    if (std::regex_search(output, match, regex)) {
        metrics.factorizationTime = std::stod(match[1]);
    }
    
    // Extrair tempo de solve (CGS/CG)
    regex = std::regex(R"(Time spent in iterative solver.*:\s*([\d.]+)\s*s)");
    if (std::regex_search(output, match, regex)) {
        metrics.solveTime = std::stod(match[1]);
        
        // Extrair número de iterações CGS
        std::regex iterRegex(R"(cgx iterations (\d+))");
        if (std::regex_search(output, match, iterRegex)) {
            metrics.cgsIterations = std::stoi(match[1]);
        }
    }
    
    // Extrair GFlops
    regex = std::regex(R"(gflop\s+for the numerical factorization:\s*([\d.]+))");
    if (std::regex_search(output, match, regex)) {
        metrics.actualFlops = std::stod(match[1]) * 1e9; // Converter GFlop para Flop
    }
    
    // Extrair não-zeros em L
    regex = std::regex(R"(number of non-zeros in L:\s*(\d+))");
    if (std::regex_search(output, match, regex)) {
        metrics.realSpaceFactors = std::stoi(match[1]);
    }
    
    // Extrair tamanho do maior supernode
    regex = std::regex(R"(size of largest supernode:\s*(\d+))");
    if (std::regex_search(output, match, regex)) {
        metrics.maxFrontalSize = std::stoi(match[1]);
    }
    
    metrics.totalTime = metrics.analysisTime + metrics.factorizationTime + metrics.solveTime;
    
    return metrics;
}
