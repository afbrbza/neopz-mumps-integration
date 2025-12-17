# Solução Completa do Problema MUMPS

## Problema Reportado
```
MUMPS solve error: INFO(1) = -3, INFO(2) = 3
MUMPS error -3: MUMPS was called with an invalid value for JOB. JOB = 3.
```

## Causa Raiz (Dupla)

### 1. Matriz Não é Positiva Definida
A equação de Darcy **gera uma matriz simétrica INDEFINIDA**, não positiva definida. 
Quando você chamava `SetDefPositive(true)`, o MUMPS tentava usar:
- **SYM=1** (simétrico positivo definido) 
- **Fatorização de Cholesky (LLT)** sem pivotamento

Mas a matriz precisava de:
- **SYM=2** (simétrico indefinido/geral)
- **Fatorização LDL^T** com pivotamento numérico

Resultado: Erro **INFO(1) = -10** (pivô zero/muito pequeno) durante a fatorização.

### 2. Copy Constructors Inadequados
As classes `TPZMumpsSolver`, `TPZSYsmpMatrixMumps` e `TPZFYsmpMatrixMumps` tinham 
copy constructors `= default`, causando:
- Cópia superficial do estado interno do MUMPS (`DMUMPS_STRUC_C fMumpsData`)
- Quando a fatorização falhava, o estado ficava inconsistente
- Na tentativa de Solve() subsequente, o MUMPS rejeitava JOB=3

## Correções Implementadas

### Correção Principal (Matriz Indefinida)
**Arquivo:** `main.cpp`

**Antes:**
```cpp
an.Assemble();
auto mCast = an.MatrixSolver<STATE>().Matrix();
mCast->SetDefPositive(true);  // ❌ ERRADO para Darcy flow!
an.Solve();
```

**Depois:**
```cpp
an.Assemble();
// Não chamar SetDefPositive() - Darcy produz matriz simétrica INDEFINITE
an.Solve();
```

**Por quê:** A matriz de Darcy flow é simétrica mas indefinida, não positiva definida.

### Correções Estruturais (Copy Constructors)

#### 1. TPZMumpsSolver.h
```cpp
// Antes (ERRADO):
TPZMumpsSolver(const TPZMumpsSolver &copy) = default;
TPZMumpsSolver &operator=(const TPZMumpsSolver &copy) = default;
TPZMumpsSolver(TPZMumpsSolver &&copy) = default;
TPZMumpsSolver &operator=(TPZMumpsSolver &&copy) = default;

// Depois (CORRETO):
TPZMumpsSolver(const TPZMumpsSolver &copy) = delete;  // Não pode copiar estado MUMPS
TPZMumpsSolver &operator=(const TPZMumpsSolver &copy) = delete;
TPZMumpsSolver(TPZMumpsSolver &&copy) noexcept;  // Move é permitido
TPZMumpsSolver &operator=(TPZMumpsSolver &&copy) noexcept;
```

#### 2. TPZMumpsSolver.cpp
Implementados move constructors que:
- Transferem ownership do estado MUMPS
- Marcam o objeto fonte como não-inicializado para evitar double-free
- Clone() agora cria novo solver não-inicializado (não pode copiar estado MUMPS)

#### 3. TPZSYSMPMumps.h/cpp e TPZYSMPMumps.h/cpp
- Copy constructor NÃO copia `fMumpsControl` (cada matriz precisa de sua própria instância)
- Move constructor transfere o `fMumpsControl` corretamente

## Como Identificar se Sua Matriz é Positiva Definida ou Indefinida

### Matrizes Positivas Definidas (usar SYM=1)
- ✅ Problemas de elasticidade linear (pequenas deformações)
- ✅ Condução de calor pura (sem convecção)
- ✅ Poisson com apenas condições de Dirichlet
- ✅ Mecânica estrutural simples

### Matrizes Simétricas Indefinidas (usar SYM=2)
- ✅ **Darcy flow** (equação de pressão-velocidade)
- ✅ Stokes (incompressível)
- ✅ Problemas mistos (pressão-deslocamento)
- ✅ Problemas com pontos de sela (saddle point)
- ✅ Helmholtz
- ✅ Mecânica com rotações finitas

### Matrizes Não-Simétricas (usar SYM=0)
- ✅ Navier-Stokes
- ✅ Convecção-difusão
- ✅ Problemas de transporte

## Resumo das Fases do MUMPS

| JOB | Fase | O que faz |
|-----|------|-----------|
| -1 | INIT | Inicializa estruturas internas |
| 1 | ANALYSIS | Análise simbólica (ordenação, estrutura) |
| 2 | FACTORIZE | Fatorização numérica (LU, LDL^T, ou Cholesky) |
| 3 | SOLVE | Resolve o sistema usando a fatorização |
| -2 | END | Libera memória |

**Importante:** JOB=3 só funciona se INIT, ANALYSIS e FACTORIZE foram bem-sucedidos!

## Teste Final
```bash
cd build && ninja && cd .. && ./build/firsttest
```

**Resultado esperado:** 
- ✅ Sem erros JOB=3
- ✅ Sem erros INFO(1) = -10
- ✅ Arquivo `postproc.0.vtk` gerado com sucesso
- ⚠️ Warnings do MKL podem aparecer (normais, internos do MUMPS)

## Lições Aprendidas

1. **Nunca use `= default` para classes que gerenciam recursos externos** (MUMPS, arquivos, sockets, etc.)

2. **Não assume que sua matriz é positiva definida** - verifique a formulação do problema

3. **Para Darcy flow:** use SYM=2 (indefinido), não SYM=1 (positivo definido)

4. **Verifique INFO(1) após cada fase do MUMPS**, especialmente após FACTORIZE

5. **Erros secundários (JOB=3) podem mascarar o erro real** - sempre verifique a primeira falha
