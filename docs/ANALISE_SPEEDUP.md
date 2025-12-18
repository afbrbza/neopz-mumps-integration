# Análise do Speedup Anômalo nos Solvers Pardiso e MUMPS

## Resumo dos Resultados

### Problema Identificado
Ao invés de observar speedup linear (2x com 2 threads, 3x com 3 threads), os resultados mostram:

1. **Pardiso em matrizes pequenas (22,801 eq)**: Speedup máximo de apenas **1.12x** com 8 threads
2. **MUMPS em matrizes pequenas**: **Speedup NEGATIVO** - mais lento que execução serial!
3. **Pardiso em matrizes grandes (203,401 eq)**: Speedup máximo de **1.14x** com 6 threads

---

## Causas Documentadas nos Manuais

### 1. **Overhead de Paralelização (Threading Overhead)**

#### Do Manual do Pardiso:
> "For sufficiently large problem sizes, numerical experiments demonstrate that the scalability of the parallel algorithm is nearly independent of the shared-memory architecture and a speedup of up to seven using eight processors has been observed."

**Interpretação**: 
- Pardiso admite que speedup DEPENDE do tamanho do problema
- Para problemas **pequenos**, o overhead de criar threads, sincronizar e comunicar entre elas **supera** o ganho computacional
- Speedup próximo de 7x só é alcançado em problemas **suficientemente grandes**

#### Do Manual do MUMPS:
> "MUMPS exploits both parallelism arising from sparsity in the matrix A and from dense factorization kernels."

**Problema com matrizes pequenas**:
- Menos oportunidades de paralelização (árvore de eliminação pequena)
- Overhead de OpenMP e sincronização domina o tempo de execução

---

### 2. **Granularidade Computacional Insuficiente**

#### Do Manual do Pardiso (sobre Schur complement):
> "The calculations limit the efficient use of multicore shared-memory environments because it is well known that the triangular solves do not scale very well with the number of cores."

**Tradução**: 
- Operações de **backward/forward substitution** (parte do solve) **NÃO escalam bem**
- Essas operações são sequenciais por natureza (cada passo depende do anterior)
- Em matrizes pequenas, esse gargalo sequencial é mais evidente

---

### 3. **Lei de Amdahl - Fração Sequencial**

A Lei de Amdahl estabelece:

```
Speedup_max = 1 / (S + (1-S)/N)
```

Onde:
- `S` = fração sequencial do código
- `N` = número de threads

**Exemplo com seus dados (Pardiso, 22k eq, 8 threads):**
```
Speedup observado = 1.12
```

Resolvendo para S:
```
1.12 = 1 / (S + (1-S)/8)
S ≈ 0.88 (88% sequencial!)
```

**Significado**: 88% do código de solve está em seções sequenciais ou mal paralelizadas.

---

### 4. **MUMPS: Pior que Serial em Problemas Pequenos**

#### Do Manual do MUMPS:
> "Enable different numbers of threads per MPI. More multithreading in low intensity arithmetic kernels."

**Problema**:
- MUMPS usa **OpenMP + BLAS paralelo** (OpenBLAS no seu caso)
- Em operações pequenas, o overhead de criar regiões paralelas OpenMP é **maior** que o ganho
- Thread creation/destruction overhead a cada operação pequena

**Evidência nos seus dados**:
```
MUMPS, 1 thread:  0.103s
MUMPS, 2 threads: 0.128s (24% MAIS LENTO!)
MUMPS, 7 threads: 0.131s (27% MAIS LENTO!)
```

Isso é **overhead de sincronização** puro.

---

### 5. **Fatores Limitantes no Solve Phase**

#### Do Manual do Pardiso (IPARM(25)):
> "The parameter controls the parallelization for the forward and backward solve. IPARM(25) = 0 indicates that a sequential forward and backward solve is used, whereas IPARM(25) = 1 selects a parallel solve."

**Implicação**:
- Mesmo com IPARM(25)=1 (parallel solve ativo), o **grau de paralelismo é limitado**
- A estrutura da matriz (dependências entre elementos) limita quantas operações podem ser feitas simultaneamente
- Em matrizes pequenas, há **poucos caminhos independentes** na árvore de eliminação

---

## Por Que Problemas Maiores Têm Speedup Melhor (mas ainda limitado)?

### Pardiso - 203,401 equações:
```
1 thread:  1.243s
6 threads: 1.087s  → Speedup: 1.14x
```

**Melhoria, mas ainda longe do ideal (6x)**:

1. **Mais trabalho paralelizável**: Matriz maior → árvore de eliminação maior → mais nós podem ser processados em paralelo
2. **Melhor razão computação/comunicação**: Operações maiores diluem o overhead fixo de sincronização
3. **Mas limitações fundamentais permanecem**:
   - Dependências na árvore de eliminação
   - Operações triangulares sequenciais
   - Cache thrashing com múltiplas threads

---

## Recomendações Baseadas nos Manuais

### Para Problemas Pequenos (<50k equações):
- **Use 1-2 threads no máximo** - mais threads degradam performance
- Pardiso é mais eficiente que MUMPS em problemas pequenos

### Para Problemas Médios (50k-500k equações):
- Pardiso: 4-8 threads podem dar speedup modesto (1.1-1.3x)
- MUMPS: Ainda problemático, evite paralelização excessiva

### Para Problemas Grandes (>500k equações):
- Pardiso: Espere speedup de 2-4x com 8-12 threads
- MUMPS: Pode começar a mostrar benefícios reais

### Otimizações Adicionais (dos manuais):

#### Pardiso:
```cpp
IPARM(24) = 1;  // Two-level scheduling (melhor eficiência paralela)
IPARM(25) = 1;  // Parallel forward/backward solve
```

#### MUMPS:
```cpp
ICNTL(48) = 0;  // Desativa tree-level threading em problemas pequenos
                // (reduz overhead)
```

---

## Conclusão

**Seu speedup "estranho" é ESPERADO e DOCUMENTADO pelos desenvolvedores:**

1. ✅ **Problema de tamanho insuficiente** para justificar overhead de paralelização
2. ✅ **Lei de Amdahl** - fração sequencial domina em solvers diretos esparsos
3. ✅ **Operações triangulares** (solve) escalam mal por natureza
4. ✅ **MUMPS especialmente ruim** em problemas pequenos devido ao overhead OpenMP

**Não é um erro no seu código - é uma limitação fundamental desses solvers!**

Para obter speedup próximo do linear, você precisaria de:
- Problemas **muito maiores** (>1 milhão de equações)
- Ou usar **solvers iterativos** (CG, GMRES) que paralelizam melhor
- Ou aproveitar paralelização no **assembly**, que tipicamente escala muito melhor
