# Como Usar os Solvers MUMPS, Pardiso e Skyline

## Escolhendo o Solver

No arquivo `main.cpp`, linha ~209, altere a variável `solverType`:

```cpp
const EnumSolvers solverType = EMumps;     // MUMPS solver
// const EnumSolvers solverType = EPardiso; // Pardiso solver  
// const EnumSolvers solverType = ESkyline; // Skyline solver
```

## Resultados Esperados

Todos os três solvers devem produzir a **mesma solução**:

```
Solution: [1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 3.0, 1.0, 2.0]
```

## Diferenças Entre os Solvers

### Skyline
- **Formato**: Skyline storage (faixas)
- **Simetria**: Explora simetria
- **Memória**: Armazena apenas a banda
- **Velocidade**: Rápido para matrizes pequenas
- **Ideal para**: Problemas pequenos a médios

### Pardiso (Intel MKL)
- **Formato**: CSR simétrico (apenas triangular superior)
- **Simetria**: Explora simetria
- **Memória**: Armazena apenas parte triangular
- **Velocidade**: Muito rápido (otimizado MKL)
- **Ideal para**: Todos os tamanhos de problema
- **Requisito**: Intel MKL instalado

### MUMPS
- **Formato**: COO não-simétrico (matriz completa)
- **Simetria**: Não explora (devido ao bug da versão simétrica)
- **Memória**: ~1.5x mais que Pardiso
- **Velocidade**: Comparável ao Pardiso
- **Ideal para**: Quando Pardiso não está disponível
- **Vantagem**: Código aberto, não requer MKL

## ⚠️ Atenção: Bug na Versão Simétrica do MUMPS

**NÃO USAR**: `TPZSSpStructMatrixMumps` (versão simétrica)

Esta versão tem um bug que produz soluções incorretas:
```
❌ [1.71, 5.25, 0, 1.75, 5.14, 0, 4.5, 1.5, 0]  // ERRADO!
```

**USAR**: `TPZSpStructMatrixMumps` (versão não-simétrica)

Esta versão funciona perfeitamente:
```
✅ [1.0, 3.0, 2.0, 1.0, 3.0, 2.0, 3.0, 1.0, 2.0]  // CORRETO!
```

**O código já está configurado corretamente** para usar a versão não-simétrica.

## Compilação e Execução

```bash
cd build
cmake --build . -j4
./firsttest
```

## Verificando a Solução

Procure por esta saída no terminal:

```
Writing matrix 'Solution' (9 x 1):
	1.00000000001  
	2.99999999999  
	1.99999999999  
	1.00000000001  
	2.99999999999  
	2.0  
	3.0  
	1.00000000001  
	2.0  
```

Os pequenos erros numéricos (e.g., 1.00000000001 em vez de 1.0) são **normais** e devidos à precisão de ponto flutuante.

## Troubleshooting

### Erro: "MUMPS factorization error: INFO(1) = -10"

**Causa**: Usando versão simétrica (`TPZSSpStructMatrixMumps`)  
**Solução**: Verifique que `main.cpp` linha ~214 usa `TPZSpStructMatrixMumps`

### Erro: "Intel oneMKL ERROR: Parameter 11"

**Causa**: Bug na versão simétrica do MUMPS  
**Solução**: Use versão não-simétrica (já está configurado corretamente)

### Solução com valores absurdos (10^9, 10^10)

**Causa**: Bug na conversão CSR→COO ou uso incorreto da versão simétrica  
**Solução**: Verifique que está usando a última versão do código

## Performance

Para o problema de teste (9 equações):

| Solver   | Tempo | Memória (NNZ) |
|----------|-------|---------------|
| Skyline  | ~1ms  | ~30           |
| Pardiso  | ~1ms  | 29            |
| MUMPS    | ~2ms  | 45            |

**NNZ** = Número de elementos não-zeros armazenados

## Documentação Adicional

- **Resumo**: [MUMPS_INTEGRATION_SUMMARY.md](MUMPS_INTEGRATION_SUMMARY.md)
- **Bug Report**: [BUG_REPORT_MUMPS_SYMMETRIC.md](BUG_REPORT_MUMPS_SYMMETRIC.md)
- **Manuais**: Ver pasta `manuals/`
  - `MUMPS_5.8.1.pdf` - Manual completo do MUMPS
  - `Pardiso.pdf` - Manual do Pardiso

## Contato

Se encontrar problemas ou tiver dúvidas, consulte os arquivos de documentação acima que contêm investigação detalhada e análise do problema.

---

**Última atualização**: 17 de Dezembro de 2025
