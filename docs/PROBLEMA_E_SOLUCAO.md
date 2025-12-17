# Resumo do Problema e Solução

## Problema Original
Erro `MUMPS was called with an invalid value for JOB. JOB = 3` (INFO(1) = -3)

## Causa Raiz
O erro JOB=3 inválido era um **sintoma secundário**. O problema real era:

1. **Matriz numericamente singular** (INFO(1) = -10): A fatorização do MUMPS falhava porque a matriz tinha um pivô zero/muito pequeno.

2. **Copy constructors inadequados**: A classe `TPZMumpsSolver` tinha copy/move constructors definidos como `= default`, o que causava:
   - Cópia superficial (shallow copy) da estrutura `DMUMPS_STRUC_C fMumpsData`
   - O MUMPS mantém estado interno que não pode ser copiado bit-a-bit
   - Quando a fatorização falhava silenciosamente, o estado ficava inconsistente
   - Na tentativa de Solve(), o MUMPS via JOB=3 sem as fases anteriores (INIT, ANALYSIS, FACTORIZE) terem sido completadas com sucesso

## Correções Implementadas

### 1. TPZMumpsSolver.h
- ❌ Deletei copy constructor e copy assignment (não podem ser copiados)
- ✅ Implementei move constructor e move assignment adequados
- ✅ Modificado o método Clone() para não usar copy constructor

### 2. TPZMumpsSolver.cpp  
- ✅ Implementado move constructor que transfere ownership do estado MUMPS
- ✅ Implementado move assignment que libera o estado antigo antes de mover
- ✅ Adicionada verificação de erro após a fase de FACTORIZE
- ✅ Clone() agora cria um novo solver não-inicializado com os mesmos parâmetros

### 3. TPZSYSMPMumps.h e TPZYSMPMumps.h
- ✅ Declarados copy/move constructors e operators explicitamente
- ✅ Copy constructor NÃO copia o fMumpsControl (cada matriz precisa de sua própria instância MUMPS)

### 4. TPZSYSMPMumps.cpp e TPZYSMPMumps.cpp
- ✅ Implementados copy constructors que criam novo fMumpsControl vazio
- ✅ Implementados move constructors que transferem o fMumpsControl
- ✅ Copy assignment reseta fDecomposed para ENoDecompose

## Problema Subjacente (Ainda Presente)
A **matriz é numericamente singular** (INFO(1) = -10). Isso indica:

- Condições de contorno podem estar mal especificadas
- Para Darcy flow com 2 Neumann BCs (top/bottom) e 2 Dirichlet BCs (left/right), pode faltar uma restrição adicional
- A permeabilidade ou outras propriedades do material podem estar causando a singularidade

## Próximos Passos Sugeridos
1. Revisar as condições de contorno do problema de Darcy
2. Considerar usar apenas Dirichlet ou adicionar uma restrição de valor médio
3. Verificar se a formulação do problema está correta para as BCs escolhidas

## Fases Corretas do MUMPS (para referência)
1. **JOB = -1 (INIT)**: Inicialização
2. **JOB = 1 (ANALYSIS)**: Análise simbólica
3. **JOB = 2 (FACTORIZE)**: Fatorização numérica  
4. **JOB = 3 (SOLVE)**: Resolução do sistema
5. **JOB = -2 (END)**: Finalização e liberação de memória

Você estava seguindo as fases corretamente, mas o erro na fatorização (passo 3) impedia o MUMPS de aceitar JOB=3.
