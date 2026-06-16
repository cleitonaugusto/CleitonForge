# CleitonForge — Brief completo para implementação

## Contexto do projeto

Este é um projeto open source em Rust chamado **CleitonForge**: uma camada
de benchmark e interoperabilidade entre simuladores de computação quântica.
O nome do crate principal no código e no crates.io será `cleitonforge`
(em minúsculas, seguindo a convenção do ecossistema Rust). O binário do
CLI será `cforge` para uso rápido no terminal.

### O problema que resolve (validado por pesquisa, não suposição)

O ecossistema de computação quântica está fragmentado entre frameworks
incompatíveis: Qiskit (IBM), Cirq (Google), Amazon Braket, PennyLane (Xanadu),
e múltiplos simuladores em Rust (QuantRS2, qoqo, q1tsim). Cada um tem sua
própria forma de medir performance e corretude.

A literatura acadêmica confirma o gap:
- Não há consolidação sobre quais benchmarks usar para avaliação empírica
  de ferramentas quânticas — um benchmark chamado "grover" pode implementar
  o algoritmo de Grover com qualquer oráculo diferente, tornando resultados
  incomparáveis entre ferramentas.
- O tooling quântico continua difícil de integrar a IDEs e pipelines CI/CD
  por falta de formatos padronizados e interoperabilidade real.
- Tentativas de arquitetura agnóstica de plataforma existem (QED-C, MQT
  Bench) mas os próprios autores admitem que não resolvem completamente
  o problema — são "um passo", não uma solução definitiva.

### Por que Rust e por que agora

- OpenQASM 3 já é o padrão de intercâmbio de fato na indústria, mas sua
  implementação de referência é em Python.
- Já existem parsers OpenQASM em Rust (crates `oq3_parser`, `oq3_lexer`,
  `oq3_semantics`, `oq3_syntax`, `oq3_semantics`) com mais de 680 mil
  downloads cada, mas sem atualização há mais de um ano — sinal de
  oportunidade de assumir/estender, não de começar do zero.
- Simuladores Rust existem mas são fragmentados e isolados: QuantRS2,
  qoqo, q1tsim, qvnt — nenhum tem camada de benchmark unificada.
- Rust dá performance superior a implementações Python para análise
  estrutural pesada de circuitos grandes.

### Estratégia de posicionamento (importante para todas as decisões técnicas)

NÃO competir com Qiskit, Cirq, IBM ou qualquer fornecedor de hardware.
A estratégia é ser a CAMADA NEUTRA que todos eles podem usar — o "juiz"
imparcial que mede e compara, sem favorecer nenhum backend. Isso significa:

- Nunca otimizar a ferramenta para favorecer um simulador específico
- Sempre dar crédito e citar a proveniência de cada framework de origem
- Arquitetura de plugin/trait para que qualquer backend novo seja
  adicionado pela comunidade sem reescrever o core

## Objetivo do dono do projeto

Analista de sistemas brasileiro, forte em Rust, entusiasta e estudioso de
computação quântica (já fez cursos, conhece e simula os algoritmos
existentes). Busca reconhecimento técnico real e contribuição de impacto
duradouro — não retorno financeiro rápido. O modelo mental é construir
algo no espírito do que ferramentas como cargo, cargo-audit ou MQT Bench
representam para seus respectivos ecossistemas: infraestrutura citada e
usada pela comunidade, não produto comercial fechado.

Licença: Apache 2.0 (mesma do OpenQASM e do ecossistema Rust em geral).
Open source desde o primeiro commit.

## Arquitetura técnica completa

```
cleitonforge/                          (workspace Rust)
├── Cargo.toml                   (workspace root)
├── cforge-core/                 (IR canônico + tipos compartilhados)
│   └── src/
│       ├── ir.rs                (grafo de circuito com proveniência)
│       ├── gate.rs               (definições de gates universais)
│       └── metrics.rs            (tipos de métricas padronizadas)
├── cforge-parser/                (camada de entrada multi-formato)
│   └── src/
│       ├── qasm2.rs               (wrapper sobre crate `qasm` ou `openqasm`)
│       ├── qasm3.rs               (wrapper sobre `oq3_parser`/`oq3_semantics`)
│       └── translate.rs          (tradução para IR canônico)
├── cforge-backends/               (trait + implementações plugáveis)
│   └── src/
│       ├── trait_def.rs           (trait `SimulationBackend`)
│       ├── statevector.rs        (backend usando QuantRS2 ou implementação própria)
│       └── tensor.rs              (backend tensor network, fase posterior)
├── cforge-metrics/                 (cálculo de métricas padronizadas)
│   └── src/
│       ├── fidelity.rs
│       ├── circuit_stats.rs       (profundidade, contagem de gates)
│       └── performance.rs         (tempo, memória)
├── cforge-cli/                     (interface de linha de comando)
│   └── src/main.rs
└── examples/
    └── compare_grover.rs           (exemplo: comparar Grover entre 2 backends)
```

## Plano de execução passo a passo

### Fase 0 — Setup do workspace (primeiro passo, fazer agora)

1. Criar workspace Cargo com os crates listados acima como membros
2. Cada crate começa vazio com `lib.rs` básico e `Cargo.toml` próprio
3. Configurar `rust-toolchain.toml` fixando uma versão estável recente
4. Criar `README.md` inicial explicando o projeto (usar o texto da seção
   "Contexto do projeto" acima como base, traduzido para inglês — projetos
   Rust open source de alcance internacional devem ter README em inglês)
5. Adicionar `LICENSE` Apache 2.0
6. Configurar `.gitignore` padrão Rust
7. Fazer commit inicial e criar repositório no GitHub

### Fase 1 — cforge-core (fundação de tipos)

Objetivo: definir o IR canônico antes de qualquer parser ou backend.

1. Definir `struct Circuit` representando um circuito quântico como grafo:
   - lista de qubits com índices
   - lista de operações (gates) em sequência, cada uma com:
     - tipo de gate (enum `GateKind`: H, X, Y, Z, CNOT, RZ, RX, RY, etc — começar
       com o conjunto universal mínimo de OpenQASM 3 `stdgates.inc`)
     - qubits afetados
     - parâmetros (para gates parametrizados como rotações)
     - metadado de proveniência: `source_framework: Option<String>` e
       `original_gate_name: Option<String>` para rastrear de onde veio
2. Definir `enum GateKind` cobrindo o gate set padrão de stdgates.inc do
   OpenQASM 3 (consultar especificação em openqasm.com se necessário)
3. Definir `struct MetricsResult` com campos: fidelity, depth, gate_count,
   execution_time_ms, memory_bytes, backend_name
4. Escrever testes unitários para construção manual de circuitos simples
   (ex: Bell state com H + CNOT) e verificar que a estrutura é coerente
5. Não se preocupar com parsing ainda — isso é só a estrutura de dados

### Fase 2 — cforge-parser (entrada OpenQASM)

Objetivo: ler arquivos .qasm e produzir um `Circuit` do cforge-core.

1. Adicionar como dependência o crate `oq3_parser` e `oq3_semantics` do
   crates.io (ou fazer fork local se estiverem desatualizados demais —
   verificar primeiro se compilam com a toolchain Rust atual)
2. Para OpenQASM 2, usar o crate `qasm` (mais simples, mais maduro)
3. Implementar função `fn parse_qasm3(source: &str) -> Result<Circuit, ParseError>`
   que usa o parser externo e traduz o AST resultante para o IR do cforge-core
4. Implementar a mesma função para QASM 2
5. Testar com os exemplos oficiais do repositório github.com/openqasm/openqasm
   (pasta `examples/`) — usar pelo menos 5 circuitos de exemplo diferentes
6. Tratar erros de parsing de forma clara, com mensagem indicando linha e
   coluna do problema (os parsers OpenQASM em Rust já fornecem isso)

### Fase 3 — cforge-backends (primeiro backend de execução)

Objetivo: executar um `Circuit` e produzir resultado de simulação.

1. Definir trait:
   ```rust
   trait SimulationBackend {
       fn name(&self) -> &str;
       fn run(&self, circuit: &Circuit, shots: usize) -> Result<SimulationResult, BackendError>;
   }
   ```
2. Implementar o PRIMEIRO backend simples: um simulador state-vector básico
   escrito do zero (não depender de QuantRS2 ainda na primeira versão,
   para ter controle total e entender o problema profundamente) —
   usar `num-complex` para números complexos e implementar aplicação de
   gates como multiplicação de matriz no vetor de estado
   - Começar suportando até 15-20 qubits (suficiente para validação,
     2^20 entradas é tratável em memória comum)
3. Depois de funcionar, adicionar um SEGUNDO backend que efetivamente
   chama QuantRS2 como dependência externa, para provar que a arquitetura
   de trait realmente permite múltiplos backends
4. Escrever teste de integração que roda o MESMO circuito (Bell state) nos
   dois backends e verifica que o resultado de fidelidade é equivalente

### Fase 4 — cforge-metrics (padronização de medição)

Objetivo: calcular métricas comparáveis entre backends diferentes.

1. `fidelity.rs`: calcular fidelidade entre o estado produzido e um estado
   de referência esperado (para circuitos onde o resultado é conhecido)
2. `circuit_stats.rs`: calcular profundidade do circuito (número de
   "camadas" de gates paralelizáveis) e contagem total de gates por tipo
3. `performance.rs`: wrapper que mede tempo de execução wall-clock e,
   se possível, uso de memória pico durante a chamada ao backend
4. Importante: a métrica de profundidade e contagem de gates deve ser
   calculada sobre o IR CANÔNICO, não sobre a representação nativa de
   cada backend — isso é o que garante comparação justa

### Fase 5 — cforge-cli (interface de uso real)

Objetivo: comando de terminal usável por qualquer pesquisador.

1. Usar `clap` para definir os argumentos de linha de comando
2. Comando básico: `cforge run --circuit bell.qasm --backends statevector,quantrs2`
3. Saída no terminal: tabela comparativa simples (pode usar crate
   `comfy-table` ou similar) mostrando lado a lado as métricas de cada
   backend para o mesmo circuito
4. Comando adicional: `cforge validate --circuit arquivo.qasm` que só
   faz parsing e mostra estatísticas do circuito sem executar simulação

### Fase 6 — Exemplo end-to-end e documentação

1. Criar `examples/compare_grover.rs`: implementação do algoritmo de
   Grover (3 qubits, oráculo simples) escrita diretamente em Rust usando
   a API do cforge-core, sem depender de arquivo .qasm externo
2. Rodar esse exemplo nos dois backends implementados e documentar o
   resultado no README como prova de conceito
3. Escrever seção "Quick Start" no README com os comandos exatos para
   reproduzir o exemplo

## Decisões técnicas já tomadas (não reabrir discussão sobre isso)

- Linguagem: Rust, edição 2021 ou mais recente
- Workspace multi-crate, não um único crate monolítico
- Licença: Apache 2.0
- Formato de entrada prioritário: OpenQASM 2 e 3 (não QIR nem Quil na v1)
- Primeiro backend: state-vector simples escrito do zero, depois QuantRS2
- Sem GPU/SIMD avançado na v1 — isso é otimização para depois que a
  arquitetura básica estiver provada
- Sem interface web na v1 — CLI primeiro

## O que fazer se travar em alguma decisão não coberta aqui

Se surgir uma decisão de arquitetura não especificada neste documento,
a regra geral é: escolher a opção mais simples que funciona, documentar
a decisão com um comentário no código explicando o porquê, e seguir.
Não pausar o progresso esperando validação externa para decisões de
implementação de baixo nível (nomes de função, organização interna de
um módulo, escolha entre duas crates equivalentes para uma tarefa
pequena).

## Primeira ação concreta

Comece pela Fase 0. Crie o workspace, os crates vazios, o README inicial
em inglês, a licença, e faça o primeiro commit. Depois siga para a Fase 1
sem esperar confirmação, a menos que encontre um bloqueio técnico real
(por exemplo, um crate que não compila).
