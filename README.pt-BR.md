# CleitonForge

[English](README.md) | **[Português]**

[![GitHub Sponsors](https://img.shields.io/github/sponsors/cleitonaugusto?label=Apoiar&logo=GitHub&color=ea4aaa)](https://github.com/sponsors/cleitonaugusto)
[![Crates.io](https://img.shields.io/crates/v/cleitonforge)](https://crates.io/crates/cleitonforge)
[![PyPI](https://img.shields.io/pypi/v/cleitonforge)](https://pypi.org/project/cleitonforge/)

**CleitonForge** é um fuzzer diferencial para compiladores e simuladores
quânticos, escrito em Rust. Ele gera circuitos aleatórios, roda cada um em mais
de uma implementação e reporta os casos em que elas discordam. Quando acha um,
encolhe o circuito até a menor versão que ainda falha. O binário CLI se chama
`cforge`.

Ele achou um bug de corretude no transpiler do Qiskit:
[issue #16594](https://github.com/Qiskit/qiskit/issues/16594). O passe
`CommutativeCancellation` cancelava `sxdg sxdg sx`, que é um X, até sobrar um
circuito vazio, do nível 2 de otimização para cima. Sem erro, sem aviso, e o
resultado medido muda. Um mantenedor do core confirmou e a correção saiu no
Qiskit 2.5.1.

---

## O problema

Simulador e compilador são difíceis de testar porque em geral não existe com o
que comparar a resposta. Uma amplitude errada continua sendo uma amplitude bem
formada. Nada levanta exceção, nada vai pro log, e o número que volta é
idêntico a um número certo. Você só percebe quando já sabe qual devia ser a
resposta, que é justamente o caso em que não precisava da ferramenta.

Teste diferencial contorna isso. Em vez de perguntar "isso está certo", ele
pergunta "essas duas concordam". Duas implementações da mesma especificação
deveriam devolver o mesmo estado para o mesmo circuito. Quando não devolvem,
uma das duas está errada, e a discordância é o sinal.

A ideia é essa. O resto é deixar isso preciso o bastante para valer a leitura:

- **O gerador** mira nos limiares numéricos em que um compilador ramifica:
  múltiplos de pi/2, ângulos logo abaixo de um corte de síntese, sequências que
  colapsam para a identidade. Ângulo aleatório uniforme quase nunca cai ali.
- **O oráculo** compara a menos de fase global, porque fase global não é
  observável e acusá-la enterraria os achados reais no ruído. As harnesses dos
  compiladores externos podem passar o veredicto para o
  [MQT QCEC](https://github.com/cda-tum/mqt-qcec), que raciocina simbolicamente
  em vez de simular.
- **O shrinker** reduz a falha a uma testemunha 1-minimal, ou seja, tirar
  qualquer gate faz o bug sumir. A #16594 saiu dele com três gates num qubit só.
- **A triagem** confere todo veredicto contra o operador exato antes de chamar
  qualquer coisa de bug, e separa classe conhecida de achado novo por rotação
  acumulada, não por nome de gate.

## A taxa de acerto, sem enfeite

Mais ou menos um bug real a cada cem mil circuitos. Essas stacks são maduras e
bem testadas, e a maioria das discordâncias acaba sendo uma transformação
declarada sobre a qual o oráculo não tinha sido avisado, não um defeito. As
campanhas contra tket, PennyLane e Cirq voltaram limpas.

Então o valor aqui não é volume. É que o que sobrevive à triagem é real,
mínimo, e reproduzível por alguém que não confia em você. Está em
[`bug-zoo/`](bug-zoo), um JSON por achado, com a lista de gates, um reprodutor
em OpenQASM 2 e o registro de como foi encontrado. As divergências entre
backends trazem também a semente do gerador e a distância em cada nível de
oráculo.

## Por que Rust

- O gerador, o vetor de estado e o shrinker estão todos no caminho quente de
  uma campanha que roda dezenas de milhares de circuitos.
- O transpiler do Qiskit está migrando para Rust, e o bug acima está em código
  Rust. Ajuda ler a linguagem daquilo que você está fuzzando.
- No lado do benchmarking, overhead baixo e previsível importa, porque ruído de
  medição é justamente o que está sendo medido.

## O lado do benchmarking

O CleitonForge começou como camada de benchmarking entre backends, e essa parte
continua funcionando e está documentada abaixo. O `cforge run` executa o mesmo
arquivo OpenQASM em mais de um backend e reporta tempo, memória, profundidade e
fidelidade lado a lado. Hoje é um componente do fuzzer, não o ponto do projeto.

Duas regras que ficaram daquele período, porque são o que faz a comparação
significar alguma coisa:

- Nunca otimizado para favorecer um simulador específico.
- Arquitetura de plugin (trait `SimulationBackend`), para que um backend seja
  adicionado sem tocar no IR do núcleo.

---

## Arquitetura

```
CleitonForge/                    workspace Rust
├── cforge-core/                 IR canônico: Circuit, GateKind, Operation
├── cforge-parser/               Parsers OpenQASM 2 + OpenQASM 3
├── cforge-backends/             Trait SimulationBackend + implementações
├── cforge-metrics/              Fidelidade, profundidade, tempo, memória
├── cforge-fuzz/                 O fuzzer
│   ├── generator.rs             Gerador ponderado de circuitos
│   ├── oracle.rs                Oráculos de divergência (N1 amplitude, N2 prob.)
│   ├── shrinker.rs              Redução gulosa a uma testemunha 1-minimal
│   ├── triage.rs                Classe conhecida vs achado novo
│   └── zoo.rs                   Saída JSON para o bug-zoo/
├── cforge-cli/                  Binário `cforge` (clap + comfy-table)
│   └── examples/
│       └── compare_grover.rs   Algoritmo de Grover, exemplo via API Rust
├── cforge-py/                   Bindings Python (PyO3)
├── tools/                       Harnesses para compiladores externos
│   ├── fuzz_qiskit*.py          Qiskit: pipeline completo, passes isolados, wide
│   ├── fuzz_pytket*.py          tket
│   ├── fuzz_pennylane*.py       Devices e transforms do PennyLane
│   ├── fuzz_cirq.py             Cirq
│   ├── oracle_qcec.py           MQT QCEC como veredicto de três vias
│   └── triage_known.py          Triagem por rotação acumulada
├── bug-zoo/                     Contraexemplos minimizados, um JSON cada
└── examples/
    └── bell.qasm               Estado de Bell em OpenQASM 2
```

### Fluxo de dados

```
arquivo .qasm
      │
      ▼
cforge-parser  ──►  Circuit (IR canônico)
                         │
              ┌──────────┴──────────┐
              ▼                     ▼
  NativeStateVector          QuantRS2Backend
     Backend                       │
              └──────────┬──────────┘
                         ▼
              cforge-metrics: fidelidade, tempo, memória
                         │
                         ▼
              cforge-cli: tabela / saída JSON
```

---

## Início Rápido

### Pré-requisitos

- Rust 1.96+ (`rustup update stable`)
- Nenhum simulador externo necessário — todas as dependências são crates Rust puros

### Compilar

```bash
git clone https://github.com/cleitonaugusto/CleitonForge.git
cd CleitonForge
cargo build --release
# binário em: target/release/cforge
```

### Rodar um circuito nos dois backends

```bash
cforge run --circuit examples/bell.qasm --backends statevector,quantrs2 --shots 1024
```

```
Circuit: 2 qubits  |  2 gates  |  depth 2  |  seed 0xdeadbeefcafebabe
┌──────────────────────┬───────────┬────────┬───────┬───────┬──────────┬───────┐
│ Backend              │ Time (ms) │ Memory │ Depth │ Gates │ Fidelity │ Shots │
╞══════════════════════╪═══════════╪════════╪═══════╪═══════╪══════════╪═══════╡
│ statevector-native   │ 0.002     │ 64 B   │ 2     │ 2     │ 1.000000 │ 1024  │
├──────────────────────┼───────────┼────────┼───────┼───────┼──────────┼───────┤
│ statevector-quantrs2 │ 0.004     │ 64 B   │ 2     │ 2     │ 1.000000 │ 1024  │
└──────────────────────┴───────────┴────────┴───────┴───────┴──────────┴───────┘
```

### Validar um circuito sem simulação

```bash
cforge validate --circuit examples/bell.qasm
```

```
File   : examples/bell.qasm
Qubits : 2
Gates  : 2
Depth  : 2
By gate:
  cx       1
  h        1
Status : OK
```

### Exportar resultados como JSON (para scripts / CI)

```bash
cforge run --circuit examples/bell.qasm --shots 1024 --format json
```

### Exemplo end-to-end: Algoritmo de Grover

Busca de 3 qubits pelo estado |101⟩, 2 iterações, ambos os backends:

```bash
cargo run --example compare_grover -p cforge-cli
```

```
Target state : |101⟩  (índice 5)
Circuit      : 43 gates  |  depth 21

Backend  : statevector-native
  Top state  : |101⟩  prob = 0.9453 (94.5 %)   [teórico: sin²(5θ) ≈ 94.8 %]
  Fidelity   : 1.00000000

Cross-backend fidelity (native vs quantrs2): 1.00000000
Both backends agree: YES ✓
```

O valor medido de 94.5% confirma a previsão teórica para N=8, M=1, k=2
iterações — validando a corretude de ambos os backends de forma independente.

---

## Referência da CLI

```
cforge run
  --circuit <path>          Arquivo OpenQASM 2 ou 3 (detecção automática)
  --backends <lista>        Vírgula-separado: statevector, quantrs2  [padrão: ambos]
  --shots <n>               Medições; 0 = somente statevector  [padrão: 0]
  --seed <u64>              Semente do PRNG para contagens reproduzíveis
  --format <table|json>     Formato de saída  [padrão: table]

cforge validate
  --circuit <path>          Analisa e exibe estatísticas; sai com código 1 se inválido
```

---

## Conjunto de Gates Suportados

CleitonForge implementa o conjunto completo do **`stdgates.inc` do OpenQASM 3**:

| Categoria | Gates |
|---|---|
| 1 qubit, sem parâmetros | `id` `x` `y` `z` `h` `s` `sdg` `t` `tdg` `sx` `sxdg` |
| 1 qubit, paramétricos | `rx(θ)` `ry(θ)` `rz(θ)` `p(θ)` `u(θ,φ,λ)` |
| 2 qubits | `cx` `cy` `cz` `ch` `csx` `crx` `cry` `crz` `cp` `cu` `swap` |
| 3 qubits | `ccx` (Toffoli) `cswap` (Fredkin) |
| Aliases | `cnot`→`cx` `u1`→`p` `u2(φ,λ)`→`u(π/2,φ,λ)` `u3`→`u` `ccnot`→`ccx` `fredkin`→`cswap` |

---

## Formatos de Entrada

| Formato | Detectado por | Observações |
|---|---|---|
| OpenQASM 2.0 | tudo que não começa com `OPENQASM 3` | `include` resolvido pelo diretório do arquivo |
| OpenQASM 3 | cabeçalho `OPENQASM 3` | Gate calls via `oq3_semantics` |

Ambos os parsers suportam **aplicação de gate em registro inteiro**:
`h q;` em um registro de 3 qubits expande para `h q[0]; h q[1]; h q[2];`

---

## API Rust

```rust
use cforge_core::{Circuit, GateKind, Operation};
use cforge_backends::{DEFAULT_SEED, NativeStateVectorBackend, SimulationBackend};
use cforge_metrics::{compute_stats, measure};

// Construir um circuito Bell
let mut circuit = Circuit::new(2);
circuit.push(Operation::new(GateKind::H,  vec![0], vec![]));
circuit.push(Operation::new(GateKind::Cx, vec![0, 1], vec![]));

// Executar e medir
let stats = compute_stats(&circuit);
let result = NativeStateVectorBackend.run(&circuit, 1024, DEFAULT_SEED)?;

println!("profundidade = {}", stats.depth);               // 2
println!("prob |00⟩ = {:.3}", result.statevector[0].norm_sqr()); // 0.500
```

---

## Adicionando um Novo Backend

Implemente o trait `SimulationBackend` (um único método):

```rust
use cforge_backends::{BackendError, SimulationBackend, SimulationResult};
use cforge_core::Circuit;

pub struct MeuBackend;

impl SimulationBackend for MeuBackend {
    fn name(&self) -> &str { "meu-backend" }

    fn run(
        &self,
        circuit: &Circuit,
        shots: usize,
        seed: u64,
    ) -> Result<SimulationResult, BackendError> {
        // ... lógica de simulação
    }
}
```

---

## Medição de Memória

Para simulações statevector o pico teórico de memória é:

```
2^n_qubits × 16 bytes  (dois f64 por amplitude Complex128)
```

No Linux, o CleitonForge mede o delta real de RSS via `/proc/self/status`
enquanto o statevector está vivo na memória.

| Qubits | Pico teórico |
|--------|-------------|
| 10     | 16 KB       |
| 20     | 16 MB       |
| 22     | 64 MB (máx) |

---

## Status do Projeto

| Fase | Crate             | Status |
|------|-------------------|--------|
| 0    | workspace setup   | ✓      |
| 1    | `cforge-core`     | ✓      |
| 2    | `cforge-parser`   | ✓      |
| 3    | `cforge-backends` | ✓      |
| 4    | `cforge-metrics`  | ✓      |
| 5    | `cforge-cli`      | ✓      |
| 6    | exemplos + docs   | ✓      |
| 7    | `cforge-py`       | ✓      |
| 8    | `cforge-fuzz`     | ✓      |

Os bindings Python já estão publicados: `pip install cleitonforge`.

**Planejado:** backends adicionais (qoqo, q1tsim), modelagem de ruído,
cobertura estendida do OpenQASM 3.

---

## Contribuindo

Issues e pull requests são bem-vindos. A arquitetura é intencionalmente
modular — adicionar um backend, uma nova métrica ou um novo formato de
entrada não requer tocar no IR do núcleo.

```bash
cargo test --workspace
cargo clippy --workspace
```

## Patrocínio

CleitonForge é desenvolvido e mantido por um pesquisador independente.
Se esta ferramenta economizou seu tempo ou ajudou a encontrar um bug, considere apoiar:

**➜ [github.com/sponsors/cleitonaugusto](https://github.com/sponsors/cleitonaugusto)**

---

## Licença

Licenciado sob a Apache License, Version 2.0. Veja [LICENSE](LICENSE).
