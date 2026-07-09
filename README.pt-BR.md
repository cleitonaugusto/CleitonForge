# CleitonForge

[English](README.md) | **[Português]**

[![GitHub Sponsors](https://img.shields.io/github/sponsors/cleitonaugusto?label=Apoiar&logo=GitHub&color=ea4aaa)](https://github.com/sponsors/cleitonaugusto)
[![Crates.io](https://img.shields.io/crates/v/cleitonforge)](https://crates.io/crates/cleitonforge)
[![PyPI](https://img.shields.io/pypi/v/cleitonforge)](https://pypi.org/project/cleitonforge/)

**CleitonForge** é uma camada neutra e open-source de benchmarking e
interoperabilidade para simuladores quânticos, escrita em Rust. O binário
CLI se chama `cforge`.

---

## O problema

O ecossistema de software quântico está fragmentado entre frameworks
incompatíveis — Qiskit (IBM), Cirq (Google), Amazon Braket, PennyLane
(Xanadu) — e vários simuladores isolados escritos em Rust (QuantRS2, qoqo,
q1tsim, qvnt). Cada framework mede desempenho de forma diferente, usa
formatos de circuito distintos e define os mesmos algoritmos de maneiras
incompatíveis.

A literatura acadêmica confirma o problema:

- Um benchmark chamado "grover" pode implementar o algoritmo de Grover com
  um oráculo completamente diferente em cada ferramenta, tornando os
  resultados incomparáveis entre frameworks.
- Ferramentas quânticas são difíceis de integrar em pipelines de CI/CD pela
  ausência de formatos de intercâmbio padronizados.
- Esforços de benchmarking agnóstico existem (QED-C, MQT Bench), mas os
  próprios autores reconhecem que são "um passo", não uma solução definitiva.

## Por que Rust, por que agora

- **OpenQASM 3** é o padrão de intercâmbio de fato da indústria, mas sua
  implementação de referência e a maior parte do tooling são em Python.
- Parsers Rust para OpenQASM (`oq3_semantics`, `oq3_parser`, `oq3_syntax`)
  existem com 680k+ downloads cada, mas estão sem manutenção há mais de um ano.
- Simuladores em Rust existem mas são silos isolados sem API comum.
- Rust entrega desempenho previsível e de baixo overhead — crítico para
  benchmarking onde ruído de medição deve ser minimizado.

## Posicionamento

O CleitonForge não compete com Qiskit, Cirq, IBM ou qualquer fornecedor de
hardware. O objetivo é ser a **camada neutra** que todos podem usar — um
árbitro imparcial que mede e compara sem favorecer nenhum backend.

- Nunca otimizado para favorecer um simulador específico
- Rastreia a proveniência de cada gate (qual framework, qual nome original)
- Arquitetura de plugin (trait `SimulationBackend`) para que qualquer backend
  seja adicionado pela comunidade sem tocar no núcleo

---

## Arquitetura

```
CleitonForge/                    workspace Rust
├── cforge-core/                 IR canônico: Circuit, GateKind, Operation
├── cforge-parser/               Parsers OpenQASM 2 + OpenQASM 3
├── cforge-backends/             Trait SimulationBackend + implementações
├── cforge-metrics/              Fidelidade, profundidade, tempo, memória
├── cforge-cli/                  Binário `cforge` (clap + comfy-table)
│   └── examples/
│       └── compare_grover.rs   Algoritmo de Grover — exemplo via API Rust
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

**Planejado:** backends adicionais (qoqo, q1tsim), modelagem de ruído,
bindings Python (PyO3), cobertura estendida do OpenQASM 3, dashboard web.

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
