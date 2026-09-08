#!/usr/bin/env python3
"""Checks a book's Rust code against itself and against this workspace.

Used by both books. Each passes its own directory and, optionally, a
`check-code.json` beside its markdown holding `chapter_builds`, `import_layout`
and `elided`.

Two checks, both of which caught real defects when this was first run:

  1. Undefined symbols — a function or type the text asks the reader to call
     and never shows. This found `require_param` (used ten times), `ParseError`
     (nine), the four parser helpers, and `fidelity_sv`, which turned out to be
     a wrong name for `state_fidelity`.

  2. Chapter assembly — for chapters that build one artifact, gather the blocks
     that make it and compile them together. This found a missing `?` on the
     `rz` arm of the gate table, which every sibling arm had.

Deliberately not attempted: compiling every block in isolation. Blocks continue
each other and lean on earlier chapters, so isolated compilation reports missing
context as failure and buries the real defects. An early version of this script
did that and reported 61 of 91 blocks "broken", nearly all of them fine.

    python3 tools/check-book-code.py ../quantum-rust-book
    python3 tools/check-book-code.py ../quantum-rust-book-B/manuscript

Exits non-zero if a check fails.
"""
import json
import pathlib
import re
import subprocess
import sys
import tempfile

if len(sys.argv) < 2:
    sys.exit(f"usage: {sys.argv[0]} <book-directory>")
HERE = pathlib.Path(sys.argv[1]).resolve()
REPO = pathlib.Path(__file__).resolve().parent.parent

# Per-book settings, if the book has any. A book with no config still gets the
# undefined-symbol check, which is the one that needs no setup.
_cfg_path = HERE / "check-code.json"
_cfg = json.loads(_cfg_path.read_text()) if _cfg_path.is_file() else {}

# Names that resolve without the book defining them: std, common crates, and
# attribute or macro spellings that look like calls to a regex.
KNOWN = {
    "main", "new", "from", "into", "len", "push", "pop", "iter", "collect",
    "map", "unwrap", "expect", "format", "println", "print", "eprintln", "vec",
    "some", "none", "clone", "to_vec", "to_string", "min", "max", "abs",
    "sqrt", "powi", "powf", "cos", "sin", "tan", "exp", "ln", "sum", "filter",
    "enumerate", "trim", "split", "parse", "matches", "assert", "write",
    "default", "zip", "rev", "chars", "bytes", "get", "insert", "remove",
    "contains", "starts_with", "ends_with", "replace", "sort", "sort_by",
    "fold", "any", "all", "count", "take", "skip", "join", "repeat", "norm",
    "norm_sqr", "conj", "re", "im", "sqrt_f64", "todo", "unimplemented",
    "panic", "dbg", "include_str", "env", "cfg", "derive", "test", "bench",
    "pyclass", "pymethods", "pymodule", "pyfunction", "pyo3", "error", "arg",
    "command", "serde", "allow", "deny", "warn", "doc", "inline", "repr",
    "should_panic", "ignore", "wrap_pyfunction", "add_class", "add_function",
    "gen", "gen_range", "thread_rng", "seed_from_u64", "abs_diff_eq",
    "assert_abs_diff_eq", "assert_eq", "assert_ne", "matches_macro",
}

# Rust keywords that a naive "identifier followed by (" also matches.
KEYWORDS = {
    "if", "for", "while", "match", "return", "fn", "let", "move", "loop",
    "else", "impl", "where", "as", "in", "mut", "ref", "use", "mod", "pub",
    "self", "super", "crate", "type", "const", "static", "struct", "enum",
    "trait", "unsafe", "async", "await", "dyn", "box", "and", "or", "not",
}

# Types and constructors the prelude and the common crates provide.
PRELUDE = {
    "Ok", "Err", "Some", "None", "Box", "Vec", "String", "Self", "Option",
    "Result", "HashMap", "HashSet", "BTreeMap", "Rc", "Arc", "Mutex", "RwLock",
    "Cow", "Path", "PathBuf", "Instant", "Duration", "Complex", "Complex64",
    "Sync", "Send", "Fn", "FnMut", "FnOnce", "Iterator", "Default", "Debug",
    "Clone", "Copy", "PartialEq", "Eq", "Hash", "Ord", "PartialOrd", "Display",
    "Error", "From", "Into", "TryFrom", "Formatter", "Write", "Read", "Bearer",
}

# Named on purpose without being shown. Each one is either elided with a note
# in the text, quoted from someone else's source, or console output that a
# regex cannot tell from code. Anything not on this list and not defined is a
# gap the reader will fall into.
ELIDED = {
    # "apply_to_q1tsim is a mechanical match on GateKind" — ch. 8, stated.
    "apply_to_q1tsim",
    # "not written out here" — ch. 14, stated.
    "apply_gate_to_density_matrix",
    "apply_single_qubit_channel_to_full_state",
    # Illustrative pseudocode of the assertion antipattern, not runnable code.
    "magnitudes",
    # Quoted from Qiskit's own source in the bug-hunting chapter.
    "StandardGate",
    # Console output inside a fenced block: "Job submitted: abc123def456".
    "Job",
}

ELIDED |= set(_cfg.get("elided", []))

CALL = re.compile(r"(?<![.\w:])([a-z_][a-z0-9_]{2,})\s*\(")
USE = re.compile(r"\buse\s+([^;]+);")
GLOB = re.compile(r"\buse\s+[^;]*::\*\s*;")
VARIANTS = re.compile(r"\benum\s+[A-Z][A-Za-z0-9_]*\s*\{(.*?)\n\}", re.S)
DEFN = re.compile(r"\bfn\s+([a-z_][a-z0-9_]*)")
CLOSURE = re.compile(r"\blet\s+(?:mut\s+)?([a-z_][a-z0-9_]*)\s*=\s*(?:move\s*)?\|")
TYPEDEF = re.compile(r"\b(?:struct|enum|trait|union)\s+([A-Z][A-Za-z0-9_]*)"
                     r"|\btype\s+([A-Z][A-Za-z0-9_]*)")
TYPEUSE = re.compile(r"(?<![\w:])([A-Z][A-Za-z0-9_]*)\s*(?:\{|::|\()")


def strip_noise(code: str) -> str:
    """Removes comments, string literals and attributes.

    All three produce identifiers that look like calls: a sentence ending in
    "(see below)", a format string, or `#[derive(Debug)]`.
    """
    code = re.sub(r"#!?\[[^\]]*\]", "", code, flags=re.S)
    code = re.sub(r"//[^\n]*", "", code)
    code = re.sub(r'"(?:[^"\\]|\\.)*"', '""', code)
    return code


def blocks():
    for f in sorted(HERE.glob("*.md")):
        for m in re.finditer(r"```rust\n(.*?)```", f.read_text(), re.S):
            yield f.name, m.group(1)


def repo_symbols() -> set:
    if not REPO.is_dir():
        return set()
    # Private items count too. These books are about this workspace, so quoting
    # a `fn` that happens not to be `pub` is ordinary — and treating it as
    # undefined sends the reader looking for a gap that is not there.
    out = subprocess.run(
        ["grep", "-rhoE",
         r"(pub )?(fn|struct|enum|trait|type) [A-Za-z_][A-Za-z0-9_]*",
         "--include=*.rs", str(REPO)],
        capture_output=True, text=True).stdout
    return {line.split()[-1] for line in out.splitlines()}


def imported(code: str) -> set:
    """Identifiers brought in by `use`.

    A name imported from a crate is defined — just not by the book. Counting
    those as missing buries the ones that really are.
    """
    names = set()
    for clause in USE.findall(code):
        names |= set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", clause))
    return names


def variants(code: str) -> set:
    """Enum variant names, which are constructors rather than types."""
    names = set()
    for body in VARIANTS.findall(code):
        # Variants are comma-separated and several often share a line, as in
        # "Rx(f64), Ry(f64), Rz(f64),", so split rather than scan line starts.
        for part in body.split(","):
            m = re.match(r"\s*([A-Z][A-Za-z0-9_]*)", part)
            if m:
                names.add(m.group(1))
    return names


def check_undefined() -> list:
    defined, used_fn, used_ty, where = set(), {}, {}, {}
    defined_ty = set()
    for name, raw in blocks():
        code = strip_noise(raw)
        defined |= set(DEFN.findall(code)) | set(CLOSURE.findall(code))
        defined_ty |= imported(code) | variants(code)
        for a, b in TYPEDEF.findall(code):
            defined_ty.add(a or b)
    defined |= defined_ty          # a tuple struct is called like a function
    defined_ty |= PRELUDE | KEYWORDS
    known = KNOWN | KEYWORDS | PRELUDE | ELIDED | repo_symbols()
    defined_ty |= ELIDED
    globbed = {name for name, raw in blocks() if GLOB.search(raw)}
    for name, raw in blocks():
        code = strip_noise(raw)
        for sym in CALL.findall(code):
            used_fn.setdefault(sym, 0)
            used_fn[sym] += 1
            where.setdefault(sym, set()).add(name)
        if name in globbed:
            continue          # `use x::*` hides what it brought in
        for sym in TYPEUSE.findall(code):
            used_ty.setdefault(sym, 0)
            used_ty[sym] += 1
            where.setdefault(sym, set()).add(name)

    missing = []
    for sym, n in used_fn.items():
        if sym not in defined and sym not in known:
            missing.append((sym, n, sorted(where[sym])))
    for sym, n in used_ty.items():
        if sym not in defined_ty and sym not in known:
            missing.append((sym, n, sorted(where[sym])))
    return sorted(missing, key=lambda r: -r[1])


CHAPTER_BUILDS = _cfg.get("chapter_builds", {})

# Chapters whose code is checked for names rather than compiled. PyO3 macros
# need the crate and a Python build to expand, so the bindings chapter cannot
# be assembled here; what it can be held to is that every path it imports
# matches the crate layout the book sets out.
IMPORT_LAYOUT = {
    chapter: {crate: set(names) for crate, names in layout.items()}
    for chapter, layout in _cfg.get("import_layout", {}).items()
}


def check_imports() -> list:
    """Flags a name imported from the wrong crate.

    Chapter 10 lays out which crate holds what. Chapter 12 imported the parser
    and the backends from cforge_core, contradicting it, and nothing noticed
    because the chapter cannot be compiled.
    """
    problems = []
    for chapter, layout in IMPORT_LAYOUT.items():
        text = (HERE / chapter).read_text()
        placed = {name: crate for crate, names in layout.items() for name in names}
        for m in re.finditer(r"\buse\s+(cforge_[a-z]+)::\{?([^;]+?)\}?;", text):
            crate, clause = m.group(1), m.group(2)
            for name in re.findall(r"[A-Za-z_][A-Za-z0-9_]*", clause):
                if name == "as":
                    continue
                want = placed.get(name)
                if want and want != crate:
                    problems.append(
                        f"{chapter}: {name} imported from {crate}, "
                        f"but Chapter 10 puts it in {want}")
    return problems


def check_assembly() -> list:
    failures = []
    for chapter, wanted in CHAPTER_BUILDS.items():
        text = (HERE / chapter).read_text()
        found = [m.group(1) for m in re.finditer(r"```rust\n(.*?)```", text, re.S)]
        # One block often defines several of the things wanted, so collect the
        # matching indices and emit each block once, in document order.
        keep = set()
        for w in wanted:
            hit = [i for i, b in enumerate(found) if w in b]
            if not hit:
                failures.append(f"{chapter}: no block defines {w!r}")
            else:
                keep.add(hit[0])
        picked = [found[i] for i in sorted(keep)]
        src = "#![allow(unused, dead_code)]\n" + "\n".join(picked)
        with tempfile.TemporaryDirectory() as d:
            f = pathlib.Path(d) / "asm.rs"
            f.write_text(src)
            r = subprocess.run(
                ["rustc", "--edition", "2021", "--crate-type", "lib",
                 "--emit=metadata", "-o", str(pathlib.Path(d) / "asm.meta"), str(f)],
                capture_output=True, text=True)
            if r.returncode != 0:
                errs = [l for l in r.stderr.splitlines() if l.startswith("error")]
                failures.append(f"{chapter}: assembled code does not compile")
                failures += [f"    {e}" for e in errs[:8]]
    return failures


def main() -> int:
    bad = 0

    missing = check_undefined()
    if missing:
        bad = 1
        print("Symbols the book uses and never defines:\n")
        for sym, n, files in missing:
            print(f"  {n:>3}x  {sym:<32} {', '.join(files[:3])}")
    else:
        print("Undefined symbols: none")

    imports = check_imports()
    if imports:
        bad = 1
        print("\nImports that contradict the book's own crate layout:\n")
        for i in imports:
            print(f"  {i}")
    elif IMPORT_LAYOUT:
        print("Crate layout: imports agree with the layout the book states")
    else:
        print("Crate layout: not configured for this book")

    failures = check_assembly()
    if failures:
        bad = 1
        print("\nChapter assembly:\n")
        for f in failures:
            print(f"  {f}")
    elif CHAPTER_BUILDS:
        print(f"Chapter assembly: {len(CHAPTER_BUILDS)} chapter(s) compile")
    else:
        # Saying "0 chapters compile" as if it were a pass is the sort of green
        # that means nothing, which is the thing these books are about.
        print("Chapter assembly: not configured for this book")

    return bad


if __name__ == "__main__":
    sys.exit(main())
