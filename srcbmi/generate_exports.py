"""
Generate the export list of the MODFLOW 6 shared library (libmf6).

Every entity in the srcbmi sources with a bind(C, name="...") attribute
is part of the public C API and is exported. The list is written in the
format of the target platform's linker:

  def      Windows module-definition file (link.exe, xilink, MinGW ld)
  map      Linux linker version script (ld --version-script)
  exports  macOS exported symbols list (ld -exported_symbols_list)

Entities inside preprocessor blocks (e.g. #if defined(__WITH_MPI__)) are
only exported when the corresponding macros are passed with --define.

Usage:
  python generate_exports.py --format def --output libmf6.def [--define X] SRC...

SRC can be source files or directories, which are searched for *.f90/*.F90.
This script is run by meson (srcbmi/meson.build) and by the Visual Studio
project (msvs/mf6bmi.vfproj) on every build.
"""

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

BIND_C = re.compile(r'bind\s*\(\s*c\s*,\s*name\s*=\s*"(\w+)"\s*\)(\s*::)?', re.I)
# bind(C, ...) with arguments; bind(C) without a name is not exported
BIND_C_ARGS = re.compile(r"bind\s*\(\s*c\s*,", re.I)
PP_IF = re.compile(
    r"^\s*#\s*(?:if\s+(!?)\s*defined\s*\(?\s*(\w+)\s*\)?|(ifdef|ifndef)\s+(\w+))\s*$"
)
PP_ELSE = re.compile(r"^\s*#\s*else\b")
PP_ENDIF = re.compile(r"^\s*#\s*endif\b")
PP_ANY = re.compile(r"^\s*#\s*(if|ifdef|ifndef|elif|else|endif)\b")


class Symbol:
    def __init__(self, name: str, is_data: bool, source: str):
        self.name = name
        self.is_data = is_data
        self.source = source


def collect_sources(paths: list[str]) -> list[Path]:
    files = []
    for p in map(Path, paths):
        if p.is_dir():
            files.extend(f for f in p.iterdir() if f.suffix.lower() == ".f90")
        else:
            files.append(p)
    return sorted(set(files), key=lambda f: f.name.lower())


def string_mask(code: str) -> list[bool]:
    """Return per character whether it is inside a string literal."""
    mask = []
    quote = None
    for ch in code:
        if quote:
            mask.append(True)
            if ch == quote:
                quote = None
        else:
            if ch in "'\"":
                quote = ch
            mask.append(quote is not None)
    return mask


def strip_comment(line: str) -> str:
    """Remove a trailing ! comment, ignoring ! inside string literals."""
    for i, (ch, in_string) in enumerate(zip(line, string_mask(line))):
        if ch == "!" and not in_string:
            return line[:i]
    return line


def parse_statement(stmt: str, where: str) -> list[Symbol]:
    mask = string_mask(stmt)
    symbols = [
        Symbol(m.group(1), m.group(2) is not None, where)
        for m in BIND_C.finditer(stmt)
        if not mask[m.start()]
    ]
    nargs = sum(not mask[m.start()] for m in BIND_C_ARGS.finditer(stmt))
    if len(symbols) != nargs:
        raise SystemExit(
            f'{where}: cannot parse bind(C, ...), expected bind(C, name="..."): '
            f"{' '.join(stmt.split())}"
        )
    return symbols


def parse_source(path: Path, defines: set[str]) -> list[Symbol]:
    """Return the bind(C) entities that are active for the given defines.

    Comments are ignored and continuation lines (&) are joined, so bind(C)
    may be written on, or split over, any line of a statement.
    """
    symbols = []
    stack = []  # one bool per open preprocessor block: is the branch active
    pending = ""  # code of a statement continued on the next line
    start = 0  # line number where the pending statement started
    for lineno, line in enumerate(path.read_text().splitlines(), 1):
        where = f"{path.name}:{lineno}"
        if pending and PP_ANY.match(line):
            raise SystemExit(f"{where}: preprocessor directive inside a statement")
        m = PP_IF.match(line)
        if m:
            negate, macro = (
                m.group(1) or m.group(3) == "ifndef",
                m.group(2) or m.group(4),
            )
            stack.append((macro in defines) != bool(negate))
            continue
        if PP_ELSE.match(line):
            if not stack:
                raise SystemExit(f"{where}: #else without #if")
            stack[-1] = not stack[-1]
            continue
        if PP_ENDIF.match(line):
            if not stack:
                raise SystemExit(f"{where}: #endif without #if")
            stack.pop()
            continue
        if PP_ANY.match(line):
            raise SystemExit(f"{where}: unsupported preprocessor directive: {line}")
        if not all(stack):
            continue
        code = strip_comment(line).strip()
        if pending:
            if not code:
                continue  # comment or blank line inside a continued statement
            code = code.removeprefix("&")
        else:
            start = lineno
        if code.endswith("&"):
            pending += code[:-1] + " "
            continue
        symbols += parse_statement(pending + code, f"{path.name}:{start}")
        pending = ""
    if pending:
        raise SystemExit(f"{path.name}:{start}: unterminated continued statement")
    if stack:
        raise SystemExit(f"{path.name}: unterminated #if block")
    return symbols


def collect_symbols(paths: list[str], defines: set[str]) -> list[Symbol]:
    symbols = []
    seen = {}
    for path in collect_sources(paths):
        for sym in parse_source(path, defines):
            if sym.name in seen:
                where = f"{seen[sym.name]}, {sym.source}"
                raise SystemExit(f"duplicate bind(C) name '{sym.name}': {where}")
            seen[sym.name] = sym.source
            symbols.append(sym)
    if not symbols:
        raise SystemExit("no bind(C, name=...) entities found")
    return symbols


HEADER = "Public C API of libmf6, generated by srcbmi/generate_exports.py"


def format_def(symbols: list[Symbol]) -> str:
    lines = [f"; {HEADER}", "EXPORTS"]
    lines += [f"    {s.name}{' DATA' if s.is_data else ''}" for s in symbols]
    return "\n".join(lines) + "\n"


def format_map(symbols: list[Symbol]) -> str:
    lines = [f"/* {HEADER} */", "{", "  global:"]
    lines += [f"    {s.name};" for s in symbols]
    lines += ["  local:", "    *;", "};"]
    return "\n".join(lines) + "\n"


def format_exports(symbols: list[Symbol]) -> str:
    lines = [f"# {HEADER}"] + [f"_{s.name}" for s in symbols]
    return "\n".join(lines) + "\n"


FORMATS = {"def": format_def, "map": format_map, "exports": format_exports}


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument("--format", choices=sorted(FORMATS), required=True)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--define", action="append", default=[], help="preprocessor macro to define"
    )
    parser.add_argument(
        "--if-changed",
        action="store_true",
        help="only write the output if its contents change",
    )
    parser.add_argument("sources", nargs="+", help="source files or directories")
    args = parser.parse_args(argv)

    text = FORMATS[args.format](collect_symbols(args.sources, set(args.define)))
    if args.if_changed and args.output.is_file() and args.output.read_text() == text:
        return 0
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
