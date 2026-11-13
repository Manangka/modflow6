"""
Tests for the export list of the libmf6 shared library.

The export list is generated on every build by srcbmi/generate_exports.py
from the bind(C, name="...") entities in srcbmi, as a module-definition file
(Windows), linker version script (Linux) or exported symbols list (macOS).
These tests check the generator and that the built library exports exactly
the generated list.
"""

import importlib.util
import platform
import re
import shutil
import struct
import subprocess
from pathlib import Path

import pytest
from conftest import project_root_path

srcbmi_path = project_root_path / "srcbmi"

spec = importlib.util.spec_from_file_location(
    "generate_exports", srcbmi_path / "generate_exports.py"
)
generate_exports = importlib.util.module_from_spec(spec)
spec.loader.exec_module(generate_exports)

# entities only available in parallel builds
PARALLEL_ONLY = {"initialize_mpi"}

# symbols that are not part of MODFLOW but may remain visible in the library:
# defined by the GNU linker, or from the GCC runtime (libgcc) when it is linked
# statically, e.g. emulated thread-local storage and nested function trampolines
LINKER_SYMBOLS = {"_init", "_fini", "_edata", "_end", "__bss_start"}
RUNTIME_PREFIXES = ("__emutls_", "__gcc_")


def is_toolchain_symbol(name: str) -> bool:
    return name in LINKER_SYMBOLS or name.startswith(RUNTIME_PREFIXES)


def symbols(parallel: bool) -> dict[str, bool]:
    defines = {"__WITH_MPI__"} if parallel else set()
    syms = generate_exports.collect_symbols([str(srcbmi_path)], defines)
    return {s.name: s.is_data for s in syms}


def test_generate_exports_symbols():
    serial = symbols(parallel=False)
    parallel = symbols(parallel=True)
    assert set(parallel) - set(serial) == PARALLEL_ONLY
    assert set(serial) <= set(parallel)

    data = sorted(n for n, is_data in serial.items() if is_data)
    assert data == [
        "BMI_LENCOMPONENTNAME",
        "BMI_LENERRMESSAGE",
        "BMI_LENGRIDTYPE",
        "BMI_LENVARADDRESS",
        "BMI_LENVARTYPE",
        "BMI_LENVERSION",
        "ISTDOUTTOFILE",
    ]


PARSER_CASES = """\
module m
  integer(c_int), bind(C, name="var_one_line") :: v1 = 1
  integer(c_int), bind(C, name="var_split") &
    :: v2 = 2
  character(len=*), parameter :: s = 'no! bind(C, name="in_string")'
  type, bind(C) :: unnamed_type
  end type
contains
  function f1() result(s) bind(C, name="same_line")
  end function
  function f2(a, &
              b) &
    result(s) &
    ! comment inside a continued statement
    bind(C, name="continuation_line")
  end function
  function f3() result(s) bind(C, &
                               & name="split_inside_bind")
  end function
  function f4() result(s) !bind(C, name="commented_out")
  end function
  ! function f5() bind(C, name="comment_line")
end module
"""


def test_generate_exports_parser(tmp_path):
    src = tmp_path / "cases.f90"
    src.write_text(PARSER_CASES)
    found = {
        s.name: s.is_data for s in generate_exports.collect_symbols([str(src)], set())
    }
    assert found == {
        "var_one_line": True,
        "var_split": True,
        "same_line": False,
        "continuation_line": False,
        "split_inside_bind": False,
    }


@pytest.mark.parametrize(
    "code",
    [
        'subroutine f() bind(C, name=trim("f"))\nend subroutine\n',
        'subroutine f() bind(C, &\n#if X\n  name="f")\n#endif\nend subroutine\n',
        'subroutine f() bind(C, name="f") &\n',
    ],
    ids=["unparsable", "directive-in-statement", "unterminated"],
)
def test_generate_exports_parser_errors(tmp_path, code):
    src = tmp_path / "src.F90"
    src.write_text(code)
    with pytest.raises(SystemExit):
        generate_exports.collect_symbols([str(src)], set())


def generate(tmp_path: Path, fmt: str, *args: str) -> str:
    output = tmp_path / f"libmf6.{fmt}"
    generate_exports.main(
        ["--format", fmt, "--output", str(output), *args, str(srcbmi_path)]
    )
    return output.read_text()


def test_generate_exports_def(tmp_path):
    text = generate(tmp_path, "def")
    lines = [l.strip() for l in text.splitlines() if l.strip() and l[0] != ";"]
    assert lines[0] == "EXPORTS"
    exported = {}
    for line in lines[1:]:
        name, *attrs = line.split()
        exported[name] = attrs == ["DATA"]
    assert exported == symbols(parallel=False)


def test_generate_exports_map(tmp_path):
    text = generate(tmp_path, "map", "--define", "__WITH_MPI__")
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)
    match = re.fullmatch(r"\s*\{\s*global:(.*?)local:\s*\*;\s*\};\s*", text, re.S)
    assert match, text
    names = {n.strip() for n in match.group(1).split(";") if n.strip()}
    assert names == set(symbols(parallel=True))


def test_generate_exports_exports(tmp_path):
    text = generate(tmp_path, "exports")
    names = [l for l in text.splitlines() if l and l[0] != "#"]
    assert all(n.startswith("_") for n in names)
    assert {n[1:] for n in names} == set(symbols(parallel=False))


def test_generate_exports_if_changed(tmp_path):
    output = tmp_path / "libmf6.def"
    args = ["--format", "def", "--output", str(output), "--if-changed"]
    generate_exports.main([*args, str(srcbmi_path)])
    mtime = output.stat().st_mtime_ns
    generate_exports.main([*args, str(srcbmi_path)])
    assert output.stat().st_mtime_ns == mtime


@pytest.mark.parametrize(
    "directive", ["#if X > 1", "#elif defined(X)", "#if defined(X) && defined(Y)"]
)
def test_generate_exports_unsupported_directive(tmp_path, directive):
    src = tmp_path / "src.F90"
    src.write_text(
        f'{directive}\nsubroutine f() bind(C, name="f")\nend subroutine\n#endif\n'
    )
    with pytest.raises(SystemExit):
        generate_exports.collect_symbols([str(src)], set())


def read_dll_exports(path: Path) -> set[str]:
    """Read the exported symbol names from the export table of a PE file."""
    data = path.read_bytes()
    pe = struct.unpack_from("<I", data, 0x3C)[0]
    assert data[pe : pe + 4] == b"PE\0\0", f"{path} is not a PE file"
    coff = pe + 4
    nsections = struct.unpack_from("<H", data, coff + 2)[0]
    opt_size = struct.unpack_from("<H", data, coff + 16)[0]
    opt = coff + 20
    magic = struct.unpack_from("<H", data, opt)[0]
    datadir = opt + (112 if magic == 0x20B else 96)
    export_rva = struct.unpack_from("<I", data, datadir)[0]
    if export_rva == 0:
        return set()

    sections = []
    for i in range(nsections):
        sec = opt + opt_size + 40 * i
        vsize, vaddr, rsize, rptr = struct.unpack_from("<IIII", data, sec + 8)
        sections.append((vaddr, max(vsize, rsize), rptr))

    def offset(rva):
        for vaddr, size, rptr in sections:
            if vaddr <= rva < vaddr + size:
                return rva - vaddr + rptr
        raise ValueError(f"RVA {rva:#x} not in any section of {path}")

    exp = offset(export_rva)
    nnames = struct.unpack_from("<I", data, exp + 24)[0]
    names_rva = struct.unpack_from("<I", data, exp + 32)[0]
    names = set()
    for i in range(nnames):
        name_off = offset(struct.unpack_from("<I", data, offset(names_rva) + 4 * i)[0])
        names.add(data[name_off : data.index(b"\0", name_off)].decode())
    return names


def read_nm_exports(path: Path, args: list[str], strip_underscore: bool) -> set[str]:
    """Read the exported (defined, global) symbol names with nm."""
    nm = shutil.which("nm")
    if nm is None:
        pytest.skip("nm not available")
    output = subprocess.check_output([nm, *args, str(path)], text=True)
    names = set()
    for line in output.splitlines():
        parts = line.split()
        if len(parts) < 3:
            continue
        name = parts[-1].split("@", 1)[0]
        if strip_underscore:
            name = name.removeprefix("_")
        names.add(name)
    return names


def libmf6_exports(targets) -> set[str]:
    lib = targets["libmf6"]
    if not lib.is_file():
        pytest.skip(f"{lib} not found")
    system = platform.system()
    if system == "Windows":
        return read_dll_exports(lib)
    if system == "Darwin":
        return read_nm_exports(lib, ["-gU"], strip_underscore=True)
    return read_nm_exports(lib, ["-D", "--defined-only"], strip_underscore=False)


def test_libmf6_dll_exports(targets):
    exports = libmf6_exports(targets)
    expected = set(symbols(parallel=False))
    missing = sorted(expected - exports)
    unexpected = sorted(
        n for n in exports - expected - PARALLEL_ONLY if not is_toolchain_symbol(n)
    )
    assert not missing, f"libmf6 does not export: {missing}"
    assert not unexpected, f"libmf6 exports names not in its C API: {unexpected}"


@pytest.mark.parallel
def test_par_libmf6_dll_exports(targets):
    exports = libmf6_exports(targets)
    missing = sorted(set(symbols(parallel=True)) - exports)
    assert not missing, f"parallel libmf6 does not export: {missing}"
