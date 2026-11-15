"""One clean-machine pass for lanl/bngsim#1694, run in a fresh venv that holds
bngsim from PyPI. Each step prints one RESULT line; nothing here edits the repo.

Usage: python pass_1694.py STEP [args]
  wheel                 how bngsim was installed (wheel tag or a local build)
  examples REF          first Python block of README.md and docs/quickstart.md at REF, in an empty dir
  runs                  plain ODE and SSA runs of a two-species .net
  sens LABEL            a sensitivity_params run; LABEL says which compiler state this is
  bngl                  Model.from_bngl on a one-rule model (after `pip install 'bngsim[bngl]'`)
  extras                Model.from_antimony_string and result.dataframe (after the antimony and pandas extras)
"""

import math
import os
import platform
import re
import subprocess
import sys
import tempfile
import traceback
import urllib.request
from pathlib import Path

NET = """begin parameters
    1 k 0.5
end parameters
begin species
    1 A() 100
    2 B() 0
end species
begin reactions
    1 1 2 k
end reactions
begin groups
    1 A 1
    2 B 2
end groups
"""
FENCE = re.compile(r"^```(?:python|py)\s*$(.*?)^```\s*$", re.M | re.S)
HOST = (
    f"{platform.platform()} {platform.machine()} Python {platform.python_version()}"
    + (f" ({os.environ['PASS_HOST_NOTE']})" if os.environ.get("PASS_HOST_NOTE") else "")
)


def result(step, ok, detail):
    print(f"RESULT {step}: {'PASS' if ok else 'FAIL'} | {HOST} | {detail}", flush=True)


def last_line(exc):
    return "".join(traceback.format_exception_only(type(exc), exc)).strip().splitlines()[-1][:300]


def model_dir():
    d = Path(tempfile.mkdtemp())
    (d / "model.net").write_text(NET, encoding="utf-8")
    return d


def step_wheel():
    import bngsim

    dist = next(Path(bngsim.__file__).parent.parent.glob("bngsim-*.dist-info"))
    tags = [ln.split(":", 1)[1].strip() for ln in (dist / "WHEEL").read_text().splitlines() if ln.startswith("Tag:")]
    result("1 install", True, f"bngsim {bngsim.__version__}, wheel tag(s) {', '.join(tags)}")


def step_examples(ref):
    base = f"https://raw.githubusercontent.com/lanl/bngsim/{ref}/"
    for page in ("README.md", "docs/quickstart.md"):
        try:
            text = urllib.request.urlopen(base + page, timeout=30).read().decode("utf-8")
            block = FENCE.search(text).group(1)
        except Exception as exc:  # noqa: BLE001
            result(f"2 example {page}@{ref}", False, f"could not fetch: {last_line(exc)}")
            continue
        work = tempfile.mkdtemp()
        try:
            cp = subprocess.run([sys.executable, "-c", block], cwd=work, capture_output=True, text=True, timeout=600)
        except subprocess.TimeoutExpired:
            result(f"2 example {page}@{ref}", False, "did not finish in 600 s")
            continue
        if cp.returncode == 0:
            out = cp.stdout.strip().splitlines()
            result(f"2 example {page}@{ref}", True, f"runs; last line printed: {out[-1] if out else '(nothing)'}")
        else:
            result(f"2 example {page}@{ref}", False, cp.stderr.strip().splitlines()[-1][:300])


def step_runs():
    import bngsim

    d = model_dir()
    try:
        m = bngsim.Model.from_net(str(d / "model.net"))
        r = bngsim.Simulator(m, method="ode").run(t_span=(0.0, 10.0), n_points=101)
        a = float(r.observables["A"][-1])
        result("3 ODE run", abs(a - 100 * math.exp(-5)) < 1e-3, f"A(10) = {a:.6f}, closed form {100 * math.exp(-5):.6f}")
    except Exception as exc:  # noqa: BLE001
        result("3 ODE run", False, last_line(exc))
    try:
        m = bngsim.Model.from_net(str(d / "model.net"))
        r = bngsim.Simulator(m, method="ssa").run(t_span=(0.0, 10.0), n_points=11, seed=1)
        a, b = float(r.observables["A"][-1]), float(r.observables["B"][-1])
        result("3 SSA run", a + b == 100 and 0 <= a <= 100, f"seed 1: A(10) = {a:g}, B(10) = {b:g} (A + B = 100 conserved)")
    except Exception as exc:  # noqa: BLE001
        result("3 SSA run", False, last_line(exc))


def step_sens(label):
    import bngsim

    d = model_dir()
    cc = os.environ.get("CC", "(unset)")
    try:
        m = bngsim.Model.from_net(str(d / "model.net"))
        r = bngsim.Simulator(m, method="ode", sensitivity_params=["k"]).run(t_span=(0.0, 10.0), n_points=11)
        s = float(r.sensitivities[-1, 0, 0])
        want = -10.0 * 100.0 * math.exp(-5.0)  # dA(10)/dk = -t A0 e^(-k t)
        ok = abs(s - want) <= 1e-3 * abs(want)
        result(f"4 sensitivities, {label}", ok, f"CC={cc}: dA(10)/dk = {s:.6f}, closed form {want:.6f}")
    except Exception as exc:  # noqa: BLE001
        result(f"4 sensitivities, {label}", False, f"CC={cc}: {last_line(exc)}")


def step_bngl():
    import bngsim

    bngl = """begin model
begin parameters
    k 0.5
end parameters
begin molecule types
    A()
    B()
end molecule types
begin seed species
    A() 100
end seed species
begin observables
    Molecules A A()
end observables
begin reaction rules
    A() -> B() k
end reaction rules
end model
"""
    d = Path(tempfile.mkdtemp())
    (d / "m.bngl").write_text(bngl, encoding="utf-8")
    try:
        m = bngsim.Model.from_bngl(str(d / "m.bngl"))
        r = bngsim.Simulator(m, method="ode").run(t_span=(0.0, 10.0), n_points=101)
        a = float(r.observables["A"][-1])
        result("5 from_bngl", abs(a - 100 * math.exp(-5)) < 1e-3, f"A(10) = {a:.6f}")
    except Exception as exc:  # noqa: BLE001
        result("5 from_bngl", False, last_line(exc))


def step_extras():
    import bngsim

    try:
        m = bngsim.Model.from_antimony_string("A -> B; k*A; k = 0.5; A = 100; B = 0")
        r = bngsim.Simulator(m, method="ode").run(t_span=(0.0, 10.0), n_points=101)
        cols = list(r.species_names) if hasattr(r, "species_names") else None
        result("6 from_antimony_string", True, f"runs; species {cols}")
    except Exception as exc:  # noqa: BLE001
        result("6 from_antimony_string", False, last_line(exc))
        return
    try:
        df = r.dataframe
        result("6 result.dataframe", True, f"shape {df.shape}, columns {list(df.columns)[:4]}")
    except Exception as exc:  # noqa: BLE001
        result("6 result.dataframe", False, last_line(exc))


if __name__ == "__main__":
    step, *args = sys.argv[1:]
    {"wheel": step_wheel, "examples": step_examples, "runs": step_runs, "sens": step_sens, "bngl": step_bngl, "extras": step_extras}[step](*args)
