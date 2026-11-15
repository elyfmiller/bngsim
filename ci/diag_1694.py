"""Diagnose the two Linux failures of the lanl/bngsim#1694 pass on python:3.13-slim.

Each check prints its full exception, not only the last line:
  sbml        import libsbml, then Model.from_sbml on a one-reaction SBML file
  antimony    Model.from_antimony_string
  bngl        Model.from_bngl on a one-rule model
  dataframe   result.dataframe on a .net run (independent of antimony)
  libs        ldd of every shared object in the libsbml and antimony packages,
              filtered to missing libraries and expat
"""

from __future__ import annotations

import glob
import importlib.util
import os
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

NET = """begin parameters
    1 k 0.5
end parameters
begin species
    1 A() 100
    2 B() 0
end species
begin reactions
    1 1 2 k #_R1
end reactions
begin groups
    1 A 1
    2 B 2
end groups
"""

BNGL = """begin model
begin parameters
    k 0.5
end parameters
begin molecule types
    A()
    B()
end molecule types
begin seed species
    A() 100
    B() 0
end seed species
begin observables
    Molecules A A()
end observables
begin reaction rules
    A() -> B() k
end reaction rules
end model
"""

SBML = """<?xml version="1.0" encoding="UTF-8"?>
<sbml xmlns="http://www.sbml.org/sbml/level3/version2/core" level="3" version="2">
  <model id="m">
    <listOfCompartments>
      <compartment id="c" size="1" constant="true"/>
    </listOfCompartments>
    <listOfSpecies>
      <species id="A" compartment="c" initialAmount="100" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
      <species id="B" compartment="c" initialAmount="0" hasOnlySubstanceUnits="true" boundaryCondition="false" constant="false"/>
    </listOfSpecies>
    <listOfParameters>
      <parameter id="k" value="0.5" constant="true"/>
    </listOfParameters>
    <listOfReactions>
      <reaction id="R1" reversible="false">
        <listOfReactants><speciesReference species="A" stoichiometry="1" constant="true"/></listOfReactants>
        <listOfProducts><speciesReference species="B" stoichiometry="1" constant="true"/></listOfProducts>
        <kineticLaw>
          <math xmlns="http://www.w3.org/1998/Math/MathML">
            <apply><times/><ci>k</ci><ci>A</ci></apply>
          </math>
        </kineticLaw>
      </reaction>
    </listOfReactions>
  </model>
</sbml>
"""

ANT = "A -> B; k*A\nk = 0.5\nA = 100\nB = 0\n"


def label() -> str:
    return os.environ.get("DIAG_LABEL", "")


def report(name, fn):
    print(f"----- DIAG {name} [{label()}]", flush=True)
    try:
        detail = fn()
        print(f"DIAG {name} [{label()}]: OK | {detail}", flush=True)
    except BaseException as exc:  # noqa: BLE001 - every failure is the finding
        traceback.print_exception(exc, file=sys.stdout)
        print(f"DIAG {name} [{label()}]: FAIL | {type(exc).__name__}", flush=True)


def a_end(result) -> float:
    return float(result.observables["A"][-1])


def check_sbml():
    import libsbml

    out = f"libsbml {libsbml.getLibSBMLDottedVersion()}"
    import bngsim

    d = Path(tempfile.mkdtemp())
    (d / "m.xml").write_text(SBML)
    model = bngsim.Model.from_sbml(str(d / "m.xml"))
    r = bngsim.Simulator(model).run(t_span=(0, 10), n_points=11)
    return f"{out}; from_sbml A(10) = {float(r.species[-1, 0]):.6f}"


def check_antimony():
    import bngsim

    model = bngsim.Model.from_antimony_string(ANT)
    r = bngsim.Simulator(model).run(t_span=(0, 10), n_points=11)
    return f"A(10) = {float(r.species[-1, 0]):.6f}"


def check_bngl():
    import bngsim

    d = Path(tempfile.mkdtemp())
    (d / "m.bngl").write_text(BNGL)
    model = bngsim.Model.from_bngl(str(d / "m.bngl"))
    return f"A(10) = {a_end(bngsim.Simulator(model).run(t_span=(0, 10), n_points=11)):.6f}"


def check_dataframe():
    import bngsim

    d = Path(tempfile.mkdtemp())
    (d / "m.net").write_text(NET)
    r = bngsim.Simulator(bngsim.Model.from_net(str(d / "m.net"))).run(t_span=(0, 10), n_points=11)
    df = r.dataframe
    return f"shape {df.shape}, columns {list(df.columns)[:4]}"


def check_libs():
    lines = []
    for pkg in ("libsbml", "antimony"):
        spec = importlib.util.find_spec(pkg)
        if spec is None:
            lines.append(f"{pkg}: not installed")
            continue
        roots = list(spec.submodule_search_locations or [os.path.dirname(spec.origin)])
        sos = sorted(
            {p for r in roots for p in glob.glob(os.path.join(r, "**", "*.so*"), recursive=True)}
        )
        sos += sorted(glob.glob(os.path.join(os.path.dirname(roots[0]), f"*{pkg}*.libs", "*.so*")))
        for so in sos:
            cp = subprocess.run(["ldd", so], capture_output=True, text=True)
            hits = [
                ln.strip() for ln in cp.stdout.splitlines() if "not found" in ln or "expat" in ln
            ]
            lines.append(
                f"{pkg}: {os.path.relpath(so, os.path.dirname(roots[0]))}: {'; '.join(hits) or 'no missing libraries, no expat'}"
            )
    print("\n".join(lines))
    return f"{len(lines)} shared objects listed above"


CHECKS = {
    "sbml": check_sbml,
    "antimony": check_antimony,
    "bngl": check_bngl,
    "dataframe": check_dataframe,
    "libs": check_libs,
}

if __name__ == "__main__":
    for name in sys.argv[1:] or list(CHECKS):
        report(name, CHECKS[name])
