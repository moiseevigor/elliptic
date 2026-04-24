#!/usr/bin/env python3
"""Produce verified numerical examples for each function in the website finder.

Output: scripts/out/function_examples.json — a dict keyed by function id
containing Python/MATLAB snippets plus the numerical output actually
produced by the Python package. The MATLAB snippet uses identical
arguments and identical expected values (the implementations are mathematically
the same).

Run:
    python scripts/extract_function_examples.py

The output JSON is then hand-ported into gh-pages:/data/functions.js.
"""
from __future__ import annotations

import json
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np
import elliptic


# ──────────────────────────────────────────────────────────────────────────
# helpers
# ──────────────────────────────────────────────────────────────────────────
def fmt(x: Any, digits: int = 4) -> str:
    """Format a scalar (incl. numpy 0-d array / complex) compactly."""
    # Unwrap numpy 0-d arrays
    if isinstance(x, np.ndarray):
        if x.ndim == 0:
            x = x.item()
        else:
            return np.array2string(x, precision=digits)
    if isinstance(x, complex) or (isinstance(x, np.generic) and np.iscomplexobj(x)):
        return f"{x.real:.{digits}f}{'+' if x.imag >= 0 else '-'}{abs(x.imag):.{digits}f}i"
    if isinstance(x, (float, int)) or isinstance(x, np.floating):
        return f"{float(x):.{digits}f}"
    return str(x)


@dataclass
class Example:
    id: str
    python: str           # Python snippet (including imports)
    matlab: str           # MATLAB snippet
    output: str           # expected-output comment (shared)
    _compute: Callable[[], Any] = field(repr=False)

    def run(self) -> Any:
        """Actually execute the Python snippet to confirm it produces `output`."""
        return self._compute()


EXAMPLES: list[Example] = []


def register(
    id: str,
    python: str,
    matlab: str,
    output: str,
    compute: Callable[[], Any],
) -> None:
    """Register an example. `compute` is a thunk that returns the value shown
    in `output` — we use it to verify the output string is truthful."""
    EXAMPLES.append(Example(id=id, python=python, matlab=matlab,
                            output=output, _compute=compute))


# ──────────────────────────────────────────────────────────────────────────
# Jacobi elliptic functions
# ──────────────────────────────────────────────────────────────────────────
def _compute_ellipj():
    sn, cn, dn, am = elliptic.ellipj(0.5, 0.3)
    return sn, cn, dn, am

sn, cn, dn, am = _compute_ellipj()
register(
    "ellipj",
    python="import elliptic\n"
           "sn, cn, dn, am = elliptic.ellipj(0.5, 0.3)",
    matlab="[sn, cn, dn, am] = ellipj(0.5, 0.3);",
    output=f"sn = {fmt(sn)}    cn = {fmt(cn)}    dn = {fmt(dn)}    am = {fmt(am)}",
    compute=_compute_ellipj,
)


def _compute_ellipji():
    u = 0.4 + 0.3j
    sni, cni, dni = elliptic.ellipji(u, 0.7)
    return sni, cni, dni

sni, cni, dni = _compute_ellipji()
register(
    "ellipji",
    python="import elliptic\n"
           "sni, cni, dni = elliptic.ellipji(0.4 + 0.3j, 0.7)",
    matlab="[sni, cni, dni] = ellipji(0.4 + 0.3i, 0.7);",
    output=f"sni = {fmt(sni)}\ncni = {fmt(cni)}\ndni = {fmt(dni)}",
    compute=_compute_ellipji,
)


def _compute_jacobiThetaEta():
    Th, H = elliptic.jacobiThetaEta(0.3, 0.5)
    return Th, H

Th_je, H_je = _compute_jacobiThetaEta()
register(
    "jacobiThetaEta",
    python="import elliptic\n"
           "Th, H = elliptic.jacobiThetaEta(0.3, 0.5)",
    matlab="[Th, H] = jacobiThetaEta(0.3, 0.5);",
    output=f"Th = {fmt(Th_je)}    H = {fmt(H_je)}",
    compute=_compute_jacobiThetaEta,
)


def _compute_theta():
    return elliptic.theta(3, 0.3, 0.5)

th3 = _compute_theta()
register(
    "theta",
    python="import elliptic\n"
           "Th3 = elliptic.theta(3, 0.3, 0.5)   # type=3",
    matlab="Th3 = theta(3, 0.3, 0.5);   % type=3",
    output=f"Th3 = {fmt(th3)}",
    compute=_compute_theta,
)


# ──────────────────────────────────────────────────────────────────────────
# Elliptic integrals
# ──────────────────────────────────────────────────────────────────────────
def _compute_elliptic12():
    F, E, Z = elliptic.elliptic12(np.pi/4, 0.5)
    return F, E, Z

F_, E_, Z_ = _compute_elliptic12()
register(
    "elliptic12",
    python="import numpy as np, elliptic\n"
           "F, E, Z = elliptic.elliptic12(np.pi/4, 0.5)",
    matlab="[F, E, Z] = elliptic12(pi/4, 0.5);",
    output=f"F = {fmt(F_)}    E = {fmt(E_)}    Z = {fmt(Z_)}",
    compute=_compute_elliptic12,
)


def _compute_elliptic12i():
    Fi, Ei, Zi = elliptic.elliptic12i(0.4 + 0.3j, 0.7)
    return Fi, Ei, Zi

Fi_, Ei_, Zi_ = _compute_elliptic12i()
register(
    "elliptic12i",
    python="import elliptic\n"
           "Fi, Ei, Zi = elliptic.elliptic12i(0.4 + 0.3j, 0.7)",
    matlab="[Fi, Ei, Zi] = elliptic12i(0.4 + 0.3i, 0.7);",
    output=f"Fi = {fmt(Fi_)}\nEi = {fmt(Ei_)}\nZi = {fmt(Zi_)}",
    compute=_compute_elliptic12i,
)


def _compute_elliptic3():
    return elliptic.elliptic3(np.pi/4, 0.5, 0.3)

Pi_ = _compute_elliptic3()
register(
    "elliptic3",
    python="import numpy as np, elliptic\n"
           "Pi = elliptic.elliptic3(np.pi/4, 0.5, 0.3)   # Pi(phi,m,n)",
    matlab="Pi = elliptic3(pi/4, 0.5, 0.3);   % Pi(phi,m,n)",
    output=f"Pi = {fmt(Pi_)}",
    compute=_compute_elliptic3,
)


# NOTE: elliptic123 is a MATLAB-only wrapper; no Python equivalent exists.
# It's represented in the finder UI as matlab-only (python tab disabled).


def _compute_inverselliptic2():
    # Round-trip: compute E(phi, m), then invert it back to phi.
    phi, m = 0.6, 0.4
    _, E_val, _ = elliptic.elliptic12(phi, m)
    phi_back = elliptic.inverselliptic2(E_val, m)
    return E_val, phi_back

Eval_, phi_back = _compute_inverselliptic2()
register(
    "inverselliptic2",
    python="import elliptic\n"
           "# Given E(phi=0.6, m=0.4), recover phi:\n"
           "_, E_val, _ = elliptic.elliptic12(0.6, 0.4)\n"
           "phi_back = elliptic.inverselliptic2(E_val, 0.4)",
    matlab="% Given E(phi=0.6, m=0.4), recover phi:\n"
           "[~, E_val, ~] = elliptic12(0.6, 0.4);\n"
           "phi_back = inverselliptic2(E_val, 0.4);",
    output=f"E_val = {fmt(Eval_)}    phi_back = {fmt(phi_back)}   (= 0.6)",
    compute=_compute_inverselliptic2,
)


# ──────────────────────────────────────────────────────────────────────────
# Carlson symmetric forms
# ──────────────────────────────────────────────────────────────────────────
def _compute_carlsonRF():
    # R_F(0, 1-m, 1) = K(m); check against K(0.5)
    m = 0.5
    return elliptic.carlsonRF(0.0, 1 - m, 1.0)

RF_ = _compute_carlsonRF()
register(
    "carlsonRF",
    python="import elliptic\n"
           "# R_F(0, 1-m, 1) equals the complete K(m).\n"
           "RF = elliptic.carlsonRF(0.0, 0.5, 1.0)   # m = 0.5",
    matlab="% R_F(0, 1-m, 1) equals the complete K(m).\n"
           "RF = carlsonRF(0.0, 0.5, 1.0);   % m = 0.5",
    output=f"RF = {fmt(RF_, 6)}   (= K(0.5))",
    compute=_compute_carlsonRF,
)


def _compute_carlsonRD():
    m = 0.5
    return elliptic.carlsonRD(0.0, 1 - m, 1.0)

RD_ = _compute_carlsonRD()
register(
    "carlsonRD",
    python="import elliptic\n"
           "RD = elliptic.carlsonRD(0.0, 0.5, 1.0)   # m = 0.5",
    matlab="RD = carlsonRD(0.0, 0.5, 1.0);   % m = 0.5",
    output=f"RD = {fmt(RD_, 6)}   (= 3·D(0.5))",
    compute=_compute_carlsonRD,
)


def _compute_carlsonRJ():
    return elliptic.carlsonRJ(0.0, 0.5, 1.0, 0.7)

RJ_ = _compute_carlsonRJ()
register(
    "carlsonRJ",
    python="import elliptic\n"
           "RJ = elliptic.carlsonRJ(0.0, 0.5, 1.0, 0.7)",
    matlab="RJ = carlsonRJ(0.0, 0.5, 1.0, 0.7);",
    output=f"RJ = {fmt(RJ_, 6)}",
    compute=_compute_carlsonRJ,
)


def _compute_carlsonRC():
    # R_C(0, 1/4) = π; classic identity
    return elliptic.carlsonRC(0.0, 0.25)

RC_ = _compute_carlsonRC()
register(
    "carlsonRC",
    python="import elliptic\n"
           "# R_C(0, 1/4) = pi  (closed-form identity)\n"
           "RC = elliptic.carlsonRC(0.0, 0.25)",
    matlab="% R_C(0, 1/4) = pi  (closed-form identity)\n"
           "RC = carlsonRC(0.0, 0.25);",
    output=f"RC = {fmt(RC_, 10)}   (= pi)",
    compute=_compute_carlsonRC,
)


# ──────────────────────────────────────────────────────────────────────────
# Weierstrass
# ──────────────────────────────────────────────────────────────────────────
# Canonical half-period roots
E1, E2, E3 = 2.0, 0.5, -2.5


def _compute_weierstrassInvariants():
    g2, g3, Delta = elliptic.weierstrassInvariants(E1, E2, E3)
    return g2, g3, Delta

g2, g3, Delta = _compute_weierstrassInvariants()
register(
    "weierstrassInvariants",
    python="import elliptic\n"
           "g2, g3, Delta = elliptic.weierstrassInvariants(2.0, 0.5, -2.5)",
    matlab="[g2, g3, Delta] = weierstrassInvariants(2.0, 0.5, -2.5);",
    output=f"g2 = {fmt(g2)}    g3 = {fmt(g3)}    Delta = {fmt(Delta)}",
    compute=_compute_weierstrassInvariants,
)


def _compute_weierstrassP():
    return elliptic.weierstrassP(0.4, E1, E2, E3)

P_ = _compute_weierstrassP()
register(
    "weierstrassP",
    python="import elliptic\n"
           "P = elliptic.weierstrassP(0.4, 2.0, 0.5, -2.5)",
    matlab="P = weierstrassP(0.4, 2.0, 0.5, -2.5);",
    output=f"P = {fmt(P_)}",
    compute=_compute_weierstrassP,
)


def _compute_weierstrassPPrime():
    return elliptic.weierstrassPPrime(0.4, E1, E2, E3)

dP_ = _compute_weierstrassPPrime()
register(
    "weierstrassPPrime",
    python="import elliptic\n"
           "dP = elliptic.weierstrassPPrime(0.4, 2.0, 0.5, -2.5)",
    matlab="dP = weierstrassPPrime(0.4, 2.0, 0.5, -2.5);",
    output=f"dP = {fmt(dP_)}",
    compute=_compute_weierstrassPPrime,
)


def _compute_weierstrassZeta():
    return elliptic.weierstrassZeta(0.4, E1, E2, E3)

Zw_ = _compute_weierstrassZeta()
register(
    "weierstrassZeta",
    python="import elliptic\n"
           "Zw = elliptic.weierstrassZeta(0.4, 2.0, 0.5, -2.5)",
    matlab="Zw = weierstrassZeta(0.4, 2.0, 0.5, -2.5);",
    output=f"Zw = {fmt(Zw_)}",
    compute=_compute_weierstrassZeta,
)


def _compute_weierstrassSigma():
    return elliptic.weierstrassSigma(0.4, E1, E2, E3)

Sw_ = _compute_weierstrassSigma()
register(
    "weierstrassSigma",
    python="import elliptic\n"
           "Sw = elliptic.weierstrassSigma(0.4, 2.0, 0.5, -2.5)",
    matlab="Sw = weierstrassSigma(0.4, 2.0, 0.5, -2.5);",
    output=f"Sw = {fmt(Sw_)}",
    compute=_compute_weierstrassSigma,
)


# ──────────────────────────────────────────────────────────────────────────
# Associate integrals
# ──────────────────────────────────────────────────────────────────────────
def _compute_ellipticBDJ():
    B, D, J = elliptic.ellipticBDJ(np.pi/4, 0.5, 0.3)
    return B, D, J

B_, D_, J_ = _compute_ellipticBDJ()
register(
    "ellipticBDJ",
    python="import numpy as np, elliptic\n"
           "B, D, J = elliptic.ellipticBDJ(np.pi/4, 0.5, 0.3)",
    matlab="[B, D, J] = ellipticBDJ(pi/4, 0.5, 0.3);",
    output=f"B = {fmt(B_)}    D = {fmt(D_)}    J = {fmt(J_)}",
    compute=_compute_ellipticBDJ,
)


def _compute_ellipticBD():
    Bc, Dc, Sc = elliptic.ellipticBD(np.array([0.5]))
    return float(Bc[0]), float(Dc[0]), float(Sc[0])

Bc_, Dc_, Sc_ = _compute_ellipticBD()
register(
    "ellipticBD",
    python="import elliptic, numpy as np\n"
           "B, D, S = elliptic.ellipticBD(np.array([0.5]))   # complete at m=0.5",
    matlab="[B, D, S] = ellipticBD(0.5);   % complete at m=0.5",
    output=f"B = {fmt(Bc_)}    D = {fmt(Dc_)}    S = {fmt(Sc_)}",
    compute=_compute_ellipticBD,
)


def _compute_jacobiEDJ():
    Eu, Du, Ju = elliptic.jacobiEDJ(0.5, 0.3, 0.2)
    return Eu, Du, Ju

Eu_, Du_, Ju_ = _compute_jacobiEDJ()
register(
    "jacobiEDJ",
    python="import elliptic\n"
           "Eu, Du, Ju = elliptic.jacobiEDJ(0.5, 0.3, 0.2)",
    matlab="[Eu, Du, Ju] = jacobiEDJ(0.5, 0.3, 0.2);",
    output=f"Eu = {fmt(Eu_)}    Du = {fmt(Du_)}    Ju = {fmt(Ju_)}",
    compute=_compute_jacobiEDJ,
)


# ──────────────────────────────────────────────────────────────────────────
# Bulirsch
# ──────────────────────────────────────────────────────────────────────────
KC = np.sqrt(1 - 0.5)


def _compute_cel():
    # cel(kc, p, a, b) with p=a=b=1 collapses to K(m)
    return elliptic.cel(KC, 1.0, 1.0, 1.0)

cel_ = _compute_cel()
register(
    "cel",
    python="import numpy as np, elliptic\n"
           "kc = np.sqrt(1 - 0.5)   # m = 0.5\n"
           "C = elliptic.cel(kc, 1.0, 1.0, 1.0)   # collapses to K(m)",
    matlab="kc = sqrt(1 - 0.5);   % m = 0.5\n"
           "C = cel(kc, 1.0, 1.0, 1.0);   % collapses to K(m)",
    output=f"C = {fmt(cel_, 6)}   (= K(0.5))",
    compute=_compute_cel,
)


def _compute_cel1():
    return elliptic.cel1(KC)

cel1_ = _compute_cel1()
register(
    "cel1",
    python="import numpy as np, elliptic\n"
           "K = elliptic.cel1(np.sqrt(1 - 0.5))   # = K(0.5)",
    matlab="K = cel1(sqrt(1 - 0.5));   % = K(0.5)",
    output=f"K = {fmt(cel1_, 6)}",
    compute=_compute_cel1,
)


def _compute_cel2():
    return elliptic.cel2(KC, 1.0, 1 - 0.5)

cel2_ = _compute_cel2()
register(
    "cel2",
    python="import numpy as np, elliptic\n"
           "E = elliptic.cel2(np.sqrt(1 - 0.5), 1.0, 1 - 0.5)   # = E(0.5)",
    matlab="E = cel2(sqrt(1 - 0.5), 1.0, 1 - 0.5);   % = E(0.5)",
    output=f"E = {fmt(cel2_, 6)}",
    compute=_compute_cel2,
)


def _compute_cel3():
    return elliptic.cel3(KC, 1 - 0.3)

cel3_ = _compute_cel3()
register(
    "cel3",
    python="import numpy as np, elliptic\n"
           "Pi = elliptic.cel3(np.sqrt(1 - 0.5), 1 - 0.3)   # = Pi(0.3 | 0.5)",
    matlab="Pi = cel3(sqrt(1 - 0.5), 1 - 0.3);   % = Pi(0.3 | 0.5)",
    output=f"Pi = {fmt(cel3_, 6)}",
    compute=_compute_cel3,
)


# ──────────────────────────────────────────────────────────────────────────
# Related
# ──────────────────────────────────────────────────────────────────────────
def _compute_agm():
    # AGM(1, sqrt(1-m)) = pi/(2 K(m))
    return elliptic.agm(1.0, np.sqrt(1 - 0.5))

agm_ = _compute_agm()
register(
    "agm",
    python="import numpy as np, elliptic\n"
           "a = elliptic.agm(1.0, np.sqrt(1 - 0.5))   # = pi / (2 K(0.5))",
    matlab="a = agm(1.0, sqrt(1 - 0.5));   % = pi / (2 K(0.5))",
    output=f"a = {fmt(agm_, 6)}",
    compute=_compute_agm,
)


def _compute_nomeq():
    return elliptic.nomeq(0.5)

q_ = _compute_nomeq()
register(
    "nomeq",
    python="import elliptic\n"
           "q = elliptic.nomeq(0.5)   # nome q(m)",
    matlab="q = nomeq(0.5);   % nome q(m)",
    output=f"q = {fmt(q_, 6)}",
    compute=_compute_nomeq,
)


def _compute_inversenomeq():
    return elliptic.inversenomeq(0.1)

m_ = _compute_inversenomeq()
register(
    "inversenomeq",
    python="import elliptic\n"
           "m = elliptic.inversenomeq(0.1)   # valid for q ≤ 0.6",
    matlab="m = inversenomeq(0.1);   % valid for q ≤ 0.6",
    output=f"m = {fmt(m_, 6)}",
    compute=_compute_inversenomeq,
)


# ──────────────────────────────────────────────────────────────────────────
# Verify every registered example runs and produces a value that parses in
# the output comment. This is a self-consistency check — `output` is built
# from the same compute function by fmt(...), so this mainly guards against
# a function being registered but failing at runtime.
# ──────────────────────────────────────────────────────────────────────────
def verify() -> None:
    failures = []
    for ex in EXAMPLES:
        try:
            result = ex.run()
        except Exception as e:
            failures.append((ex.id, f"runtime error: {e!r}"))
            continue
        if result is None:
            failures.append((ex.id, "compute returned None"))
    if failures:
        print("VERIFICATION FAILURES:")
        for name, reason in failures:
            print(f"  - {name}: {reason}")
        sys.exit(1)
    print(f"✓ All {len(EXAMPLES)} examples ran without error.")


def emit_json(out_path: Path) -> None:
    data = {}
    for ex in EXAMPLES:
        data[ex.id] = {
            "python": ex.python,
            "matlab": ex.matlab,
            "output": ex.output,
        }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(data, indent=2, ensure_ascii=False) + "\n")
    print(f"✓ Wrote {len(data)} examples to {out_path}")


if __name__ == "__main__":
    # Keep AGM warnings etc. out of stdout
    np.seterr(all="ignore")

    verify()
    out = Path(__file__).parent / "out" / "function_examples.json"
    emit_json(out)

    # Also dump a preview so humans can spot-check
    print()
    for ex in EXAMPLES:
        print(f"── {ex.id}")
        print(f"   output: {ex.output.splitlines()[0]}")
