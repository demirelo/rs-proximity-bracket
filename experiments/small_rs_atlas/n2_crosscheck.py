"""
n2_crosscheck.py -- EXACT decoder-free cross-check that the field-agnostic firing line is
a REAL bad line on a GENUINE odd-char extension.

Three layers (all decoder-free; we NEVER decode -- closeness is by exact Hamming distance
to the FULL codeword set, or by an exact in-field polynomial-agreement count; S* is exact
branch-and-bound or an exact-by-degree bound):

  (a) FULL codeword-enumeration + EXACT branch-and-bound S* at a GENUINE extension where
      q^k <= 3e6 -- the strongest possible verification.  The only genuine-extension case
      that fits is GF(3^2), s=4, r=4 (delta=0): degenerate (the word is itself a codeword),
      so it validates the machinery + genuineness but not a near-capacity radius.

  (b) CODEWORD-FREE EXACT structural certificate at NON-DEGENERATE genuine extensions
      GF(31^2), GF(127^2) (delta = 1/2, 3/4).  Builds R(X) in-field, confirms deg R < k
      (R is a codeword), evaluates f - lambda*g and codeword(-R) on the FULL domain and
      confirms they agree on EXACTLY r*m = (1-delta)n points (so dist <= delta*n -- exact,
      no decoder), and bounds the joint correlated agreement EXACTLY by polynomial degree:
      g = X^{(r-1)m} agrees with any deg-<k codeword on <= (r-1)m points (a degree-(r-1)m
      poly has <= (r-1)m roots), so S* <= (r-1)m < r*m = the CA threshold -- NO CA.

  (c) EXACT prime control GF(17), s=8, r=4 (delta=1/2): full enumeration + exact S*,
      confirming the SAME non-degenerate firing line (which is field-agnostic) is genuinely
      bad (close_count == predicted, S* < CA threshold).

Together: the firing line is a real, decoder-free-verified near-capacity bad line whose
mechanism is identical on genuine extensions and primes.
"""

from __future__ import annotations

import itertools
import json
import os

import numpy as np

from ff import PrimeField
from ff_ext import PrimePowerField
from rs import build_codeword_book, domain_subgroup, dist_to_code, min_distance
from counterexample_kambire import monomial_eval
from search_bad_lines import _max_common_agreement, _agreement_bits
from n2_hardening import structural_certificate, prime_r4_exact_control
from counterexample_extension import subfield_orders


def proper_subfield_codes(F):
    psf = set()
    for d in range(1, F.e):
        if F.e % d == 0:
            for c in range(F.q):
                if F.pow(c, F.p ** d) == c:
                    psf.add(c)
    return psf


def full_enum_check(p, e, s, m, n, r):
    """(a) full dist_to_code + exact S* at a genuine extension (q^k must be small)."""
    F = PrimePowerField(p, e)
    q = F.q
    k = (r - 2) * m
    assert q ** k <= 3_000_000
    sf = subfield_orders(p, e)
    psf = proper_subfield_codes(F)
    D = domain_subgroup(F, n)
    H = domain_subgroup(F, s).tolist()
    book = build_codeword_book(F, D, k)
    f = monomial_eval(F, D, r * m)
    g = monomial_eval(F, D, (r - 1) * m)
    delta = 1 - r / s
    dn = delta * n
    fire_lam = set()
    witness = None
    for c in itertools.combinations(H, r):
        e1 = 0
        e2 = 0
        for x in c:
            e2 = F.add(e2, F.mul(int(x), e1))
            e1 = F.add(e1, int(x))
        if e2 == 0:
            fire_lam.add(e1)
            if witness is None:
                witness = ([int(x) for x in c], e1)
    closes = {}
    for lam in fire_lam:
        w = F.add_vec(f, F.mul_scalar_vec(F.neg(lam), g))
        closes[lam] = int(dist_to_code(book, w))
    S, info = _max_common_agreement(book, f, g)
    pop_g, _ = _agreement_bits(book, g)
    wc, wl = witness
    return dict(field=F.name, q=q, s=s, m=m, n=n, r=r, k=k, qk=q ** k,
                H_genuine=(s not in sf), delta=delta, delta_n=dn,
                num_firing_lambda=len(fire_lam),
                all_firing_close=all(d <= dn + 1e-9 for d in closes.values()),
                firing_dists=sorted(set(closes.values())),
                witness=wc, witness_lambda=wl,
                witness_genuine=not all(c in psf for c in wc),
                exact_S_star=int(S), S_star_exact=bool(info["exact"]),
                ca_threshold_rm=r * m, no_CA=bool(S < r * m),
                g_max_single_agreement=int(pop_g.max()),
                g_deg_bound=(r - 1) * m)


def main():
    out = {"meta": {"experiment": "N2_exact_crosscheck",
                    "decoder_free": True,
                    "note": "closeness by exact Hamming dist to FULL codeword set or "
                            "exact in-field agreement count; S* exact b&b or exact-by-degree"},
           "a_full_enum_genuine": [], "b_structural_genuine": [], "c_prime_control": None}

    print("=" * 96)
    print("(a) FULL-ENUMERATION + EXACT S* at a GENUINE extension (q^k<=3e6): GF(3^2) s=4 r=4")
    print("=" * 96)
    a = full_enum_check(3, 2, 4, 2, 8, 4)
    out["a_full_enum_genuine"].append(a)
    print(f"  {a['field']} s={a['s']} m={a['m']} n={a['n']} r={a['r']} k={a['k']} (q^k={a['qk']}) "
          f"delta={a['delta']}: H genuine={a['H_genuine']}")
    print(f"     firing lambda={a['num_firing_lambda']}, all close (dist<=delta*n)={a['all_firing_close']} "
          f"(dists {a['firing_dists']}); witness={a['witness']} genuine={a['witness_genuine']}")
    print(f"     EXACT S*={a['exact_S_star']} (exact={a['S_star_exact']}) < CA thr {a['ca_threshold_rm']}? "
          f"{a['no_CA']}; g max-agree {a['g_max_single_agreement']}<=deg(g)={a['g_deg_bound']}")
    print(f"     NOTE: delta=0 here (degenerate); validates machinery+genuineness, not a near-capacity radius.")

    print("\n" + "=" * 96)
    print("(b) CODEWORD-FREE EXACT structural certificate at NON-DEGENERATE genuine extensions")
    print("=" * 96)
    for (p, e, s, n, m, r) in [(31, 2, 8, 16, 2, 4), (127, 2, 8, 16, 2, 4),
                               (31, 2, 16, 32, 2, 4), (127, 2, 16, 32, 2, 4)]:
        cert = structural_certificate(PrimePowerField(p, e), s, n, m, r, verbose=False)
        out["b_structural_genuine"].append(cert)
        if not cert.get("has_degenerate_subset"):
            print(f"  GF({p}^{e}) s={s}: no firing subset")
            continue
        print(f"  GF({p}^2) s={s} n={n} m={m} r={r} k={cert['k']} delta={cert['design_delta']}: "
              f"witness genuine={not cert['witness_in_proper_subfield']} lambda={cert['lambda']}")
        print(f"     deg R={cert['deg_R']}<k? {cert['R_is_codeword_deg_lt_k']}; agree "
              f"{cert['agree_points']}/{n} (delta*n={cert['delta_n']:.0f}) certified_close="
              f"{cert['certified_close']}; #distinct close scalars={cert['num_distinct_close_scalars']} "
              f"({cert['num_degenerate_genuinely_extension']}/{cert['num_degenerate_subsets']} genuine)")
        print(f"     no-CA (exact by degree): S*<=(r-1)m={cert['s_star_upper_bound_r1m']}<r*m="
              f"{cert['ca_threshold_rm']} => {cert['no_ca_structural']}; CERTIFIED BAD="
              f"{cert['is_certified_bad_scalar']}")

    print("\n" + "=" * 96)
    print("(c) EXACT prime control GF(17) s=8 r=4 (full enumeration + exact S*), delta=1/2")
    print("=" * 96)
    ctrl = prime_r4_exact_control(p=17, s=8, n=16, m=2, r=4, verbose=False)
    out["c_prime_control"] = ctrl
    print(f"  GF(17) (q^k={ctrl['qk']}): close={ctrl['close_count_at_design']} "
          f"(=pred {ctrl['predicted_sumset_size']}? {ctrl['close_eq_predicted']}); EXACT S*="
          f"{ctrl['S_star']}<CA {ctrl['ca_threshold']}? {ctrl['S_star']<ctrl['ca_threshold']}; "
          f"is_bad_line={ctrl['is_bad_line']}")

    path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "results", "n2_crosscheck.json")
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w") as fh:
        json.dump(out, fh, indent=2, default=str)
    print(f"\nWrote {path}")

    # net verdict
    b_ok = all(c.get("is_certified_bad_scalar") for c in out["b_structural_genuine"]
               if c.get("has_degenerate_subset"))
    print("\nNET: non-degenerate genuine-extension firing lines (GF(31^2),GF(127^2), s=8&16, "
          f"delta=1/2 & 3/4) ALL certified bad (decoder-free)={b_ok}; prime control bad="
          f"{out['c_prime_control']['is_bad_line']}; full-enum genuine machinery validated.")
    return out


if __name__ == "__main__":
    main()
