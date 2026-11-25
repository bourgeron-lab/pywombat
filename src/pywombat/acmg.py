"""ACMG/AMP variant-level evidence annotation.

Computes evidence criteria that need only the variant and the family:
frequency, PP3/BP4 (missense predictor), splice PP3/BP4, BP7,
PM4/BP3 candidates. Each criterion becomes an acmg_* column.
"""

import importlib.resources
import re
from pathlib import Path
from typing import Optional

import polars as pl
import yaml


_POINTS = {"supporting": 1, "moderate": 2, "strong": 4, "very_strong": 8}

_SPLICE_CANONICAL = {"splice_donor_variant", "splice_acceptor_variant"}

_MISSENSE_CONSEQUENCES = {"missense_variant"}

_SYNONYMOUS_CONSEQUENCES = {"synonymous_variant"}

_INFRAME_CONSEQUENCES = {"inframe_insertion", "inframe_deletion", "stop_lost"}

_INTRONIC_CONSEQUENCES = {"intron_variant"}


def load_ruleset(path_or_name: str | Path) -> dict:
    """Load a ruleset YAML by path or built-in name.

    If *path_or_name* has no path separator and no `.yml`/`.yaml` extension,
    it is looked up in the package's ``rulesets/`` directory.
    """
    p = Path(path_or_name)
    if p.exists():
        with open(p) as f:
            return yaml.safe_load(f)

    stem = str(path_or_name)
    if "/" not in stem and "\\" not in stem:
        if not stem.endswith((".yml", ".yaml")):
            stem += ".yml"
        pkg_dir = Path(__file__).parent / "rulesets"
        candidate = pkg_dir / stem
        if candidate.exists():
            with open(candidate) as f:
                return yaml.safe_load(f)

    raise FileNotFoundError(f"Ruleset not found: {path_or_name}")


def _get_version() -> str:
    try:
        from importlib.metadata import version
        return version("pywombat")
    except Exception:
        return "dev"


def _missing_columns(df: pl.DataFrame, ruleset: dict) -> list[str]:
    """Return annotation columns the ruleset expects but the DataFrame lacks."""
    missing = []
    cols = set(df.columns)

    freq = ruleset.get("frequency", {})
    if freq.get("column") and freq["column"] not in cols:
        missing.append(freq["column"])
    if freq.get("ac_column") and freq["ac_column"] not in cols:
        missing.append(freq["ac_column"])
    if freq.get("an_column") and freq["an_column"] not in cols:
        missing.append(freq["an_column"])

    pred = ruleset.get("predictor", {})
    if pred.get("column") and pred["column"] not in cols:
        missing.append(pred["name"] or pred["column"])

    spl = ruleset.get("splice", {})
    if spl.get("column") and spl["column"] not in cols:
        missing.append(spl["column"])

    if "VEP_Consequence" not in cols:
        missing.append("VEP_Consequence")

    if "LCR" not in cols:
        missing.append("LCR")

    return missing


def _freq_status(df: pl.DataFrame, ruleset: dict) -> pl.Series:
    """Classify each row as absent, absent_low_coverage, rare, or BA1."""
    freq_cfg = ruleset.get("frequency", {})
    freq_col = freq_cfg.get("column", "fafmax_faf95_max_joint")
    ba1_threshold = freq_cfg.get("ba1", 0.05)
    ac_col = freq_cfg.get("ac_column", "AC_joint")
    an_col = freq_cfg.get("an_column", "AN_joint")
    an_threshold = freq_cfg.get("absent_needs_an", 20000)

    n = df.height
    cols = set(df.columns)

    if freq_col not in cols:
        return pl.Series("acmg_freq", [None] * n, dtype=pl.Utf8)

    freq_val = df[freq_col].cast(pl.Float64, strict=False)

    has_ac = ac_col in cols
    has_an = an_col in cols

    result = []
    for i in range(n):
        fv = freq_val[i]

        if fv is not None and fv > ba1_threshold:
            result.append("BA1")
        elif fv is None or fv == 0.0:
            if has_ac and has_an:
                ac_raw = df[ac_col][i]
                an_raw = df[an_col][i]
                try:
                    ac_v = float(str(ac_raw)) if ac_raw is not None else None
                except (ValueError, TypeError):
                    ac_v = None
                try:
                    an_v = float(str(an_raw)) if an_raw is not None else None
                except (ValueError, TypeError):
                    an_v = None

                if ac_v is not None and ac_v == 0:
                    if an_v is not None and an_v >= an_threshold:
                        result.append("absent")
                    else:
                        result.append("absent_low_coverage")
                elif ac_v is None:
                    result.append("absent_low_coverage")
                else:
                    result.append("rare")
            else:
                result.append("absent")
        else:
            result.append("rare")

    return pl.Series("acmg_freq", result, dtype=pl.Utf8)


def _pp3_bp4_missense(df: pl.DataFrame, ruleset: dict) -> pl.Series:
    """PP3/BP4 from the missense predictor (REVEL by default)."""
    pred_cfg = ruleset.get("predictor", {})
    pred_col = pred_cfg.get("column", "REVEL")
    pred_name = pred_cfg.get("name", pred_col)
    pp3_cuts = pred_cfg.get("pp3", {})
    bp4_cuts = pred_cfg.get("bp4", {})

    n = df.height
    cols = set(df.columns)

    if pred_col not in cols or "VEP_Consequence" not in cols:
        return pl.Series("acmg_pp3_bp4", [None] * n, dtype=pl.Utf8)

    score = df[pred_col].cast(pl.Float64, strict=False)
    consequence = df["VEP_Consequence"].cast(pl.Utf8)

    pp3_strong = pp3_cuts.get("strong", 0.932)
    pp3_moderate = pp3_cuts.get("moderate", 0.773)
    pp3_supporting = pp3_cuts.get("supporting", 0.644)
    bp4_supporting = bp4_cuts.get("supporting", 0.290)
    bp4_moderate = bp4_cuts.get("moderate", 0.183)
    bp4_strong = bp4_cuts.get("strong", 0.016)
    bp4_very_strong = bp4_cuts.get("very_strong", 0.003)

    result = []
    for i in range(n):
        csq = consequence[i]
        sv = score[i]

        if csq is None or not _has_consequence(csq, _MISSENSE_CONSEQUENCES):
            result.append(None)
            continue

        if sv is None:
            result.append(None)
            continue

        if sv >= pp3_strong:
            result.append(f"PP3_Strong:{pred_name}={sv}")
        elif sv >= pp3_moderate:
            result.append(f"PP3_Moderate:{pred_name}={sv}")
        elif sv >= pp3_supporting:
            result.append(f"PP3_Supporting:{pred_name}={sv}")
        elif sv <= bp4_very_strong:
            result.append(f"BP4_VeryStrong:{pred_name}={sv}")
        elif sv <= bp4_strong:
            result.append(f"BP4_Strong:{pred_name}={sv}")
        elif sv <= bp4_moderate:
            result.append(f"BP4_Moderate:{pred_name}={sv}")
        elif sv <= bp4_supporting:
            result.append(f"BP4_Supporting:{pred_name}={sv}")
        else:
            result.append(None)

    return pl.Series("acmg_pp3_bp4", result, dtype=pl.Utf8)


def _splice_evidence(df: pl.DataFrame, ruleset: dict) -> pl.Series:
    """Splice PP3/BP4 from SpliceAI. Not applied to canonical ±1/2 sites."""
    spl_cfg = ruleset.get("splice", {})
    spl_col = spl_cfg.get("column", "SpliceAI_maxDS")
    pp3_min = spl_cfg.get("pp3_min", 0.2)
    bp4_max = spl_cfg.get("bp4_max", 0.1)

    n = df.height
    cols = set(df.columns)

    if spl_col not in cols:
        return pl.Series("acmg_splice", [None] * n, dtype=pl.Utf8)

    score = df[spl_col].cast(pl.Float64, strict=False)

    has_csq = "VEP_Consequence" in cols
    consequence = df["VEP_Consequence"].cast(pl.Utf8) if has_csq else None

    result = []
    for i in range(n):
        sv = score[i]
        if sv is None:
            result.append(None)
            continue

        if has_csq and consequence[i] is not None:
            if _has_consequence(consequence[i], _SPLICE_CANONICAL):
                result.append(None)
                continue

        if sv >= pp3_min:
            result.append(f"PP3_Supporting:SpliceAI={sv}")
        elif sv <= bp4_max:
            result.append(f"BP4_Supporting:SpliceAI={sv}")
        else:
            result.append(None)

    return pl.Series("acmg_splice", result, dtype=pl.Utf8)


def _bp7(df: pl.DataFrame, ruleset: dict) -> pl.Series:
    """BP7: synonymous or deep intronic with low SpliceAI."""
    spl_cfg = ruleset.get("splice", {})
    spl_col = spl_cfg.get("column", "SpliceAI_maxDS")
    bp4_max = spl_cfg.get("bp4_max", 0.1)
    offsets = spl_cfg.get("bp7_intron_offsets", {})
    donor_offset = offsets.get("donor", 7)
    acceptor_offset = offsets.get("acceptor", -21)

    n = df.height
    cols = set(df.columns)
    has_spl = spl_col in cols
    has_csq = "VEP_Consequence" in cols
    has_hgvsc = "VEP_HGVSc" in cols

    if not has_csq:
        return pl.Series("acmg_bp7", [None] * n, dtype=pl.Utf8)

    consequence = df["VEP_Consequence"].cast(pl.Utf8)
    splice_score = df[spl_col].cast(pl.Float64, strict=False) if has_spl else None
    hgvsc = df["VEP_HGVSc"].cast(pl.Utf8) if has_hgvsc else None

    result = []
    for i in range(n):
        csq = consequence[i]
        if csq is None:
            result.append(None)
            continue

        is_syn = _has_consequence(csq, _SYNONYMOUS_CONSEQUENCES)
        is_intronic = _has_consequence(csq, _INTRONIC_CONSEQUENCES)

        if not is_syn and not is_intronic:
            result.append(None)
            continue

        if has_spl:
            sv = splice_score[i] if splice_score is not None else None
            if sv is not None and sv > bp4_max:
                result.append(None)
                continue
        else:
            result.append(None)
            continue

        if is_syn:
            result.append("BP7")
            continue

        if is_intronic and has_hgvsc:
            offset = _parse_intron_offset(hgvsc[i])
            if offset is not None:
                if offset > 0 and offset >= donor_offset:
                    result.append("BP7")
                elif offset < 0 and offset <= acceptor_offset:
                    result.append("BP7")
                else:
                    result.append(None)
            else:
                result.append(None)
        else:
            result.append(None)

    return pl.Series("acmg_bp7", result, dtype=pl.Utf8)


def _pm4_bp3(df: pl.DataFrame) -> pl.Series:
    """PM4/BP3 candidates from in-frame indels and stop_lost."""
    n = df.height
    cols = set(df.columns)

    if "VEP_Consequence" not in cols:
        return pl.Series("acmg_pm4_bp3", [None] * n, dtype=pl.Utf8)

    consequence = df["VEP_Consequence"].cast(pl.Utf8)
    has_lcr = "LCR" in cols
    lcr = df["LCR"] if has_lcr else None

    result = []
    for i in range(n):
        csq = consequence[i]
        if csq is None:
            result.append(None)
            continue

        if _has_consequence(csq, _INFRAME_CONSEQUENCES):
            if has_lcr and lcr is not None:
                lcr_val = lcr[i]
                if _is_truthy(lcr_val):
                    result.append("BP3_candidate")
                    continue
            result.append("PM4_candidate")
        else:
            result.append(None)

    return pl.Series("acmg_pm4_bp3", result, dtype=pl.Utf8)


def _has_consequence(csq_string: str, target_set: set) -> bool:
    """Check if any consequence in a &-separated list is in the target set."""
    for csq in csq_string.split("&"):
        if csq.strip().lower() in {c.lower() for c in target_set}:
            return True
    return False


def _is_truthy(val) -> bool:
    """Check if a value is truthy (not null, not empty, not '.', not False)."""
    if val is None:
        return False
    if isinstance(val, bool):
        return val
    s = str(val).strip().lower()
    return s not in ("", ".", "false", "0", "none", "null")


def _parse_intron_offset(hgvsc: Optional[str]) -> Optional[int]:
    """Extract the intronic offset from a VEP HGVSc notation.

    E.g. 'ENST00000454634.6:c.123+7A>G' → 7
         'ENST00000454634.6:c.124-21A>G' → -21
    """
    if hgvsc is None:
        return None
    m = re.search(r"[+-]\d+", hgvsc)
    if m:
        try:
            return int(m.group())
        except ValueError:
            return None
    return None


def _resolve_criteria(row: dict, points_map: dict) -> tuple[list[str], int]:
    """Combine individual acmg_* columns into criteria list and point sum.

    PP3 fires at most once (stronger of missense and splice).
    BP4 fires at most once (stronger of missense and splice).
    """
    criteria = []
    total = 0

    # Frequency: not scored here (PM2/BS1 need mode of inheritance)
    # but BA1 is standalone
    freq = row.get("acmg_freq")
    if freq == "BA1":
        criteria.append("BA1")
        total -= points_map.get("very_strong", 8)

    # PP3/BP4: resolve once-rule across missense and splice
    pp3_bp4 = row.get("acmg_pp3_bp4")
    splice = row.get("acmg_splice")

    missense_criterion = _extract_criterion(pp3_bp4)
    splice_criterion = _extract_criterion(splice)

    pp3_entries = []
    bp4_entries = []

    for crit in [missense_criterion, splice_criterion]:
        if crit is None:
            continue
        name, strength = crit
        if name == "PP3":
            pp3_entries.append((strength, crit))
        elif name == "BP4":
            bp4_entries.append((strength, crit))

    if pp3_entries:
        best = max(pp3_entries, key=lambda x: points_map.get(x[0], 0))
        strength = best[0]
        criteria.append(f"PP3_{_capitalize_strength(strength)}")
        total += points_map.get(strength, 0)

    if bp4_entries:
        best = max(bp4_entries, key=lambda x: points_map.get(x[0], 0))
        strength = best[0]
        criteria.append(f"BP4_{_capitalize_strength(strength)}")
        total -= points_map.get(strength, 0)

    # BP7
    bp7 = row.get("acmg_bp7")
    if bp7 == "BP7":
        criteria.append("BP7")
        total -= points_map.get("supporting", 1)

    # PM4/BP3
    pm4_bp3 = row.get("acmg_pm4_bp3")
    if pm4_bp3 == "PM4_candidate":
        criteria.append("PM4_candidate")
    elif pm4_bp3 == "BP3_candidate":
        criteria.append("BP3_candidate")

    return criteria, total


def _extract_criterion(value: Optional[str]) -> Optional[tuple[str, str]]:
    """Parse 'PP3_Moderate:REVEL=0.81' → ('PP3', 'moderate')."""
    if value is None:
        return None
    m = re.match(r"(PP3|BP4)_(\w+)", value)
    if m:
        return (m.group(1), m.group(2).lower())
    return None


def _capitalize_strength(s: str) -> str:
    if s == "very_strong":
        return "VeryStrong"
    return s.capitalize()


def annotate_acmg(
    df: pl.DataFrame,
    ruleset: dict,
) -> pl.DataFrame:
    """Add ACMG evidence columns to a filtered DataFrame.

    Args:
        df: Filtered variant DataFrame (one row per variant×sample)
        ruleset: Parsed ruleset dict (from load_ruleset)

    Returns:
        DataFrame with acmg_* columns appended
    """
    ruleset_name = ruleset.get("name", "unknown")
    version = _get_version()
    points_map = ruleset.get("points", _POINTS)

    missing = _missing_columns(df, ruleset)

    # Compute each criterion
    freq_series = _freq_status(df, ruleset)
    pp3_bp4_series = _pp3_bp4_missense(df, ruleset)
    splice_series = _splice_evidence(df, ruleset)
    bp7_series = _bp7(df, ruleset)
    pm4_bp3_series = _pm4_bp3(df)

    df = df.with_columns([
        freq_series,
        pp3_bp4_series,
        splice_series,
        bp7_series,
        pm4_bp3_series,
    ])

    # Resolve criteria and points per row
    criteria_list = []
    points_list = []

    for i in range(df.height):
        row = {
            "acmg_freq": df["acmg_freq"][i],
            "acmg_pp3_bp4": df["acmg_pp3_bp4"][i],
            "acmg_splice": df["acmg_splice"][i],
            "acmg_bp7": df["acmg_bp7"][i],
            "acmg_pm4_bp3": df["acmg_pm4_bp3"][i],
        }
        crit, pts = _resolve_criteria(row, points_map)
        criteria_list.append(";".join(crit) if crit else None)
        points_list.append(pts if pts != 0 else None)

    missing_str = ",".join(missing) if missing else None

    df = df.with_columns([
        pl.lit(f"{ruleset_name}/pywombat-{version}").alias("acmg_ruleset"),
        pl.Series("acmg_criteria", criteria_list, dtype=pl.Utf8),
        pl.Series("acmg_points_auto", points_list, dtype=pl.Int32),
        pl.lit(missing_str).alias("acmg_missing"),
    ])

    return df
