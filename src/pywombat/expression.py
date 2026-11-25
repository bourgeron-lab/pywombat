"""Expression-based filter parsing: YAML expressions → Polars expressions."""

import re

import polars as pl


def parse_impact_filter_expression(expression: str, df: pl.DataFrame) -> pl.Expr:
    """Parse a filter expression string into a Polars expression.

    Supported operators:
        - Comparison: =, !=, <=, >=, <, >
        - Logical: &, |, ()
        - Special: contains, is_empty, not_empty, is_snv, is_indel

    Column names can contain dots (e.g., cadd_v1.7, faf95.joint).
    """
    expr_str = expression.strip()

    tokens = re.split(r"(\s*[&|]\s*|\(|\))", expr_str)
    tokens = [t.strip() for t in tokens if t.strip()]

    if tokens and tokens[-1] in ("&", "|"):
        raise ValueError(
            f"Expression ends with a dangling '{tokens[-1]}': {expression!r}"
        )
    for i, tok in enumerate(tokens[:-1]):
        if tok in ("&", "|") and tokens[i + 1] in ("&", "|", ")"):
            raise ValueError(
                f"Operator '{tok}' is not followed by a condition: {expression!r}"
            )
        if tok == "(" and tokens[i + 1] in ("&", "|"):
            raise ValueError(
                f"Operator '{tokens[i + 1]}' right after '(': {expression!r}"
            )

    def parse_condition(condition: str) -> pl.Expr:
        """Parse a single condition into a Polars expression."""
        condition = condition.strip()

        if condition.lower() == "is_snv":
            if "REF" not in df.columns or "ALT" not in df.columns:
                raise ValueError("is_snv requires REF and ALT columns in dataframe")

            nucleotides = ["a", "c", "g", "t", "A", "C", "G", "T"]
            return (
                (pl.col("REF").str.len_chars() == 1)
                & (pl.col("ALT").str.len_chars() == 1)
                & (pl.col("REF").is_in(nucleotides))
                & (pl.col("ALT").is_in(nucleotides))
            )

        if condition.lower() == "is_indel":
            if "REF" not in df.columns or "ALT" not in df.columns:
                raise ValueError("is_indel requires REF and ALT columns in dataframe")

            nucleotides = ["a", "c", "g", "t", "A", "C", "G", "T"]
            is_snv = (
                (pl.col("REF").str.len_chars() == 1)
                & (pl.col("ALT").str.len_chars() == 1)
                & (pl.col("REF").is_in(nucleotides))
                & (pl.col("ALT").is_in(nucleotides))
            )
            return ~is_snv

        is_empty_match = re.match(r"^([\w.]+)\s+is_empty$", condition, re.IGNORECASE)
        if is_empty_match:
            col_name = is_empty_match.group(1)
            if col_name not in df.columns:
                raise ValueError(f"Column '{col_name}' not found in dataframe")
            col_expr = pl.col(col_name)
            return col_expr.is_null() | (col_expr == "") | (col_expr == ".")

        not_empty_match = re.match(r"^([\w.]+)\s+not_empty$", condition, re.IGNORECASE)
        if not_empty_match:
            col_name = not_empty_match.group(1)
            if col_name not in df.columns:
                raise ValueError(f"Column '{col_name}' not found in dataframe")
            col_expr = pl.col(col_name)
            return ~(col_expr.is_null() | (col_expr == "") | (col_expr == "."))

        contains_match = re.match(r"^([\w.]+)\s+contains\s+(.+)$", condition, re.IGNORECASE)
        if contains_match:
            col_name = contains_match.group(1)
            value = contains_match.group(2).strip().strip("'\"")
            if col_name not in df.columns:
                raise ValueError(f"Column '{col_name}' not found in dataframe")
            col_expr = pl.col(col_name).str.to_lowercase()
            return col_expr.str.contains(value.lower())

        for op in ["<=", ">=", "!=", "=", "<", ">"]:
            if op in condition:
                parts = condition.split(op, 1)
                if len(parts) == 2:
                    col_name = parts[0].strip()
                    value = parts[1].strip()

                    if col_name not in df.columns:
                        raise ValueError(f"Column '{col_name}' not found in dataframe")

                    if value.upper() == "NULL":
                        col_expr = pl.col(col_name)
                        if op == "=":
                            return col_expr.is_null()
                        elif op == "!=":
                            return ~col_expr.is_null()
                        else:
                            raise ValueError(
                                f"Operator '{op}' not supported for null comparison, use = or !="
                            )

                    if value.upper() == "NAN":
                        col_expr = pl.col(col_name).cast(pl.Float64, strict=False)
                        if op == "=":
                            return col_expr.is_nan()
                        elif op == "!=":
                            return ~col_expr.is_nan()
                        else:
                            raise ValueError(
                                f"Operator '{op}' not supported for NaN comparison, use = or !="
                            )

                    try:
                        value_num = float(value)
                        col_expr = pl.col(col_name).cast(pl.Float64, strict=False)

                        if op == "=":
                            return col_expr == value_num
                        elif op == "!=":
                            return col_expr != value_num
                        elif op == "<=":
                            return col_expr <= value_num
                        elif op == ">=":
                            return col_expr >= value_num
                        elif op == "<":
                            return col_expr < value_num
                        elif op == ">":
                            return col_expr > value_num
                    except ValueError:
                        value = value.strip("'\"")
                        col_expr = pl.col(col_name).str.to_lowercase()
                        value_lower = value.lower()

                        if op == "=":
                            return col_expr == value_lower
                        elif op == "!=":
                            return col_expr != value_lower
                        else:
                            raise ValueError(
                                f"Operator '{op}' not supported for string comparison"
                            )
                break

        raise ValueError(f"Could not parse condition: {condition}")

    def build_expression(tokens: list, idx: int = 0) -> tuple[pl.Expr, int]:
        """Recursively build expression from tokens."""
        if idx >= len(tokens):
            return None, idx

        result = None
        i = idx

        while i < len(tokens):
            token = tokens[i]

            if token == "(":
                sub_expr, new_i = build_expression(tokens, i + 1)
                if result is None:
                    result = sub_expr
                else:
                    result = result & sub_expr
                i = new_i - 1
            elif token == ")":
                return result, i + 1
            elif token == "&":
                next_expr, new_i = build_expression(tokens, i + 1)
                if next_expr is not None:
                    result = result & next_expr if result is not None else next_expr
                return result, new_i
            elif token == "|":
                next_expr, new_i = build_expression(tokens, i + 1)
                if next_expr is not None:
                    result = result | next_expr if result is not None else next_expr
                return result, new_i
            else:
                cond_expr = parse_condition(token)
                if result is None:
                    result = cond_expr
                else:
                    result = result & cond_expr

            i += 1

        return result, i

    expr, _ = build_expression(tokens)
    return expr
