"""`clinval held-out`: held_out of an existing run recomputed from a VCF reference."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from typer.testing import CliRunner

from insilico_clinval.cli import app
from insilico_clinval.variants.held_out import read_held_out
from tests.conftest import EREPO_SUITE


def test_held_out_command(tmp_path: Path) -> None:
    runner = CliRunner()
    out = tmp_path / "run"
    assert (
        runner.invoke(app, ["select", "--suite", str(EREPO_SUITE), "--out", str(out)]).exit_code
        == 0
    )
    variants = json.loads((out / "variants.selected.json").read_text())["variants"]
    erepo = [v for v in variants if v["clingen"] is not None]
    known = erepo[0]
    vcf = tmp_path / "benchmark.vcf"
    vcf.write_text(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=10000>\n##contig=<ID=chr2,length=10000>\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        f"{known['chrom']}\t{known['pos']}\t.\t{known['ref']}\t{known['alt']}\t.\tPASS\t.\n"
    )
    args = ["held-out", "--suite", str(EREPO_SUITE), "--out", str(out), "--reference", str(vcf),
            "--release", "v2.2.0", "--approved-after", "2025-01-01"]  # fmt: skip
    result = runner.invoke(app, args)
    assert result.exit_code == 0, result.output
    assert f"{len(erepo) - 1} of {len(erepo)} eRepo variants held out" in result.output
    held = read_held_out(out / "held_out.tsv")
    assert held[known["variant_id"]] is False
    assert sum(held.values()) == len(erepo) - 1
    with (out / "held_out.tsv").open() as fh:
        rows = {r["variant_id"]: r for r in csv.DictReader(fh, delimiter="\t")}
    # every fixture record was approved after the cutoff, so only the known allele disagrees
    assert rows[known["variant_id"]]["agrees_with_date"] == "false"
    assert {r["agrees_with_date"] for k, r in rows.items() if k != known["variant_id"]} == {"true"}
    manifest = json.loads((out / "run_manifest.json").read_text())
    assert any(w.startswith("held-out: 1 of") for w in manifest["warnings"])
    assert {"held_out"} <= {o["role"] for o in manifest["outputs"]}
    assert any(d["name"] == "ClinGen eRepo (held-out reference)" for d in manifest["databases"])
    assert read_held_out(tmp_path / "missing.tsv") == {}
