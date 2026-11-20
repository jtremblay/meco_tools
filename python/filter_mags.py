#!/usr/bin/env python3
"""Keep a subset of the MAGs of a ShotgunMG run, by DAS Tool score or CheckM2 quality.

Reads the MAG directory of a ShotgunMG run (<outdir>/binning/mags) and writes the same files,
restricted to the MAGs that pass every filter given, to a new directory:

  <prefix>.summary.tsv              <prefix>_checkm2_report.tsv
  <prefix>.bin2contig.tsv           <prefix>.bins_list.txt
  <prefix>_genes.link.tsv           <prefix>.feature_table.tsv
  <prefix>.feature_table_normalized.tsv
  bat/<prefix>.bin2classification.txt
  bat/<prefix>.bin2classification_with_names.tsv
  bat/<prefix>.ORF2LCA.txt          bins/*.fa
  filter_mags.log                   the command line and what was kept

The abundance of a MAG is the sum of the counts of its own contigs, so the filtered feature
tables are exact. alpha_div/ and beta_div/ describe the whole set of MAGs: they are not copied
and have to be computed again from the filtered feature tables.

DAS Tool keeps selecting bins below its score threshold and only stops writing them out. The
MAGs of a run at --dastool_score_threshold 0.3 that score at least 0.5 are therefore those of a
run at 0.5, with the same contigs.

Examples:
  filter_mags.py -i results/binning/mags -o mags_score0.5 --min_dastool_score 0.5
  filter_mags.py -i results/binning/mags -o mags_hq_mq --quality high,medium

Uses the Python standard library only.
"""

import argparse
import os
import shutil
import sys
from collections import Counter

QUALITIES = ("high", "medium", "low")

# Tables to subset: (path relative to the MAG directory, column holding the MAG name,
# first line is a header, the MAG name carries the bin fasta suffix). Lines starting with
# '#' are always kept as header lines.
TABLES = [
    ("{prefix}_checkm2_report.tsv", 0, True, False),
    ("{prefix}.bin2contig.tsv", 0, False, False),
    ("{prefix}.bins_list.txt", 0, False, False),
    ("{prefix}_genes.link.tsv", 0, False, False),
    ("{prefix}.feature_table.tsv", 0, False, False),
    ("{prefix}.feature_table_normalized.tsv", 0, False, False),
    ("bat/{prefix}.bin2classification.txt", 0, False, True),
    ("bat/{prefix}.bin2classification_with_names.tsv", 0, False, True),
    ("bat/{prefix}.ORF2LCA.txt", 1, False, True),
]
NOT_COPIED = ("alpha_div", "beta_div")


def parse_args(arguments):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-i", "--indir", required=True, help="MAG directory of a ShotgunMG run (<outdir>/binning/mags)")
    parser.add_argument("-o", "--outdir", required=True, help="New directory for the filtered MAGs (must not exist or be empty)")
    parser.add_argument("--min_dastool_score", type=float, help="Keep MAGs whose DAS Tool score is at least this value")
    parser.add_argument("--min_completeness", type=float, help="Keep MAGs whose CheckM2 completeness (%%) is at least this value")
    parser.add_argument("--max_contamination", type=float, help="Keep MAGs whose CheckM2 contamination (%%) is at most this value")
    parser.add_argument("--quality", help="Keep MAGs of these quality categories, comma-separated among " + ", ".join(QUALITIES))
    parser.add_argument("--prefix", default="mags", help="Prefix of the MAG files [mags]")
    parser.add_argument("--bin_suffix", default=".fa", help="Suffix of the MAG fasta files in bins/ [.fa]")
    parser.add_argument("--symlink_bins", action="store_true", help="Link the MAG fasta files instead of copying them")
    args = parser.parse_args(arguments)

    if all(v is None for v in (args.min_dastool_score, args.min_completeness, args.max_contamination, args.quality)):
        parser.error("give at least one filter (--min_dastool_score, --min_completeness, --max_contamination, --quality)")
    if args.quality is not None:
        args.quality = {q.strip().lower() for q in args.quality.split(",") if q.strip()}
        unknown = args.quality - set(QUALITIES)
        if unknown or not args.quality:
            parser.error("--quality takes a comma-separated list among {0}".format(", ".join(QUALITIES)))
    if not os.path.isdir(args.indir):
        parser.error("input directory not found: {0}".format(args.indir))
    if os.path.realpath(args.indir) == os.path.realpath(args.outdir):
        parser.error("the output directory must differ from the input directory")
    if os.path.isdir(args.outdir) and os.listdir(args.outdir):
        parser.error("the output directory is not empty: {0}".format(args.outdir))
    return args


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def select_mags(summary_path, args):
    """Return (all MAGs, kept MAGs, number of MAGs dropped for a missing value, header, rows)."""
    with open(summary_path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not header or header[0] != "Bin_Id":
        sys.exit("ERROR: {0} does not start with a Bin_Id column".format(summary_path))

    # (column, test) for each filter given
    filters = []
    if args.min_dastool_score is not None:
        filters.append(("DASTool_score", lambda v: v >= args.min_dastool_score))
    if args.min_completeness is not None:
        filters.append(("Completeness", lambda v: v >= args.min_completeness))
    if args.max_contamination is not None:
        filters.append(("Contamination", lambda v: v <= args.max_contamination))
    missing = [col for col, _ in filters if col not in header]
    if args.quality is not None and "Quality" not in header:
        missing.append("Quality")
    if missing:
        sys.exit("ERROR: column(s) missing from {0}: {1}".format(summary_path, ", ".join(missing)))

    kept, no_value = [], 0
    for row in rows:
        values = dict(zip(header, row))
        ok = True
        for col, test in filters:
            value = as_float(values.get(col))
            if value is None:
                no_value += 1
                ok = False
                break
            if not test(value):
                ok = False
                break
        if ok and args.quality is not None and values.get("Quality", "").lower() not in args.quality:
            ok = False
        if ok:
            kept.append(row[0])
    return [row[0] for row in rows], kept, no_value, header, rows


def subset_table(src, dst, keep, key_col, has_header, suffix):
    """Copy the header lines and the lines of the kept MAGs; return (kept, total) data lines."""
    kept = total = 0
    with open(src) as fin, open(dst, "w") as fout:
        for number, line in enumerate(fin):
            if (has_header and number == 0) or line.startswith("#"):
                fout.write(line)
                continue
            if not line.strip():
                continue
            total += 1
            fields = line.rstrip("\n").split("\t")
            key = fields[key_col] if key_col < len(fields) else ""
            if suffix and key.endswith(suffix):
                key = key[: -len(suffix)]
            if key in keep:
                fout.write(line)
                kept += 1
    return kept, total


def main(arguments):
    args = parse_args(arguments)
    indir, outdir, prefix = args.indir, args.outdir, args.prefix

    summary_name = "{0}.summary.tsv".format(prefix)
    summary_path = os.path.join(indir, summary_name)
    if not os.path.isfile(summary_path):
        sys.exit("ERROR: {0} not found; is {1} the binning/mags directory of a ShotgunMG run?".format(summary_path, indir))
    all_mags, kept_mags, no_value, header, rows = select_mags(summary_path, args)
    keep = set(kept_mags)

    os.makedirs(outdir, exist_ok=True)
    report = ["filter_mags.py " + " ".join(arguments), ""]
    report.append("MAGs kept: {0} of {1}".format(len(kept_mags), len(all_mags)))
    if no_value:
        report.append("MAGs dropped because a filtered value is missing (NA): {0}".format(no_value))
    if "Quality" in header:
        quality_col = header.index("Quality")
        counts = Counter(row[quality_col] for row in rows if row[0] in keep)
        report.append("Quality of the kept MAGs: " + ", ".join("{0} {1}".format(q, counts.get(q, 0)) for q in QUALITIES))
    report.append("")

    with open(summary_path) as fin, open(os.path.join(outdir, summary_name), "w") as fout:
        fout.write(fin.readline())
        for line in fin:
            if line.split("\t", 1)[0] in keep:
                fout.write(line)
    report.append("{0}: {1} MAGs".format(summary_name, len(kept_mags)))

    for template, key_col, has_header, strip_suffix in TABLES:
        name = template.format(prefix=prefix)
        src = os.path.join(indir, name)
        if not os.path.isfile(src):
            report.append("{0}: not found, skipped".format(name))
            continue
        dst = os.path.join(outdir, name)
        os.makedirs(os.path.dirname(dst), exist_ok=True)
        kept, total = subset_table(src, dst, keep, key_col, has_header, args.bin_suffix if strip_suffix else "")
        report.append("{0}: {1} of {2} lines".format(name, kept, total))

    bins_in = os.path.join(indir, "bins")
    if os.path.isdir(bins_in):
        bins_out = os.path.join(outdir, "bins")
        os.makedirs(bins_out, exist_ok=True)
        absent = []
        for mag in kept_mags:
            src = os.path.join(bins_in, mag + args.bin_suffix)
            if not os.path.isfile(src):
                absent.append(mag)
                continue
            dst = os.path.join(bins_out, mag + args.bin_suffix)
            if args.symlink_bins:
                os.symlink(os.path.realpath(src), dst)
            else:
                shutil.copy2(src, dst)
        report.append("bins/: {0} fasta files {1}".format(len(kept_mags) - len(absent), "linked" if args.symlink_bins else "copied"))
        if absent:
            report.append("WARNING: no fasta file in bins/ for {0} MAGs: {1}".format(len(absent), ", ".join(absent[:10])))
    else:
        report.append("bins/: not found, skipped")

    skipped = [d for d in NOT_COPIED if os.path.isdir(os.path.join(indir, d))]
    if skipped:
        report.append("")
        report.append("Not copied: {0} (computed on all the MAGs; compute them again from the filtered feature tables)".format(
            ", ".join(d + "/" for d in skipped)))

    text = "\n".join(report) + "\n"
    with open(os.path.join(outdir, "filter_mags.log"), "w") as out:
        out.write(text)
    sys.stdout.write(text)
    if not kept_mags:
        sys.stderr.write("WARNING: no MAG passed the filters\n")


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
