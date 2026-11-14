#!/usr/bin/env bash
# Smoke test for the meco_tools container image. Run inside the image:
#   docker run --rm meco_tools:<tag> bash /opt/meco_tools/tests/smoke_test.sh
# Checks that every script's language dependencies resolve; it does not exercise
# external bioinformatics tools, which are intentionally not part of the image.
set -uo pipefail

ROOT=${MECO_TOOLS_HOME:-/opt/meco_tools}
fail=0
pass() { printf '  ok    %s\n' "$1"; }
bad()  { printf '  FAIL  %s\n' "$1"; fail=1; }

echo "== runtimes"
if command -v python2 >/dev/null || ls /opt/conda/bin/python2* >/dev/null 2>&1; then bad "python2 present"; else pass "no python2"; fi
python --version 2>&1 | grep -q '^Python 3\.' && pass "$(python --version 2>&1)" || bad "python is not python3"
pass "$(perl -e 'printf "perl %vd", $^V')"
pass "$(R --version | head -1)"
[ "$(command -v perl)" = /usr/bin/perl ] && pass "perl is Debian's /usr/bin/perl" || bad "perl resolves to $(command -v perl), expected /usr/bin/perl"
for b in ps pigz wget gzip awk bash R Rscript; do command -v "$b" >/dev/null && pass "$b" || bad "missing $b"; done
echo "== external tools must NOT be in the image"
for b in samtools blastn blastp diamond vsearch usearch bwa bowtie2 hmmscan prodigal megahit java; do
  command -v "$b" >/dev/null && bad "unexpected external tool: $b ($(command -v "$b"))"
done
pass "no external bioinformatics tools found"

echo "== perl: compile every script (resolves all 'use' statements)"
n=0
for f in "$ROOT"/AmpliconTagger/*.pl "$ROOT"/perl/scripts/*.pl; do
  if out=$(perl -c "$f" 2>&1); then n=$((n+1))
  else bad "$(basename "$f"): $(echo "$out" | grep -m1 -E "Can't locate|error|syntax|Global symbol" | cut -c1-160)"; fi
done
pass "$n perl scripts compile"

echo "== python: compile every script and import its third-party modules"
python - "$ROOT" <<'EOF' || fail=1
import ast, importlib, os, sys, warnings
root = sys.argv[1]
sys.path.insert(0, os.path.join(root, "python"))  # for the local 'nrc' package
files = [os.path.join(d, f) for sub in ("AmpliconTagger", "python", "python/nrc")
         for d in [os.path.join(root, sub)] for f in sorted(os.listdir(d)) if f.endswith(".py")]
mods, rc = set(), 0
for f in files:
    src = open(f).read()
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("error", SyntaxWarning)
            tree = ast.parse(src, f)
    except (SyntaxError, SyntaxWarning) as e:
        print("  FAIL  %s: %s" % (f, e)); rc = 1; continue
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            mods.update(a.name for a in node.names)
        elif isinstance(node, ast.ImportFrom) and node.level == 0 and node.module:
            mods.add(node.module)
mods.discard("__future__")
for m in sorted(mods):
    try:
        importlib.import_module(m)
    except Exception as e:
        print("  FAIL  import %s: %s" % (m, e)); rc = 1
print("  ok    %d python files parse, %d modules import" % (len(files), len(mods)))
sys.exit(rc)
EOF

echo "== R: parse every script and load every library it uses"
Rscript --vanilla - "$ROOT" <<'EOF' || fail=1
root  <- commandArgs(trailingOnly = TRUE)[1]
files <- c(list.files(file.path(root, "R"), "\\.R$", full.names = TRUE),
           list.files(file.path(root, "AmpliconTagger"), "\\.R$", full.names = TRUE))
rc <- 0; pkgs <- character()
for (f in files) {
  ok <- tryCatch({ parse(f); TRUE }, error = function(e) { cat("  FAIL ", f, conditionMessage(e), "\n"); FALSE })
  if (!ok) rc <- 1
  src <- readLines(f, warn = FALSE)
  m <- regmatches(src, regexpr("(library|require)\\(\\s*['\"]?[A-Za-z0-9.]+", src))
  pkgs <- c(pkgs, sub("^.*\\(\\s*['\"]?", "", m))
}
pkgs <- sort(unique(pkgs))
for (p in pkgs) {
  if (!suppressWarnings(suppressPackageStartupMessages(requireNamespace(p, quietly = TRUE)))) {
    cat("  FAIL  R package", p, "\n"); rc <- 1
  }
}
cat(sprintf("  ok    %d R scripts parse, %d packages load\n", length(files), length(pkgs)))
quit(status = rc)
EOF

echo
if [ "$fail" -eq 0 ]; then echo "SMOKE TEST PASSED"; else echo "SMOKE TEST FAILED"; fi
exit "$fail"
