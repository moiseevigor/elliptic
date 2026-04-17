#!/usr/bin/env bash
# update_gh_pages.sh — Regenerate the GitHub Pages site from docs/gh-pages-body.md
#
# Usage:
#   ./scripts/update_gh_pages.sh              # build and push
#   ./scripts/update_gh_pages.sh --dry-run    # build only, inspect /tmp/gh-pages-work
#
# Requirements: git, pandoc, python3
# The markdown source lives in docs/gh-pages-body.md — edit that file, then run
# this script to publish.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "$0")/.." && pwd)"
WORKTREE_DIR="/tmp/gh-pages-work"
SOURCE_MD="$REPO_ROOT/docs/gh-pages-body.md"
DRY_RUN=0

# ── parse args ────────────────────────────────────────────────────────────────
for arg in "$@"; do
  case "$arg" in
    --dry-run|-n) DRY_RUN=1 ;;
    --help|-h)
      echo "Usage: $0 [--dry-run]"
      echo "  Edit docs/gh-pages-body.md, then run this script to update the site."
      exit 0 ;;
    *) echo "Unknown argument: $arg"; exit 1 ;;
  esac
done

# ── prerequisites ─────────────────────────────────────────────────────────────
for cmd in git pandoc python3; do
  command -v "$cmd" >/dev/null || { echo "ERROR: $cmd not found"; exit 1; }
done

[[ -f "$SOURCE_MD" ]] || { echo "ERROR: $SOURCE_MD not found"; exit 1; }

# ── prepare worktree ─────────────────────────────────────────────────────────
cd "$REPO_ROOT"

if [[ -d "$WORKTREE_DIR" ]]; then
  echo "Removing existing worktree at $WORKTREE_DIR"
  git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || rm -rf "$WORKTREE_DIR"
fi

echo "Creating worktree for gh-pages at $WORKTREE_DIR"
git fetch origin gh-pages
git worktree add "$WORKTREE_DIR" origin/gh-pages

# ── convert markdown → HTML fragment ─────────────────────────────────────────
echo "Converting $SOURCE_MD to HTML"
pandoc --from markdown --to html5 --no-highlight "$SOURCE_MD" -o /tmp/gh-pages-fragment.html

# ── fix heading anchors to match GitHub Pages Dinky theme ────────────────────
# pandoc:        <h1 id="foo">Text</h1>
# Dinky theme:   <h1>\n<a id="foo" class="anchor" ...></a>Text</h1>
python3 - <<'PYEOF'
import re, sys

with open('/tmp/gh-pages-fragment.html', 'r') as f:
    body = f.read()

def fix_heading(m):
    tag  = m.group(1)
    attrs = m.group(2)
    text  = m.group(3)
    id_m  = re.search(r'id="([^"]+)"', attrs)
    if id_m:
        aid = id_m.group(1)
        anchor = (f'<a id="{aid}" class="anchor" href="#{aid}" aria-hidden="true">'
                  f'<span aria-hidden="true" class="octicon octicon-link"></span></a>')
        return f'<{tag}>\n{anchor}{text}</{tag}>'
    return m.group(0)

body = re.sub(r'<(h[1-6])([^>]*)>(.*?)</\1>', fix_heading, body, flags=re.DOTALL)

with open('/tmp/gh-pages-fragment-fixed.html', 'w') as f:
    f.write(body)
PYEOF

# ── assemble full index.html ─────────────────────────────────────────────────
echo "Assembling index.html"
python3 - <<PYEOF
with open('$WORKTREE_DIR/index.html', 'r') as f:
    template = f.read()

marker_open  = '<section>'
marker_close = '</section>'
header = template[:template.index(marker_open) + len(marker_open)]
footer = template[template.index(marker_close):]

with open('/tmp/gh-pages-fragment-fixed.html', 'r') as f:
    body = f.read()

with open('$WORKTREE_DIR/index.html', 'w') as f:
    f.write(header + '\n' + body + '\n' + footer)

print(f"index.html: {len(header)+len(body)+len(footer)} chars")
PYEOF

# ── update params.json ────────────────────────────────────────────────────────
echo "Updating params.json"
python3 - <<PYEOF
import json

with open('$WORKTREE_DIR/params.json', 'r') as f:
    params = json.load(f)

with open('$SOURCE_MD', 'r') as f:
    md = f.read()

params['body'] = md.replace('\n', '\r\n')

with open('$WORKTREE_DIR/params.json', 'w') as f:
    json.dump(params, f, ensure_ascii=False)

print("params.json updated")
PYEOF

# ── copy examples directory ──────────────────────────────────────────────────
echo "Copying examples/"
EXAMPLES_SRC="$REPO_ROOT/examples"
EXAMPLES_DST="$WORKTREE_DIR/examples"
if [[ -d "$EXAMPLES_SRC" ]]; then
  rm -rf "$EXAMPLES_DST"
  cp -r "$EXAMPLES_SRC" "$EXAMPLES_DST"
  echo "  Copied $(find "$EXAMPLES_DST" -type f | wc -l) files to examples/"
else
  echo "  No examples/ directory found — skipping."
fi

# ── sanity check ─────────────────────────────────────────────────────────────
echo "Sanity check:"
python3 -c "
import re
with open('$WORKTREE_DIR/index.html') as f:
    html = f.read()
h1 = len(re.findall(r'<h1', html))
h2 = len(re.findall(r'<h2', html))
new_fns = ['weierstrassP','carlsonRF','ellipticBDJ','jacobiEDJ','cel']
missing = [fn for fn in new_fns if fn not in html]
print(f'  h1={h1}  h2={h2}  missing={missing or \"none\"}')
if missing:
    raise SystemExit('ERROR: expected functions missing from output')
"

# ── commit and push (or dry-run) ─────────────────────────────────────────────
cd "$WORKTREE_DIR"
git add index.html params.json
[[ -d "$EXAMPLES_DST" ]] && git add examples/

if [[ $DRY_RUN -eq 1 ]]; then
  echo "DRY RUN — skipping commit/push.  Inspect $WORKTREE_DIR to review changes."
  git diff --stat HEAD
else
  if git diff --cached --quiet; then
    echo "Nothing changed — gh-pages is already up to date."
  else
    MSG="docs: update gh-pages from docs/gh-pages-body.md

Source: $(cd "$REPO_ROOT" && git log -1 --format='%h %s' HEAD)"
    git commit -m "$MSG"
    git push origin HEAD:gh-pages
    echo "✓ Pushed to gh-pages — https://moiseevigor.github.io/elliptic/"
  fi
fi

# ── cleanup ──────────────────────────────────────────────────────────────────
cd "$REPO_ROOT"
git worktree remove --force "$WORKTREE_DIR" 2>/dev/null || true
rm -f /tmp/gh-pages-fragment.html /tmp/gh-pages-fragment-fixed.html
echo "Done."
