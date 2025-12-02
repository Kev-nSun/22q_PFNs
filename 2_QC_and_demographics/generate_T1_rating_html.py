#!/usr/bin/env python3
"""
Generate a self-contained HTML page for T1 QC ratings from slice PNGs.

Features:
- Scan an image directory for `sub-<ID>_ses-<ID>_<VIEW>.png` images
- Group images by (subject, session) for expected views (default: S1, S3, A2, A3)
- Produce an HTML with a simple UI to rate each view per (sub, ses)
- Export ratings to CSV client-side (no backend required)
- Optional portability: copy referenced images next to the HTML (assets/) and rewrite paths

Example:
  python /cbica/projects/bbl_22q/data/scripts/T1_QC_scripts/generate_T1_rating_html.py \
    --root /cbica/projects/bbl_22q/data/T1_QC/slices \
    --out /cbica/projects/bbl_22q/data/T1_QC/22q_T1_QC_ratings.html \
    --portable \
    --allow-missing

"""

import argparse
import base64
import glob
import os
import re
import shutil
from typing import Dict, List, Tuple
import json


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Generate T1 QC rating HTML from slice images"
    )
    parser.add_argument(
        "--root",
        default=None,
        help=(
            "Directory containing PNG images (sub-*_ses-*_VIEW.png). "
            "Required unless --generate-from-nifti is used; in that case, "
            "defaults to --png-outdir."
        ),
    )
    parser.add_argument("--out", required=True, help="Output HTML file path")
    parser.add_argument(
        "--views",
        nargs="+",
        default=["S1", "S3", "A2", "A3"],
        help="List of expected views per row",
    )
    parser.add_argument(
        "--pattern",
        default="sub-*_ses-*_*.png",
        help="Glob pattern (relative to --root) for candidate images",
    )
    parser.add_argument(
        "--allow-missing",
        action="store_true",
        help="Include rows even if not all expected views are present",
    )
    parser.add_argument(
        "--portable",
        action="store_true",
        help="Copy used images into assets/ next to HTML and rewrite paths",
    )
    parser.add_argument(
        "--embed-small",
        type=int,
        default=0,
        help=(
            "If >0, embed images smaller than this many KB as base64 in HTML. "
            "Larger images will be linked/copied as files."
        ),
    )
    return parser.parse_args()


FilenameRegex = re.compile(r"^(sub-[^_]+)_(ses-[^_]+)_([^.]+)\.png$")


def natural_sort_key(sub: str, ses: str) -> Tuple:
    m_sub = re.search(r"sub-(\d+)$", sub)
    m_ses = re.search(r"ses-(\d+)$", ses)
    return (
        int(m_sub.group(1)) if m_sub else sub,
        int(m_ses.group(1)) if m_ses else ses,
    )


def collect_rows(
    root: str, pattern: str, views: List[str], allow_missing: bool
) -> Tuple[List[Tuple[str, str]], Dict[Tuple[str, str], Dict[str, str]]]:
    pairs: Dict[Tuple[str, str], Dict[str, str]] = {}
    for path in glob.glob(os.path.join(root, pattern)):
        fn = os.path.basename(path)
        m = FilenameRegex.match(fn)
        if not m:
            continue
        sub, ses, view = m.groups()
        if view not in views:
            continue
        pairs.setdefault((sub, ses), {})[view] = path

    keys = sorted(pairs.keys(), key=lambda k: natural_sort_key(k[0], k[1]))
    if not allow_missing:
        keys = [k for k in keys if all(v in pairs[k] for v in views)]
    return keys, pairs


def ensure_portable_assets(
    keys: List[Tuple[str, str]],
    pairs: Dict[Tuple[str, str], Dict[str, str]],
    out_html: str,
    embed_small_kb: int,
) -> Tuple[Dict[Tuple[str, str], Dict[str, str]], List[str]]:
    """
    Prepare image references for portability:
    - Create assets/ next to HTML
    - Copy images there and rewrite paths relative to HTML
    - For small images (<= embed_small_kb), embed as data URIs
    Returns updated mapping and list of copied file paths
    """
    copied: List[str] = []
    assets_dir = os.path.join(os.path.dirname(os.path.abspath(out_html)), "assets")
    os.makedirs(assets_dir, exist_ok=True)

    updated: Dict[Tuple[str, str], Dict[str, str]] = {}
    for key in keys:
        updated[key] = {}
        for view, src in pairs[key].items():
            try:
                size_kb = os.path.getsize(src) / 1024.0
            except OSError:
                size_kb = float("inf")

            if embed_small_kb > 0 and size_kb <= embed_small_kb:
                try:
                    with open(src, "rb") as f:
                        b64 = base64.b64encode(f.read()).decode("ascii")
                    updated[key][view] = f"data:image/png;base64,{b64}"
                except Exception:
                    # Fallback to copying
                    dst = os.path.join(assets_dir, os.path.basename(src))
                    shutil.copy2(src, dst)
                    copied.append(dst)
                    updated[key][view] = os.path.relpath(dst, os.path.dirname(out_html))
            else:
                dst = os.path.join(assets_dir, os.path.basename(src))
                shutil.copy2(src, dst)
                copied.append(dst)
                updated[key][view] = os.path.relpath(dst, os.path.dirname(out_html))

    return updated, copied


def make_relative_refs(
    keys: List[Tuple[str, str]],
    pairs: Dict[Tuple[str, str], Dict[str, str]],
    out_html: str,
) -> Dict[Tuple[str, str], Dict[str, str]]:
    updated: Dict[Tuple[str, str], Dict[str, str]] = {}
    out_dir = os.path.dirname(os.path.abspath(out_html))
    for key in keys:
        updated[key] = {}
        for view, src in pairs[key].items():
            updated[key][view] = os.path.relpath(os.path.abspath(src), out_dir)
    return updated


def render_html(
    keys: List[Tuple[str, str]],
    pairs: Dict[Tuple[str, str], Dict[str, str]],
    views: List[str],
    out_path: str,
):
    # Minimal CSS + JS for table layout and CSV export
    css = """
    body { font-family: system-ui, -apple-system, Segoe UI, Roboto, Arial, sans-serif; margin: 16px; }
    h1 { font-size: 20px; margin: 0 0 12px 0; }
    .controls { display: flex; gap: 8px; align-items: center; margin-bottom: 12px; flex-wrap: wrap; }
    .grid { display: grid; grid-template-columns: 1fr; gap: 16px; }
    .card { border: 1px solid #ddd; border-radius: 8px; padding: 8px; background: #fff; }
    .header { display: flex; align-items: baseline; justify-content: space-between; gap: 8px; margin-bottom: 8px; }
    .views { display: grid; grid-template-columns: repeat(4, minmax(0, 1fr)); gap: 12px; }
    .view { display: flex; flex-direction: column; gap: 6px; }
    .view img { width: 100%; height: auto; display: block; border: 1px solid #eee; border-radius: 4px; }
    .rating { display: flex; gap: 6px; align-items: center; flex-wrap: wrap; }
    .muted { color: #666; font-size: 12px; }
    .search { padding: 6px 10px; border: 1px solid #ccc; border-radius: 6px; min-width: 260px; }
    button { padding: 6px 10px; border: 1px solid #aaa; border-radius: 6px; background: #f7f7f7; cursor: pointer; }
    button.primary { background: #0b5; color: #fff; border-color: #0a4; }
    button:disabled { opacity: 0.6; cursor: not-allowed; }
    """

    # Data for client-side JS
    # Build rows: [{sub, ses, images: {view: src}}]
    data_rows = []
    for sub, ses in keys:
        images = {v: pairs[(sub, ses)].get(v, "") for v in views}
        data_rows.append(
            {
                "sub": sub,
                "ses": ses,
                "images": images,
            }
        )

    # Serialize minimal JSON safely
    json_data = json.dumps(data_rows)
    json_views = json.dumps(views)

    js = f"""
    const DATA = {json_data};
    const VIEWS = {json_views};
    let PREFILL = {{}}; // key: `${{sub}}|${{ses}}` -> {{ view: score }}
    let RATINGS = {{}}; // live user selections, same keying as PREFILL

    const makeKey = (sub, ses) => sub + '|' + ses;

    function escapeCsv(val) {{
      if (val == null) return '';
      const s = String(val);
      if (s.includes('"') || s.includes(',') || s.includes('\\n')) {{
        return '"' + s.replace(/"/g, '""') + '"';
      }}
      return s;
    }}

    function toCSV() {{
      const headers = ['subid','sesid', ...VIEWS.map(v => v + '_score')];
      const lines = [headers.join(',')];
      for (const row of DATA) {{
        const key = makeKey(row.sub, row.ses);
        const scores = VIEWS.map(v => {{
          if (RATINGS[key] && RATINGS[key][v] !== undefined) return RATINGS[key][v];
          if (PREFILL[key] && PREFILL[key][v] !== undefined) return PREFILL[key][v];
          return '';
        }});
        const line = [row.sub, row.ses, ...scores].map(escapeCsv).join(',');
        lines.push(line);
      }}
      return lines.join('\\n');
    }}

    function downloadCSV() {{
      const csv = toCSV();
      const blob = new Blob([csv], {{ type: 'text/csv;charset=utf-8;' }});
      const url = URL.createObjectURL(blob);
      const a = document.createElement('a');
      a.href = url;
      a.download = 'T1-ratings.csv';
      document.body.appendChild(a);
      a.click();
      document.body.removeChild(a);
      URL.revokeObjectURL(url);
    }}

    function parseCSV(text) {{
      const rows = [];
      let i = 0, field = '', inQuotes = false, row = [];
      while (i < text.length) {{
        const c = text[i++];
        if (inQuotes) {{
          if (c === '"') {{
            if (text[i] === '"') {{ field += '"'; i++; }} else {{ inQuotes = false; }}
          }} else {{
            field += c;
          }}
        }} else {{
          if (c === '"') {{ inQuotes = true; }}
          else if (c === ',') {{ row.push(field); field=''; }}
          else if (c === '\\n' || c === '\\r') {{
            if (c === '\\r' && text[i] === '\\n') i++;
            row.push(field); rows.push(row); row = []; field='';
          }} else {{ field += c; }}
        }}
      }}
      if (field.length || row.length) {{ row.push(field); rows.push(row); }}
      return rows.filter(r => r.length && r.some(x => x !== ''));
    }}

    function buildPrefillFromCSV(rows) {{
      if (!rows || rows.length < 2) return {{}};
      const header = rows[0];
      const subIdx = header.indexOf('subid');
      const sesIdx = header.indexOf('sesid');
      const viewIdx = Object.fromEntries(VIEWS.map(v => [v, header.indexOf(v + '_score')]));
      const map = {{}};
      for (let r = 1; r < rows.length; r++) {{
        const row = rows[r];
        const sub = row[subIdx];
        const ses = row[sesIdx];
        if (!sub || !ses) continue;
        const key = `${{sub}}|${{ses}}`;
        map[key] = map[key] || {{}};
        for (const v of VIEWS) {{
          const idx = viewIdx[v];
          if (idx >= 0 && row[idx] !== undefined && row[idx] !== '') {{
            map[key][v] = String(row[idx]);
          }}
        }}
      }}
      return map;
    }}

    function applyPrefill(prefill) {{
      const cards = document.querySelectorAll('[data-card]');
      let applied = 0;
      for (const card of cards) {{
        const sub = card.getAttribute('data-sub');
        const ses = card.getAttribute('data-ses');
        const key = `${{sub}}|${{ses}}`;
        const views = prefill[key];
        if (!views) continue;
        for (const v of VIEWS) {{
          const val = views[v];
          if (val === undefined) continue;
          const sel = card.querySelector(`[data-score="${{v}}"]`);
          if (sel) {{ sel.value = val; applied++; }}
        }}
      }}
      return applied;
    }}

    function onCsvSelected(file) {{
      const reader = new FileReader();
      reader.onload = (e) => {{
        const text = e.target.result;
        const rows = parseCSV(text);
        PREFILL = buildPrefillFromCSV(rows);
        const n = applyPrefill(PREFILL);
        const msg = document.getElementById('resumeMsg');
        if (msg) msg.textContent = `Applied ${{n}} ratings from CSV`;
      }};
      reader.readAsText(file);
    }}

    function render(filter='') {{
      const root = document.getElementById('root');
      root.innerHTML = '';
      const q = filter.trim().toLowerCase();
      const rows = DATA.filter(r => (r.sub + ' ' + r.ses).toLowerCase().includes(q));
      for (const row of rows) {{
        const card = document.createElement('div');
        card.className = 'card';
        card.setAttribute('data-card', '');
        card.setAttribute('data-sub', row.sub);
        card.setAttribute('data-ses', row.ses);

        const header = document.createElement('div');
        header.className = 'header';
        const title = document.createElement('div');
        title.innerHTML = `<strong>${{row.sub}}</strong> &nbsp; <span class="muted">${{row.ses}}</span>`;
        const hint = document.createElement('div');
        hint.className = 'muted';
        hint.textContent = 'Rate each view: 0=Fail, 1=Pass';
        header.appendChild(title);
        header.appendChild(hint);
        card.appendChild(header);

        const views = document.createElement('div');
        views.className = 'views';
        for (const v of VIEWS) {{
          const box = document.createElement('div');
          box.className = 'view';
          const img = document.createElement('img');
          const src = row.images[v] || '';
          if (src) img.src = src;
          img.alt = `${{row.sub}} ${{row.ses}} ${{v}}`;
          const label = document.createElement('label');
          label.textContent = v;
          const select = document.createElement('select');
          select.setAttribute('data-score', v);
          for (const [val, name] of [['', ''], ['0','Fail'], ['1','Pass']]) {{
            const opt = document.createElement('option');
            opt.value = val; opt.textContent = name;
            select.appendChild(opt);
          }}
          // initialize from RATINGS or PREFILL
          const key = makeKey(row.sub, row.ses);
          let initVal = '';
          if (RATINGS[key] && RATINGS[key][v] !== undefined) {{ initVal = RATINGS[key][v]; }}
          else if (PREFILL[key] && PREFILL[key][v] !== undefined) {{ initVal = PREFILL[key][v]; }}
          select.value = initVal;
          // persist on change
          select.addEventListener('change', () => {{
            RATINGS[key] = RATINGS[key] || {{}};
            RATINGS[key][v] = select.value;
          }});
          const rating = document.createElement('div');
          rating.className = 'rating';
          rating.appendChild(label);
          rating.appendChild(select);
          box.appendChild(img);
          box.appendChild(rating);
          views.appendChild(box);
        }}
        card.appendChild(views);
        root.appendChild(card);
      }}
      document.getElementById('count').textContent = rows.length + ' rows';
    }}

    window.addEventListener('DOMContentLoaded', () => {{
      const search = document.getElementById('search');
      search.addEventListener('input', () => render(search.value));
      const csv = document.getElementById('csvUpload');
      csv.addEventListener('change', () => {{ if (csv.files && csv.files[0]) onCsvSelected(csv.files[0]); }});
      render('');
    }});
    """

    html_doc = f"""
<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <meta name="viewport" content="width=device-width, initial-scale=1" />
  <title>T1 QC Ratings</title>
  <style>{css}</style>
  <script>{js}</script>
</head>
<body>
  <h1>T1 QC Ratings</h1>
  <div class="controls">
    <input id="search" class="search" type="search" placeholder="Filter by subject/session (e.g., sub-123 ses-1)" />
    <button class="primary" onclick="downloadCSV()">Download CSV</button>
    <input id="csvUpload" type="file" accept=".csv,text/csv" />
    <span id="resumeMsg" class="muted"></span>
    <span id="count" class="muted"></span>
  </div>
  <div id="root" class="grid"></div>
</body>
</html>
"""

    with open(out_path, "w", encoding="utf-8") as f:
        f.write(html_doc)


def main() -> None:
    args = parse_args()

    keys, pairs = collect_rows(args.root, args.pattern, args.views, args.allow_missing)

    if not keys:
        raise SystemExit(
            "No rows found. Check --root and --pattern, and ensure filenames match 'sub-*_ses-*_VIEW.png'."
        )

    # Rewrite image references according to portability settings
    if args.portable:
        pairs, copied = ensure_portable_assets(keys, pairs, args.out, args.embed_small)
        print(
            f"Prepared portable assets in 'assets/' next to HTML. Files copied: {len(copied)}"
        )
    else:
        pairs = make_relative_refs(keys, pairs, args.out)

    os.makedirs(os.path.dirname(os.path.abspath(args.out)) or ".", exist_ok=True)
    render_html(keys, pairs, args.views, args.out)
    print(f"Wrote HTML: {args.out}")
    print(f"Rows: {len(keys)} | Views per row: {len(args.views)}")
    if args.portable:
        print(
            "Note: For sharing, distribute the HTML and the 'assets/' folder together."
        )
    else:
        print(
            "Note: HTML references images by relative paths; receivers need the same directory structure."
        )


if __name__ == "__main__":
    main()
