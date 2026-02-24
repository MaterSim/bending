#!/usr/bin/env python3
"""Generate HTML gallery for a timestamped folder and update index.html."""

from __future__ import annotations

import argparse
from pathlib import Path
from datetime import datetime


def generate_gallery(outdir: Path) -> Path:
    pngs = sorted(outdir.glob("*.png"))
    title = f"Figures: {outdir.name}"
    html = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "  <meta charset='utf-8'>",
        f"  <title>{title}</title>",
        "  <style>",
        "    body { font-family: Arial, sans-serif; margin: 20px; }",
        "    .grid { display: grid; grid-template-columns: repeat(auto-fill, minmax(360px, 1fr)); gap: 16px; }",
        "    figure { margin: 0; padding: 8px; border: 1px solid #ddd; border-radius: 8px; }",
        "    img { width: 100%; height: auto; }",
        "    figcaption { margin-top: 6px; font-size: 0.9em; color: #333; }",
        "  </style>",
        "</head>",
        "<body>",
        f"  <h1>{title}</h1>",
        "  <div class='grid'>",
    ]

    for png in pngs:
        html.extend([
            "    <figure>",
            f"      <img src='{png.name}' alt='{png.name}'>",
            f"      <figcaption>{png.name}</figcaption>",
            "    </figure>",
        ])

    html.extend([
        "  </div>",
        "</body>",
        "</html>",
    ])

    out_path = outdir / "index.html"
    out_path.write_text("\n".join(html), encoding="utf-8")
    return out_path


def generate_root_index(fig_root: Path) -> Path:
    folders = sorted([p for p in fig_root.iterdir() if p.is_dir()], reverse=True)
    html = [
        "<!DOCTYPE html>",
        "<html>",
        "<head>",
        "  <meta charset='utf-8'>",
        "  <title>Figure Index</title>",
        "  <style>",
        "    body { font-family: Arial, sans-serif; margin: 20px; }",
        "    ul { line-height: 1.8; }",
        "  </style>",
        "</head>",
        "<body>",
        "  <h1>Figure Index</h1>",
        "  <ul>",
    ]

    for folder in folders:
        try:
            ts = datetime.strptime(folder.name, "%Y%m%d_%H%M%S")
            label = ts.strftime("%Y-%m-%d %H:%M:%S")
        except ValueError:
            label = folder.name
        html.append(f"    <li><a href='{folder.name}/index.html'>{label}</a></li>")

    html.extend([
        "  </ul>",
        "</body>",
        "</html>",
    ])

    out_path = fig_root / "index.html"
    out_path.write_text("\n".join(html), encoding="utf-8")
    return out_path


def main():
    parser = argparse.ArgumentParser(description="Generate HTML for figures")
    parser.add_argument("--outdir", required=True, help="Timestamped folder with PNGs")
    parser.add_argument("--index", action="store_true", help="Update root index.html")
    args = parser.parse_args()

    outdir = Path(args.outdir).resolve()
    if not outdir.exists():
        raise SystemExit(f"Outdir not found: {outdir}")

    gallery = generate_gallery(outdir)
    print(f"Wrote {gallery}")

    if args.index:
        fig_root = outdir.parent
        index = generate_root_index(fig_root)
        print(f"Wrote {index}")


if __name__ == "__main__":
    main()
