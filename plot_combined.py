#!/usr/bin/env python3
"""Plot energy and neighbor distances in a single figure.

Usage:
  python3 plot_combined.py --dump bend.dump [--log log.lammps]

Outputs: combined_plot.png
"""

import sys
import re
import argparse
from pathlib import Path
import tempfile
import os
import subprocess

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import image as mpimg


def smooth_data(data, window=5):
    """Apply moving average smoothing."""
    if len(data) < window:
        return data
    kernel = np.ones(window) / window
    return np.convolve(data, kernel, mode='same')


def parse_log(path: Path):
    """Parse log.lammps for energy data."""
    step = []
    toteng = []
    poteng = []
    zmin = []

    header_re = re.compile(r"\bStep\b.*\bTotEng\b.*\bPotEng\b")

    with path.open(encoding='utf-8', errors='ignore') as f:
        in_block = False
        col_idx = {}
        for line in f:
            line = line.strip()
            if not line:
                in_block = False
                col_idx = {}
                continue

            if header_re.search(line):
                headers = line.split()
                col_idx = {name: i for i, name in enumerate(headers)}
                in_block = True
                continue

            if in_block:
                parts = line.split()
                if len(parts) < 3:
                    in_block = False
                    continue
                try:
                    s = float(parts[col_idx["Step"]])
                    te = float(parts[col_idx["TotEng"]])
                    pe = float(parts[col_idx["PotEng"]])
                    zm = float(parts[col_idx["c_zmin"]])
                except Exception:
                    in_block = False
                    continue

                step.append(s)
                toteng.append(te)
                poteng.append(pe)
                zmin.append(zm)

    return step, toteng, poteng, zmin


def parse_dump(path: Path):
    """Parse dump file for neighbor distances."""
    with path.open(encoding='utf-8', errors='ignore') as f:
        while True:
            line = f.readline()
            if not line:
                break
            if not line.startswith("ITEM: TIMESTEP"):
                continue
            timestep = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: NUMBER OF ATOMS"):
                raise ValueError("Unexpected dump format: missing NUMBER OF ATOMS")
            natoms = int(f.readline().strip())

            line = f.readline()
            if not line.startswith("ITEM: BOX BOUNDS"):
                raise ValueError("Unexpected dump format: missing BOX BOUNDS")
            f.readline(); f.readline(); f.readline()

            line = f.readline()
            if not line.startswith("ITEM: ATOMS"):
                raise ValueError("Unexpected dump format: missing ATOMS header")

            cols = line.split()[2:]
            col_idx = {name: i for i, name in enumerate(cols)}
            required = ["id", "type", "x", "y", "z"]
            for r in required:
                if r not in col_idx:
                    raise ValueError(f"Dump missing required column: {r}")

            atoms = {}
            for _ in range(natoms):
                parts = f.readline().split()
                if not parts:
                    continue
                aid = int(parts[col_idx["id"]])
                atype = int(parts[col_idx["type"]])
                x = float(parts[col_idx["x"]])
                y = float(parts[col_idx["y"]])
                z = float(parts[col_idx["z"]])
                atoms[aid] = (atype, x, y, z)

            yield timestep, atoms


def dist(a, b):
    """Compute distance between two atoms."""
    dx = a[1] - b[1]
    dy = a[2] - b[2]
    dz = a[3] - b[3]
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def compute_neighbor_distances(dump_path):
    """Extract mean neighbor distances from dump file."""
    steps = []
    data_by_layer = {}  # layer -> list of (step_idx, mean_distance)
    first_id_by_layer = {}
    all_layers = set()

    for step_idx, (step, atoms) in enumerate(parse_dump(dump_path)):
        steps.append(step)

        # Map: layer(type) -> atom id with minimum z in that layer
        target_ids = {}
        for aid, (atype, x, y, z) in atoms.items():
            all_layers.add(atype)
            if atype not in target_ids or z < atoms[target_ids[atype]][3]:
                target_ids[atype] = aid

        for layer_key, target_id in target_ids.items():
            prev_id = target_id - 1
            next_id = target_id + 1

            if target_id not in atoms or prev_id not in atoms or next_id not in atoms:
                continue

            target = atoms[target_id]
            prev_atom = atoms[prev_id]
            next_atom = atoms[next_id]

            if prev_atom[0] != target[0] or next_atom[0] != target[0]:
                continue

            d1 = dist(target, prev_atom)
            d2 = dist(target, next_atom)
            mean_d = max(d1, d2)

            data_by_layer.setdefault(layer_key, []).append((step_idx, mean_d))
            first_id_by_layer[layer_key] = target_id

    # Create full arrays with 3.0 for missing data
    n_steps = len(steps)
    means_by_layer = {}
    for layer in all_layers:
        means_by_layer[layer] = np.full(n_steps, 3.0)
        if layer in data_by_layer:
            for step_idx, mean_d in data_by_layer[layer]:
                means_by_layer[layer][step_idx] = mean_d

    return steps, means_by_layer, first_id_by_layer


def render_last_frame(dump_path, output_path, width=1600, height=800):
    """Render the last frame from dump file using external render script."""
    print(f"Rendering last frame from {dump_path}...")
    
    # Find render_frames.py script
    script_dir = Path(__file__).parent
    render_script = script_dir / "render_frames.py"
    
    if not render_script.exists():
        print(f"Warning: render_frames.py not found at {render_script}")
        # Create a placeholder image
        fig, ax = plt.subplots(figsize=(width/100, height/100), dpi=100)
        ax.text(0.5, 0.5, "render_frames.py not found", 
                ha='center', va='center', fontsize=12)
        ax.axis('off')
        fig.savefig(output_path, bbox_inches='tight')
        plt.close(fig)
        return output_path
    
    try:
        # Create temp directory for frame output
        temp_dir = tempfile.mkdtemp()
        
        # Render with a reasonable step to ensure we get the last frame
        # Use step=100 to render a few frames including the last one
        result = subprocess.run(
            [sys.executable, str(render_script), str(dump_path), temp_dir, 
             str(width), str(height), "100"],  # Render every 100th frame
            capture_output=True,
            text=True,
            timeout=60
        )
        
        if result.returncode != 0:
            raise Exception(f"render_frames.py failed: {result.stderr}")
        
        # Find all rendered frames and use the last one
        frame_files = sorted(Path(temp_dir).glob("frame_*.png"))
        if not frame_files:
            raise Exception("No frames were rendered")
        
        # Use the last frame (which should be the last frame of the simulation)
        last_frame_path = frame_files[-1]
        
        # Copy to output path
        import shutil
        shutil.copy(last_frame_path, output_path)
        
        # Cleanup temp directory
        shutil.rmtree(temp_dir, ignore_errors=True)
        
        print(f"Last frame rendered to {output_path}")
        return output_path
        
    except Exception as e:
        print(f"Warning: Could not render frame: {e}")
        # Create a placeholder image
        fig, ax = plt.subplots(figsize=(width/100, height/100), dpi=100)
        ax.text(0.5, 0.5, f"Frame rendering failed\n{str(e)[:50]}", 
                ha='center', va='center', fontsize=10, wrap=True)
        ax.axis('off')
        fig.savefig(output_path, bbox_inches='tight')
        plt.close(fig)
        return output_path


def main():
    parser = argparse.ArgumentParser(description="Plot combined energy and neighbor distances")
    parser.add_argument("--dump", "-d", required=True, help="LAMMPS dump file path")
    parser.add_argument("--log", "-l", default="log.lammps", help="LAMMPS log file (default: log.lammps)")

    args = parser.parse_args()

    dump_path = Path(args.dump)
    log_path = Path(args.log)

    if not dump_path.exists():
        raise SystemExit(f"Dump file not found: {dump_path}")
    if not log_path.exists():
        raise SystemExit(f"Log file not found: {log_path}")

    # Parse both data sources
    step_log, toteng, poteng, zmin = parse_log(log_path)
    steps_dump, means_by_layer, first_id_by_layer = compute_neighbor_distances(dump_path)

    if not step_log:
        raise SystemExit("No thermo data found (Step/TotEng/PotEng).")
    if not steps_dump:
        raise SystemExit("No dump data found.")

    # Render last frame to temporary file
    with tempfile.NamedTemporaryFile(suffix='.png', delete=False) as tmp:
        frame_img_path = tmp.name
    render_last_frame(dump_path, frame_img_path, width=1280, height=800)
    
    # Create figure with 3x1 subplots with height ratios 1:1:1.8
    fig, (ax1, ax3, ax5) = plt.subplots(3, 1, figsize=(8, 8), 
                                        gridspec_kw={'height_ratios': [1, 1, 1.2]})

    # === Top panel: Energy ===
    ax1.plot(step_log, smooth_data(toteng, window=1), label="TotEng", linewidth=1.5)
    ax1.plot(step_log, smooth_data(poteng, window=1), label="PotEng", linewidth=1.5)
    ax1.set_ylabel("Energy (kcal/mol)")
    ax1.set_title("Energy Evolution")

    ax2 = ax1.twinx()
    ax2.plot(step_log, smooth_data(zmin, window=1), color="tab:green", label="c_zmin", linewidth=1.5)
    ax2.set_ylabel("c_zmin (Å)", color="tab:green")

    lines1, labels1 = ax1.get_legend_handles_labels()
    lines2, labels2 = ax2.get_legend_handles_labels()
    ax1.legend(lines1 + lines2, labels1 + labels2, loc=2)
    ax1.grid(True, alpha=0.3)

    # === Bottom panel: Neighbor distances ===
    for layer, means in sorted(means_by_layer.items()):
        atom_id = first_id_by_layer.get(layer, "?")
        if layer <=3 or layer >= len(means_by_layer)-4:
            ax3.plot(steps_dump, means, label=f"layer {layer} (id {atom_id})",
                 linewidth=1.5, alpha=0.8)

    ax3.set_xlabel("Step")
    ax3.set_ylabel("Mean distance (Å)")
    ax3.set_title("Mean Neighbor Distance per Layer")
    ax3.legend(loc=2)
    ax3.grid(True, alpha=0.3)
    ax3.set_title(f"{args.dump}-distance")
    ax3.set_ylim(2.5, 4.5)
    ax4 = ax3.twinx()
    ax4.plot(step_log, smooth_data(zmin, window=1), color="tab:green", label="c_zmin", linewidth=1.5)
    ax4.set_ylabel("c_zmin (Å)", color="tab:green")
    
    # === Bottom panel: Last frame visualization ===
    img = mpimg.imread(frame_img_path)
    ax5.imshow(img, aspect='auto')
    ax5.axis('off')
    ax5.set_title(f"Last Frame from {dump_path.name}")
    ax5.margins(0)
    
    fig.tight_layout()
    fig.savefig("combined_plot.png", dpi=200)
    print("Saved combined_plot.png")
    
    # Clean up temporary file
    try:
        os.unlink(frame_img_path)
    except:
        pass


if __name__ == "__main__":
    main()
