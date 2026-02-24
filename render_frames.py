#!/usr/bin/env python3
"""
Script to render LAMMPS dump files to PNG images using OVITO
"""

import sys
import os
from pathlib import Path
from ovito.io import import_file
from ovito.vis import Viewport, OSPRayRenderer


def render_dump_to_images(dump_file, output_dir=None, width=1920, height=1080, frame_step=1):
    """
    Load a LAMMPS dump file and render frames to PNG images
    
    Parameters:
    -----------
    dump_file : str
        Path to the dump file
    output_dir : str, optional
        Output directory for images. If None, creates 'frames' subdirectory
    width : int
        Image width in pixels
    height : int
        Image height in pixels
    frame_step : int
        Render every Nth frame (default: 1 = all frames)
    """
    
    # Import the dump file
    print(f"Loading dump file: {dump_file}")
    pipeline = import_file(dump_file)
    
    # Add pipeline to scene so it will be rendered
    pipeline.add_to_scene()
    
    # Generate output directory if not provided
    if output_dir is None:
        base_name = Path(dump_file).stem
        output_dir = f"{base_name}_frames"
    
    # Create output directory
    os.makedirs(output_dir, exist_ok=True)
    
    print(f"Output directory: {output_dir}")
    print(f"Number of frames: {pipeline.source.num_frames}")
    print(f"Rendering every {frame_step} frame(s)...")
    
    # Set up viewport for rendering - Left view (XZ plane)
    vp = Viewport()
    vp.type = Viewport.Type.Ortho  # Use orthographic projection
    vp.camera_pos = (0, 100, 0)  # Position camera along Y-axis
    vp.camera_dir = (0, -1, 0)   # Look towards the structure (negative Y)
    vp.fov = 60.0
    
    # Zoom to fit the scene
    vp.zoom_all()
    
    # Use OSPRay renderer (software-based, works without display)
    renderer = OSPRayRenderer()
    
    # Render frames
    frames_to_render = list(range(0, pipeline.source.num_frames, frame_step))
    for i, frame_num in enumerate(frames_to_render):
        output_file = os.path.join(output_dir, f"frame_{frame_num:05d}.png")
        
        try:
            vp.render_image(
                filename=output_file,
                size=(width, height),
                frame=frame_num,
                renderer=renderer,
                background=(1, 1, 1),  # White background
                alpha=False
            )
            
            if (i + 1) % 10 == 0 or i == len(frames_to_render) - 1:
                print(f"  Rendered {i + 1}/{len(frames_to_render)} frames")
                
        except Exception as e:
            print(f"  Error rendering frame {frame_num}: {e}")
            continue
    
    print(f"\nDone! Images saved to: {output_dir}/")
    return output_dir


def main():
    """Main function to handle command line arguments"""
    
    if len(sys.argv) < 2:
        print("Usage: python render_frames.py <dump_file> [output_dir] [width] [height] [frame_step]")
        print("\nExamples:")
        print("  python render_frames.py bend-L10-s50-b3-l1.dump")
        print("  python render_frames.py 10/bend-L10-s50-b3-l1.dump frames_output")
        print("  python render_frames.py dump.dump frames_output 1280 720 5")
        print("\nParameters:")
        print("  frame_step: Render every Nth frame (e.g., 5 = every 5th frame)")
        sys.exit(1)
    
    dump_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None
    width = int(sys.argv[3]) if len(sys.argv) > 3 else 1280
    height = int(sys.argv[4]) if len(sys.argv) > 4 else 720
    frame_step = int(sys.argv[5]) if len(sys.argv) > 5 else 1
    
    # Check if dump file exists
    if not Path(dump_file).exists():
        print(f"Error: Dump file not found: {dump_file}")
        sys.exit(1)
    
    render_dump_to_images(dump_file, output_dir, width, height, frame_step)


if __name__ == "__main__":
    main()
