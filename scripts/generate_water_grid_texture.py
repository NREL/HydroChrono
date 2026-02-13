#!/usr/bin/env python3
"""
Generate a simple grid texture for the water surface visualization.
Creates a seamless tileable grid pattern.
"""

import sys

try:
    from PIL import Image, ImageDraw
except ImportError:
    print("Pillow not installed. Installing...")
    import subprocess
    subprocess.check_call([sys.executable, "-m", "pip", "install", "Pillow"])
    from PIL import Image, ImageDraw


def create_grid_texture(size=256, line_width=3, num_cells=4):
    """
    Create a seamless tileable grid texture.
    
    Args:
        size: Texture size in pixels (square)
        line_width: Width of grid lines in pixels
        num_cells: Number of grid cells per texture tile
    
    Returns:
        PIL Image with grid pattern
    """
    # Water color (light blue-gray, similar to the PBR material)
    water_color = (20, 80, 120, 180)  # RGBA - semi-transparent blue
    
    # Grid line color (darker, more visible)
    line_color = (5, 30, 50, 220)  # RGBA - dark blue, mostly opaque
    
    # Create image with alpha channel
    img = Image.new('RGBA', (size, size), water_color)
    draw = ImageDraw.Draw(img)
    
    # Calculate cell size
    cell_size = size // num_cells
    
    # Draw vertical lines
    for i in range(num_cells + 1):
        x = i * cell_size
        draw.rectangle([x - line_width//2, 0, x + line_width//2, size - 1], fill=line_color)
    
    # Draw horizontal lines
    for i in range(num_cells + 1):
        y = i * cell_size
        draw.rectangle([0, y - line_width//2, size - 1, y + line_width//2], fill=line_color)
    
    return img


def main():
    # Generate the texture
    img = create_grid_texture(size=256, line_width=4, num_cells=4)
    
    # Save to data/textures/water_grid.png
    output_path = "data/textures/water_grid.png"
    img.save(output_path)
    print(f"Created grid texture: {output_path}")
    print(f"  Size: {img.size[0]}x{img.size[1]}")
    print(f"  Format: RGBA PNG")


if __name__ == "__main__":
    main()







