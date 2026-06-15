#!/usr/bin/env python3
"""
Taxonomic Sankey Diagram Generator with Enhanced Color Options
Creates interactive sankey diagrams from taxonomic classification data.
"""

import argparse
import pandas as pd
import plotly.graph_objects as go
from collections import defaultdict
import sys
import os
import numpy as np
from typing import List, Dict, Tuple, Optional
import colorsys
import math

# 预定义颜色调色板 - 完全独立版本，不依赖Plotly
COLOR_PALETTES = {
    # 分类调色板 (适用于分类数据)
    'set2': [
        '#66c2a5', '#fc8d62', '#8da0cb', '#e78ac3',
        '#a6d854', '#ffd92f', '#e5c494', '#b3b3b3'
    ],
    'set3': [
        '#8dd3c7', '#ffffb3', '#bebada', '#fb8072',
        '#80b1d3', '#fdb462', '#b3de69', '#fccde5',
        '#d9d9d9', '#bc80bd', '#ccebc5', '#ffed6f'
    ],
    'pastel1': [
        '#fbb4ae', '#b3cde3', '#ccebc5', '#decbe4',
        '#fed9a6', '#ffffcc', '#e5d8bd', '#fddaec',
        '#f2f2f2'
    ],
    'pastel2': [
        '#b3e2cd', '#fdcdac', '#cbd5e8', '#f4cae4',
        '#e6f5c9', '#fff2ae', '#f1e2cc', '#cccccc'
    ],
    'dark2': [
        '#1b9e77', '#d95f02', '#7570b3', '#e7298a',
        '#66a61e', '#e6ab02', '#a6761d', '#666666'
    ],
    'bold': [
        '#7F3C8D', '#11A579', '#3969AC', '#F2B701',
        '#E73F74', '#80BA5A', '#E68310', '#008695',
        '#CF1C90', '#f97b72', '#4b4b8f', '#A5AA99'
    ],
    'vivid': [
        '#E58606', '#5D69B1', '#52BCA3', '#99C945',
        '#CC61B0', '#24796C', '#DAA51B', '#2F8AC4',
        '#764E9F', '#ED645A', '#CC3A8E', '#A5AA99'
    ],
    'prism': [
        '#5F4690', '#1D6996', '#38A6A5', '#0F8554',
        '#73AF48', '#EDAD08', '#E17C05', '#CC503E',
        '#94346E', '#6F4070', '#994E95', '#666666'
    ],
    'paired': [
        '#A6CEE3', '#1F78B4', '#B2DF8A', '#33A02C',
        '#FB9A99', '#E31A1C', '#FDBF6F', '#FF7F00',
        '#CAB2D6', '#6A3D9A', '#FFFF99', '#B15928'
    ],
    'set1': [
        '#E41A1C', '#377EB8', '#4DAF4A', '#984EA3',
        '#FF7F00', '#FFFF33', '#A65628', '#F781BF',
        '#999999'
    ],
    'tab10': [
        '#1f77b4', '#ff7f0e', '#2ca02c', '#d62728',
        '#9467bd', '#8c564b', '#e377c2', '#7f7f7f',
        '#bcbd22', '#17becf'
    ],
    'tab20': [
        '#1f77b4', '#aec7e8', '#ff7f0e', '#ffbb78',
        '#2ca02c', '#98df8a', '#d62728', '#ff9896',
        '#9467bd', '#c5b0d5', '#8c564b', '#c49c94',
        '#e377c2', '#f7b6d2', '#7f7f7f', '#c7c7c7',
        '#bcbd22', '#dbdb8d', '#17becf', '#9edae5'
    ],
    
    # 连续调色板 (适用于连续数据)
    'viridis': [
        '#440154', '#482878', '#3e4989', '#31688e',
        '#26828e', '#1f9e89', '#35b779', '#6ece58',
        '#b5de2b', '#fde725'
    ],
    'plasma': [
        '#0d0887', '#46039f', '#7201a8', '#9c179e',
        '#bd3786', '#d8576b', '#ed7953', '#fb9f3a',
        '#fdca26', '#f0f921'
    ],
    'inferno': [
        '#000004', '#1b0c41', '#4a0c6b', '#781c6d',
        '#a52c60', '#cf4446', '#ed6925', '#fb9b06',
        '#f7d03c', '#fcffa4'
    ],
    'magma': [
        '#000004', '#180f3d', '#440f76', '#721f81',
        '#9e2f7f', '#cd4071', '#f1605d', '#fd9668',
        '#feca8d', '#fcfdbf'
    ],
    'cividis': [
        '#00204d', '#00336f', '#39486b', '#575c6d',
        '#707173', '#8a8779', '#a69d78', '#c4b56a',
        '#e4cf5b', '#ffea46'
    ],
    'turbo': [
        '#30123b', '#4145ab', '#4675ed', '#39a2fc',
        '#1bcfd4', '#24eca6', '#61fc6c', '#a4fc3b',
        '#d1e834', '#f8c72a', '#ff9820', '#ff641e',
        '#f93a1a', '#d90b0c', '#a70403', '#7a0403'
    ],
    
    # 彩虹调色板
    'rainbow': [
        '#6e40aa', '#7a3fac', '#863ead', '#923dae',
        '#9e3cae', '#aa3bae', '#b63aad', '#c239ac',
        '#cd38ab', '#d837a9', '#e336a7', '#ee36a4',
        '#f836a1', '#ff389e', '#ff3c9b', '#ff4197',
        '#ff4793', '#ff4e8f', '#ff568b', '#ff5e86',
        '#ff6681', '#ff6f7c', '#ff7877', '#ff8272',
        '#ff8b6d', '#ff9567', '#ff9f62', '#ffa95d',
        '#ffb258', '#ffbc53', '#ffc64e', '#ffd049',
        '#ffda45', '#ffe440', '#ffee3c', '#f9f83a',
        '#f0ff38'
    ],
    'jet': [
        '#00008F', '#0000FF', '#007FFF', '#00FFFF',
        '#7FFF7F', '#FFFF00', '#FF7F00', '#FF0000',
        '#7F0000'
    ],
    
    # 生态/生物调色板
    'greens': [
        '#f7fcf5', '#e5f5e0', '#c7e9c0', '#a1d99b',
        '#74c476', '#41ab5d', '#238b45', '#006d2c',
        '#00441b'
    ],
    'blues': [
        '#f7fbff', '#deebf7', '#c6dbef', '#9ecae1',
        '#6baed6', '#4292c6', '#2171b5', '#08519c',
        '#08306b'
    ],
    'reds': [
        '#fff5f0', '#fee0d2', '#fcbba1', '#fc9272',
        '#fb6a4a', '#ef3b2c', '#cb181d', '#a50f15',
        '#67000d'
    ],
    'purples': [
        '#fcfbfd', '#efedf5', '#dadaeb', '#bcbddc',
        '#9e9ac8', '#807dba', '#6a51a3', '#54278f',
        '#3f007d'
    ],
    'oranges': [
        '#fff5eb', '#fee6ce', '#fdd0a2', '#fdae6b',
        '#fd8d3c', '#f16913', '#d94801', '#a63603',
        '#7f2704'
    ],
    'greys': [
        '#ffffff', '#f0f0f0', '#d9d9d9', '#bdbdbd',
        '#969696', '#737373', '#525252', '#252525',
        '#000000'
    ],
    
    # 地球色调色板
    'earth': [
        '#8B4513', '#A0522D', '#D2691E', '#CD853F',
        '#F4A460', '#DEB887', '#D2B48C', '#BC8F8F',
        '#8B7355', '#6B8E23'
    ],
    'forest': [
        '#228B22', '#32CD32', '#3CB371', '#66CDAA',
        '#8FBC8F', '#98FB98', '#90EE90', '#00FA9A',
        '#00FF7F', '#2E8B57'
    ],
    'ocean': [
        '#000080', '#0000CD', '#1E90FF', '#00BFFF',
        '#87CEEB', '#87CEFA', '#B0E0E6', '#E0FFFF',
        '#AFEEEE', '#48D1CC'
    ],
    'warm': [
        '#FF6B6B', '#FF8E72', '#FFAA7A', '#FFC285',
        '#FFD699', '#FFE5AD', '#FFF0C1', '#FFFAE6',
        '#E8F4EA', '#D4E9D7'
    ],
    'cool': [
        '#6A8EAE', '#7BA0C0', '#8CB2D3', '#9DC4E6',
        '#AED6F1', '#BFE9FC', '#D0ECFF', '#E1F5FE',
        '#F2F9FF', '#FFFFFF'
    ],
    'diverging': [
        '#003f5c', '#2f4b7c', '#665191', '#a05195',
        '#d45087', '#f95d6a', '#ff7c43', '#ffa600',
        '#ffff00', '#ffffff'
    ],
    
    # 分类学调色板 - 为分类学数据特别设计
    'taxonomy': [
        '#1f77b4',  # 蓝色 - Realm
        '#ff7f0e',  # 橙色 - Kingdom
        '#2ca02c',  # 绿色 - Phylum
        '#d62728',  # 红色 - Class
        '#9467bd',  # 紫色 - Order
        '#8c564b',  # 棕色 - Family
        '#e377c2',  # 粉色 - Genus
        '#7f7f7f',  # 灰色 - Species
        '#bcbd22',  # 黄绿色
        '#17becf',  # 青色
    ],
    'taxonomy2': [
        '#3d5a80',  # 深蓝
        '#98c1d9',  # 浅蓝
        '#e0fbfc',  # 淡蓝
        '#ee6c4d',  # 珊瑚
        '#a8dadc',  # 青绿
        '#457b9d',  # 中蓝
        '#e63946',  # 红
        '#f1faee',  # 米白
        '#ffafcc',  # 粉红
        '#a2d2ff',  # 淡蓝
    ],
    'nature': [
        '#3d5a80',  # 深蓝
        '#98c1d9',  # 浅蓝
        '#e0fbfc',  # 淡蓝
        '#ee6c4d',  # 珊瑚
        '#293241',  # 深灰蓝
        '#a8dadc',  # 青绿
        '#457b9d',  # 中蓝
        '#1d3557',  # 海军蓝
        '#e63946',  # 红
        '#f1faee',  # 米白
    ],
    
    # 病毒学专用调色板
    'virology': [
        '#4B0082',  # 靛蓝 - DNA病毒
        '#8A2BE2',  # 蓝紫 - 逆转录病毒
        '#00CED1',  # 深青 - RNA病毒
        '#20B2AA',  # 浅海绿 - 单链病毒
        '#32CD32',  # 石灰绿 - 双链病毒
        '#FFD700',  # 金色 - 有包膜病毒
        '#FF8C00',  # 深橙 - 无包膜病毒
        '#DC143C',  # 猩红 - 动物病毒
        '#228B22',  # 森林绿 - 植物病毒
        '#1E90FF',  # 道奇蓝 - 细菌病毒
    ],
}

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate interactive sankey diagrams from taxonomic classification data',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Color Palette Examples:
  # Use Set3 color palette (default)
  python taxonomic_sankey.py -i input.tsv --palette set3
  
  # Use Viridis sequential palette
  python taxonomic_sankey.py -i input.tsv --palette viridis
  
  # Use taxonomy-specific palette
  python taxonomic_sankey.py -i input.tsv --palette taxonomy
  
  # Use warm colors
  python taxonomic_sankey.py -i input.tsv --palette warm
  
  # Custom node and link colors
  python taxonomic_sankey.py -i input.tsv --node-color "#2E86AB" --link-color "#A23B72"
  
  # Show available color palettes
  python taxonomic_sankey.py --list-palettes
  
Color Scheme Options:
  level_based:    Each taxonomic level gets a unique color (default)
  flow_based:     Colors based on flow direction
  hierarchical:   Colors reflect hierarchical relationships
  random:         Random colors for each taxon
  uniform:        All nodes same color (with --node-color)
  
Examples:
  python taxonomic_sankey.py -i input.tsv --palette set3 --color-scheme level_based
  python taxonomic_sankey.py -i input.tsv --palette viridis
  python taxonomic_sankey.py -i input.tsv --palette virology --levels "Realm,Kingdom,Phylum,Class"
        """
    )
    
    # Required arguments
    parser.add_argument('-i', '--input', required=False,
                       help='Input TSV file with taxonomic classification')
    
    # Output arguments
    parser.add_argument('-o', '--output', default=None,
                       help='Output file name')
    parser.add_argument('--format', choices=['html', 'png', 'jpg', 'pdf', 'svg'],
                       default='html', help='Output format')
    
    # Layout arguments
    parser.add_argument('--width', type=int, default=1200,
                       help='Figure width in pixels')
    parser.add_argument('--height', type=int, default=1000,
                       help='Figure height in pixels')
    parser.add_argument('--title', default='Taxonomic Classification Sankey Diagram',
                       help='Plot title')
    parser.add_argument('--font-size', type=int, default=10,
                       help='Font size for labels')
    parser.add_argument('--font-family', default='Arial, sans-serif',
                       help='Font family for text')
    
    # Data arguments
    parser.add_argument('--levels', default='Realm,Kingdom,Phylum,Class,Order,Family,Genus',
                       help='Comma-separated taxonomic levels')
    parser.add_argument('--missing-label', default='Unknown',
                       help='Label for missing values')
    parser.add_argument('--min-flow', type=int, default=1,
                       help='Minimum flow value to display (applies to links and genus nodes)')
    parser.add_argument('--min-genus-flow', type=int, default=None,
                       help='Minimum flow value specifically for genus nodes (overrides --min-flow for genus)')
    
    # Color arguments - Core
    parser.add_argument('--palette', default='set3',
                       help='Color palette name (default: set3). Use --list-palettes to see options')
    parser.add_argument('--color-scheme', 
                       choices=['level_based', 'flow_based', 'hierarchical', 'random', 'uniform'],
                       default='level_based',
                       help='Color scheme for nodes')
    
    # Color arguments - Advanced
    parser.add_argument('--node-color', default=None,
                       help='Uniform color for all nodes (overrides palette)')
    parser.add_argument('--link-color', default=None,
                       help='Uniform color for all links')
    parser.add_argument('--reverse-palette', action='store_true',
                       help='Reverse the color palette')
    parser.add_argument('--opacity', type=float, default=0.85,
                       help='Opacity for nodes and links (0.0 to 1.0)')
    parser.add_argument('--brightness', type=float, default=1.0,
                       help='Brightness multiplier for colors (0.0 to 2.0)')
    parser.add_argument('--saturation', type=float, default=1.0,
                       help='Saturation multiplier for colors (0.0 to 2.0)')
    
    # Color arguments - Level-specific
    parser.add_argument('--level-colors', default=None,
                       help='Comma-separated colors for each level (overrides palette)')
    parser.add_argument('--link-opacity', type=float, default=None,
                       help='Opacity for links specifically')
    
    # Display arguments
    parser.add_argument('--simple', action='store_true',
                       help='Use simplified version')
    parser.add_argument('--no-stats', action='store_true',
                       help='Do not print statistics')
    parser.add_argument('--no-show', action='store_true',
                       help='Do not display the plot')
    parser.add_argument('--verbose', action='store_true',
                       help='Print detailed progress')
    parser.add_argument('--list-palettes', action='store_true',
                       help='List available color palettes and exit')
    
    # Advanced visualization
    parser.add_argument('--node-thickness', type=int, default=20,
                       help='Node thickness (default: 20)')
    parser.add_argument('--node-pad', type=int, default=25,
                       help='Padding between nodes (default: 25)')
    parser.add_argument('--show-percentages', action='store_true',
                       help='Show percentages in node labels')
    parser.add_argument('--label-truncate', type=int, default=30,
                       help='Truncate labels to this length (default: 30)')
    parser.add_argument('--label-angle', type=int, default=0,
                       help='Label angle for genus nodes (0-90, default: 0)')
    parser.add_argument('--genus-label-font-size', type=int, default=None,
                       help='Font size specifically for genus labels')
    parser.add_argument('--title-font-size', type=int, default=18,
                       help='Font size for title (default: 18)')
    
    return parser.parse_args()

def list_color_palettes():
    """List all available color palettes with visual preview."""
    print("Available Color Palettes:")
    print("=" * 80)
    
    categories = {
        '分类调色板 (Categorical)': ['set2', 'set3', 'pastel1', 'pastel2', 'dark2', 
                                  'bold', 'vivid', 'prism', 'paired', 'set1',
                                  'tab10', 'tab20'],
        '连续调色板 (Sequential)': ['viridis', 'plasma', 'inferno', 'magma', 'cividis',
                                  'turbo', 'rainbow', 'jet', 'greens', 'blues', 
                                  'reds', 'purples', 'oranges', 'greys'],
        '主题调色板 (Thematic)': ['earth', 'forest', 'ocean', 'warm', 'cool', 'diverging'],
        '分类学调色板 (Taxonomy)': ['taxonomy', 'taxonomy2', 'nature', 'virology']
    }
    
    for category, palettes in categories.items():
        print(f"\n{category}:")
        for palette in palettes:
            if palette in COLOR_PALETTES:
                colors = COLOR_PALETTES[palette]
                # Display palette name and first 3 colors
                color_samples = []
                for i, color in enumerate(colors[:3]):
                    # Create ANSI color codes for terminal preview
                    if color.startswith('#'):
                        r = int(color[1:3], 16)
                        g = int(color[3:5], 16)
                        b = int(color[5:7], 16)
                        color_samples.append(f"\033[48;2;{r};{g};{b}m  \033[0m")
                    else:
                        color_samples.append("██")
                
                print(f"  {palette:12} - {''.join(color_samples)} ... ({len(colors)} colors)")
    
    print("\n" + "=" * 80)
    print("Usage examples:")
    print("  python taxonomic_sankey.py --list-palettes")
    print("  python taxonomic_sankey.py -i input.tsv --palette set3")
    print("  python taxonomic_sankey.py -i input.tsv --palette taxonomy --color-scheme level_based")
    print("  python taxonomic_sankey.py -i input.tsv --palette virology --levels \"Realm,Kingdom,Phylum,Class\"")
    print("  python taxonomic_sankey.py -i input.tsv --palette warm --brightness 1.2")
    sys.exit(0)

def validate_color_arguments(args):
    """Validate color-related arguments."""
    # Check palette exists
    if args.palette.lower() not in COLOR_PALETTES:
        print(f"Error: Palette '{args.palette}' not found.", file=sys.stderr)
        print(f"Available palettes: {', '.join(sorted(COLOR_PALETTES.keys()))}", file=sys.stderr)
        sys.exit(1)
    
    # Validate opacity
    if not 0.0 <= args.opacity <= 1.0:
        print(f"Error: Opacity must be between 0.0 and 1.0, got {args.opacity}", file=sys.stderr)
        sys.exit(1)
    
    # Validate brightness
    if not 0.0 <= args.brightness <= 2.0:
        print(f"Error: Brightness must be between 0.0 and 2.0, got {args.brightness}", file=sys.stderr)
        sys.exit(1)
    
    # Validate saturation
    if not 0.0 <= args.saturation <= 2.0:
        print(f"Error: Saturation must be between 0.0 and 2.0, got {args.saturation}", file=sys.stderr)
        sys.exit(1)
    
    # Validate node thickness
    if args.node_thickness < 1:
        print(f"Error: Node thickness must be at least 1, got {args.node_thickness}", file=sys.stderr)
        sys.exit(1)
    
    # Validate label angle
    if not 0 <= args.label_angle <= 90:
        print(f"Error: Label angle must be between 0 and 90, got {args.label_angle}", file=sys.stderr)
        sys.exit(1)
    
    # Parse level colors if provided
    if args.level_colors:
        try:
            level_colors = [c.strip() for c in args.level_colors.split(',')]
            # Basic validation of color format
            for color in level_colors:
                if not (color.startswith('#') and len(color) == 7):
                    print(f"Warning: Color '{color}' may not be in valid hex format (#RRGGBB)", file=sys.stderr)
            args.parsed_level_colors = level_colors
        except Exception as e:
            print(f"Error parsing level colors: {e}", file=sys.stderr)
            sys.exit(1)

def hex_to_rgb(hex_color):
    """Convert hex color to RGB tuple."""
    hex_color = hex_color.lstrip('#')
    if len(hex_color) == 6:
        return tuple(int(hex_color[i:i+2], 16) for i in (0, 2, 4))
    elif len(hex_color) == 3:
        return tuple(int(hex_color[i:i+1]*2, 16) for i in (0, 1, 2))
    else:
        return (128, 128, 128)  # Default gray

def rgb_to_hex(rgb):
    """Convert RGB tuple to hex color."""
    return '#{:02x}{:02x}{:02x}'.format(*rgb)

def adjust_color(color: str, brightness: float = 1.0, saturation: float = 1.0) -> str:
    """Adjust the brightness and saturation of a color."""
    if brightness == 1.0 and saturation == 1.0:
        return color
    
    # Convert hex to RGB
    if color.startswith('#'):
        r, g, b = hex_to_rgb(color)
    elif color.startswith('rgb('):
        # Parse rgb string
        parts = color.replace('rgb(', '').replace(')', '').split(',')
        r, g, b = map(int, map(str.strip, parts))
    elif color.startswith('rgba('):
        # Parse rgba string, ignore alpha
        parts = color.replace('rgba(', '').replace(')', '').split(',')
        r, g, b = map(int, map(str.strip, parts[:3]))
    else:
        return color  # Return as-is if not recognized
    
    # Convert RGB to HSV
    h, s, v = colorsys.rgb_to_hsv(r/255.0, g/255.0, b/255.0)
    
    # Adjust saturation and brightness (value in HSV)
    s = max(0.0, min(1.0, s * saturation))
    v = max(0.0, min(1.0, v * brightness))
    
    # Convert back to RGB
    r_new, g_new, b_new = colorsys.hsv_to_rgb(h, s, v)
    
    # Convert to 0-255 range
    r_new = int(r_new * 255)
    g_new = int(g_new * 255)
    b_new = int(b_new * 255)
    
    return rgb_to_hex((r_new, g_new, b_new))

def apply_opacity(color: str, opacity: float) -> str:
    """Apply opacity to a color."""
    if opacity >= 1.0:
        return color
    
    if color.startswith('#'):
        # Convert hex to rgba
        r, g, b = hex_to_rgb(color)
        return f'rgba({r}, {g}, {b}, {opacity})'
    elif color.startswith('rgb('):
        # Convert rgb to rgba
        parts = color.replace('rgb(', '').replace(')', '').split(',')
        r, g, b = map(int, map(str.strip, parts))
        return f'rgba({r}, {g}, {b}, {opacity})'
    elif color.startswith('rgba('):
        # Replace existing opacity
        parts = color.replace('rgba(', '').replace(')', '').split(',')
        r, g, b = map(int, map(str.strip, parts[:3]))
        return f'rgba({r}, {g}, {b}, {opacity})'
    else:
        # Default to gray
        return f'rgba(128, 128, 128, {opacity})'

def get_palette_colors(palette_name: str, n_colors: int, reverse: bool = False) -> List[str]:
    """Get colors from a palette."""
    palette = COLOR_PALETTES[palette_name.lower()]
    
    if n_colors <= len(palette):
        colors = palette[:n_colors]
    else:
        # Repeat colors if needed
        colors = []
        for i in range(n_colors):
            color_idx = i % len(palette)
            colors.append(palette[color_idx])
    
    if reverse:
        colors = list(reversed(colors))
    
    return colors

def create_taxon_color_mapping(df, taxonomic_levels, args):
    """Create color mapping for all unique taxa."""
    color_info = {
        'scheme': args.color_scheme,
        'palette': args.palette,
        'taxon_to_color': {}
    }
    
    # Get base colors based on scheme
    if args.color_scheme == 'level_based':
        # Different color for each taxonomic level
        n_levels = len(taxonomic_levels)
        
        if args.level_colors and hasattr(args, 'parsed_level_colors'):
            level_colors = args.parsed_level_colors
            # Ensure we have enough colors
            while len(level_colors) < n_levels:
                level_colors.append(level_colors[-1])  # Repeat last color
        else:
            level_colors = get_palette_colors(args.palette, n_levels, args.reverse_palette)
        
        # Assign colors to each level
        level_color_map = {}
        for i, level in enumerate(taxonomic_levels):
            base_color = level_colors[i % len(level_colors)]
            adjusted_color = adjust_color(base_color, args.brightness, args.saturation)
            final_color = apply_opacity(adjusted_color, args.opacity)
            level_color_map[level] = final_color
        
        # Map each taxon to its level's color
        for level in taxonomic_levels:
            level_color = level_color_map[level]
            for taxon in df[level].unique():
                color_info['taxon_to_color'][(level, taxon)] = level_color
        
        color_info['level_colors'] = level_color_map
        
    elif args.color_scheme == 'uniform':
        # Uniform color for all taxa
        if args.node_color:
            base_color = args.node_color
        else:
            base_color = get_palette_colors(args.palette, 1, args.reverse_palette)[0]
        
        adjusted_color = adjust_color(base_color, args.brightness, args.saturation)
        final_color = apply_opacity(adjusted_color, args.opacity)
        
        for level in taxonomic_levels:
            for taxon in df[level].unique():
                color_info['taxon_to_color'][(level, taxon)] = final_color
        
        color_info['uniform_color'] = final_color
        
    elif args.color_scheme == 'random':
        # Random but consistent colors for each unique taxon
        palette_colors = get_palette_colors(args.palette, 100, args.reverse_palette)
        
        # Create deterministic random mapping based on taxon name
        for level in taxonomic_levels:
            for taxon in df[level].unique():
                # Create a hash from the taxon name
                hash_val = hash(f"{level}_{taxon}") % len(palette_colors)
                base_color = palette_colors[hash_val]
                adjusted_color = adjust_color(base_color, args.brightness, args.saturation)
                final_color = apply_opacity(adjusted_color, args.opacity)
                color_info['taxon_to_color'][(level, taxon)] = final_color
        
        color_info['unique_taxa_colored'] = len(color_info['taxon_to_color'])
        
    elif args.color_scheme == 'hierarchical':
        # Colors based on hierarchical path
        palette_colors = get_palette_colors(args.palette, 100, args.reverse_palette)
        
        # Group by full taxonomic path
        path_to_color = {}
        for idx, row in df.iterrows():
            path = tuple(str(row[level]) for level in taxonomic_levels)
            if path not in path_to_color:
                # Assign a color based on hash of the path
                hash_val = hash(path) % len(palette_colors)
                base_color = palette_colors[hash_val]
                adjusted_color = adjust_color(base_color, args.brightness, args.saturation)
                final_color = apply_opacity(adjusted_color, args.opacity)
                path_to_color[path] = final_color
            
            # Assign this color to each taxon in the path
            for level, taxon in zip(taxonomic_levels, path):
                color_info['taxon_to_color'][(level, taxon)] = path_to_color[path]
        
        color_info['unique_paths'] = len(path_to_color)
        
    else:  # flow_based (default)
        # Colors based on flow (simplified: color by first non-unknown level)
        palette_colors = get_palette_colors(args.palette, len(taxonomic_levels), args.reverse_palette)
        
        for level_idx, level in enumerate(taxonomic_levels):
            base_color = palette_colors[level_idx % len(palette_colors)]
            adjusted_color = adjust_color(base_color, args.brightness, args.saturation)
            final_color = apply_opacity(adjusted_color, args.opacity)
            
            for taxon in df[level].unique():
                color_info['taxon_to_color'][(level, taxon)] = final_color
        
        color_info['level_colors'] = {level: color_info['taxon_to_color'][(level, list(df[level].unique())[0])] 
                                      for level in taxonomic_levels}
    
    return color_info

def generate_sankey_diagram(df, taxonomic_levels, args):
    """Generate sankey diagram with enhanced color options and improved layout."""
    if args.verbose:
        print(f"Generating sankey diagram with {args.color_scheme} color scheme...")
        print(f"Using palette: {args.palette}")
        print(f"Processing {len(df)} sequences across {len(taxonomic_levels)} taxonomic levels")
    
    # 检查taxonomic_levels是否为空
    if not taxonomic_levels:
        print("Error: No valid taxonomic levels found in the data.", file=sys.stderr)
        print(f"Requested levels: {args.levels}", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Create color mapping
    color_info = create_taxon_color_mapping(df, taxonomic_levels, args)
    
    # Build sankey data structure
    nodes = []
    links = []
    node_indices = {}
    
    # Calculate total counts for size scaling
    total_sequences = len(df)
    
    # 确定Genus级别
    genus_level = taxonomic_levels[-1] if taxonomic_levels else 'Genus'
    
    # 设置Genus级别的min-flow阈值
    genus_min_flow = args.min_genus_flow if args.min_genus_flow is not None else args.min_flow
    
    if args.verbose:
        print(f"Applying min-flow={args.min_flow} for links and higher-level nodes")
        print(f"Applying min-flow={genus_min_flow} for {genus_level} nodes")
    
    # Create nodes for each unique taxon at each level
    node_counter = 0
    
    for level_idx, level in enumerate(taxonomic_levels):
        # Count occurrences of each taxon at this level
        taxon_counts = df[level].value_counts().to_dict()
        
        for taxon, count in taxon_counts.items():
            # 确定是否过滤此节点
            # 对于Genus级别，使用特定的min-flow阈值
            # 对于高级别分类（前几个级别），不过滤，或者使用更低的阈值
            if level == genus_level:
                # Genus级别：应用min-flow过滤
                if count < genus_min_flow:
                    continue
            else:
                # 高级别分类：不过滤，或者可以使用更低的阈值
                # 这里我们设置高级别分类的最小阈值为1，即不过滤
                pass
            
            # 创建节点标签，根据长度进行截断
            if len(taxon) > args.label_truncate:
                display_taxon = taxon[:args.label_truncate-3] + "..."
            else:
                display_taxon = taxon
            
            if args.show_percentages:
                percentage = (count / total_sequences) * 100
                if percentage >= 1.0:
                    label = f"{display_taxon}\n{count} ({percentage:.1f}%)"
                else:
                    label = f"{display_taxon}\n{count} ({percentage:.2f}%)"
            else:
                label = f"{display_taxon}\n({count})"
            
            # Get color for this taxon
            color_key = (level, taxon)
            if color_key in color_info['taxon_to_color']:
                color = color_info['taxon_to_color'][color_key]
            else:
                # Default color if not found
                default_color = adjust_color('#CCCCCC', args.brightness, args.saturation)
                color = apply_opacity(default_color, args.opacity)
            
            # 为Genus级别设置不同的字体大小
            font_size = args.genus_label_font_size if (level == genus_level and args.genus_label_font_size) else args.font_size
            
            # Create node
            nodes.append({
                'id': node_counter,
                'name': taxon,
                'level': level,
                'label': label,
                'color': color,
                'count': count,
                'font_size': font_size
            })
            
            # Store index for linking
            node_indices[(level, taxon)] = node_counter
            node_counter += 1
    
    # 检查是否创建了任何节点
    if not nodes:
        print("Error: No nodes were created. This might be due to min-flow thresholds being too high.", file=sys.stderr)
        print(f"  Total sequences: {total_sequences}", file=sys.stderr)
        print(f"  Min-flow for genus: {genus_min_flow}", file=sys.stderr)
        print(f"  Min-flow for links: {args.min_flow}", file=sys.stderr)
        sys.exit(1)
    
    # Create links between levels
    for i in range(len(taxonomic_levels) - 1):
        level_from = taxonomic_levels[i]
        level_to = taxonomic_levels[i + 1]
        
        # Count transitions between taxa
        transition_counts = {}
        
        for idx, row in df.iterrows():
            from_taxon = row[level_from]
            to_taxon = row[level_to]
            
            # Skip if either taxon is unknown
            if (from_taxon == f'{args.missing_label}_{level_from}' or 
                to_taxon == f'{args.missing_label}_{level_to}'):
                continue
            
            key = (from_taxon, to_taxon)
            transition_counts[key] = transition_counts.get(key, 0) + 1
        
        # Create links for transitions above threshold
        for (from_taxon, to_taxon), count in transition_counts.items():
            if count >= args.min_flow:
                from_key = (level_from, from_taxon)
                to_key = (level_to, to_taxon)
                
                if from_key in node_indices and to_key in node_indices:
                    source_idx = node_indices[from_key]
                    target_idx = node_indices[to_key]
                    
                    # Get link color
                    if args.link_color:
                        base_color = args.link_color
                    else:
                        # Use source node color
                        base_color = nodes[source_idx]['color']
                        # Extract RGB values if it's rgba
                        if 'rgba(' in base_color:
                            # Convert rgba to hex for adjustment
                            parts = base_color.replace('rgba(', '').replace(')', '').split(',')
                            r, g, b = map(int, map(str.strip, parts[:3]))
                            base_color = rgb_to_hex((r, g, b))
                    
                    # Adjust and apply opacity
                    adjusted_color = adjust_color(base_color, args.brightness, args.saturation)
                    link_opacity = args.link_opacity if args.link_opacity is not None else args.opacity * 0.7
                    link_color = apply_opacity(adjusted_color, link_opacity)
                    
                    links.append({
                        'source': source_idx,
                        'target': target_idx,
                        'value': count,
                        'color': link_color
                    })
    
    # Prepare data for Plotly
    node_labels = [node['label'] for node in nodes]
    node_colors = [node['color'] for node in nodes]
    
    link_sources = [link['source'] for link in links]
    link_targets = [link['target'] for link in links]
    link_values = [link['value'] for link in links]
    link_colors = [link['color'] for link in links]
    
    if args.verbose:
        print(f"Created {len(nodes)} nodes and {len(links)} links")
        print(f"Genus level: {genus_level}")
        print(f"Genus nodes: {len([n for n in nodes if n['level'] == genus_level])}")
    
    # Create sankey diagram with improved layout
    sankey_trace = go.Sankey(
        arrangement='perpendicular',
        node=dict(
            pad=args.node_pad,
            thickness=args.node_thickness,
            line=dict(color="rgba(0,0,0,0.3)", width=0.8),
            label=node_labels,
            color=node_colors,
            hovertemplate='%{label}<extra></extra>',
            align='justify',
        ),
        link=dict(
            source=link_sources,
            target=link_targets,
            value=link_values,
            color=link_colors,
            hovertemplate='%{value} sequences<extra></extra>',
            line=dict(width=0.5, color='rgba(0,0,0,0.2)')
        )
    )
    
    # 如果有标签角度，调整布局
    if args.label_angle > 0:
        sankey_trace.node.labelangle = args.label_angle
    
    fig = go.Figure(data=[sankey_trace])
    
    return fig, len(nodes), len(links), color_info, genus_level

def configure_layout(fig, args, num_nodes, num_links, color_info, taxonomic_levels, genus_level):
    """Configure plot layout with color information and improved spacing."""
    # Add color scheme info to title
    title_suffix = f"<br><sup>Palette: {args.palette} | Scheme: {color_info['scheme']}"
    if args.reverse_palette:
        title_suffix += " | Reversed"
    if args.brightness != 1.0:
        title_suffix += f" | Brightness: {args.brightness}"
    if args.saturation != 1.0:
        title_suffix += f" | Saturation: {args.saturation}"
    title_suffix += "</sup>"

    # 根据分类级别数量调整顶部边距
    num_levels = len(taxonomic_levels)
    # 动态计算顶部边距：基础值 + 每个级别的高度
    top_margin = 50 + (num_levels * 5)  # 减少顶部边距，因为图例移到底部

    # 增加底部边距，为图例留出空间
    bottom_margin = 50 + (num_levels * 15)

    # 标题: 空标题时跳过 (避免 HTML 嵌入时重叠)
    if args.title.strip():
        fig.update_layout(
            title_text=args.title + title_suffix,
            title_font_size=args.title_font_size,
            title_x=0.5,
            title_y=0.98,
        )
        # 统计信息注解 - 有标题时放在标题下方
        fig.add_annotation(
            text=f"Total sequences: {args.total_sequences:,} | Nodes: {num_nodes:,} | Links: {num_links:,}",
            showarrow=False,
            xref="paper", yref="paper",
            x=0.5, y=1.08,
            font_size=12,
            font_color="#333333",
            bgcolor="rgba(245,245,245,0.8)",
            bordercolor="#dddddd",
            borderwidth=1,
            borderpad=6,
            opacity=0.9
        )
    else:
        top_margin = 20  # 无标题时缩小顶部边距
    fig.update_layout(
        font_size=args.font_size,
        font_family=args.font_family,
        height=args.height,
        width=args.width,
        margin=dict(t=top_margin, b=bottom_margin, l=120, r=80),
        paper_bgcolor='white',
        plot_bgcolor='white',
    )

    # Add level annotations with matching colors - 移动到图表最下方
    level_colors = []
    if color_info['scheme'] == 'level_based' and 'level_colors' in color_info:
        level_colors = [color_info['level_colors'][level] for level in taxonomic_levels]
    elif color_info['scheme'] == 'uniform' and 'uniform_color' in color_info:
        level_colors = [color_info['uniform_color']] * len(taxonomic_levels)
    else:
        # Generate colors for annotations
        palette_colors = get_palette_colors(args.palette, len(taxonomic_levels), args.reverse_palette)
        level_colors = []
        for i, color in enumerate(palette_colors):
            adjusted = adjust_color(color, args.brightness, args.saturation)
            level_colors.append(apply_opacity(adjusted, 0.9))

    # 为每个分类级别添加标签 - 移动到图表最下方
    for i, level in enumerate(taxonomic_levels):
        if i < len(level_colors):
            color = level_colors[i]
        else:
            color = level_colors[-1] if level_colors else '#CCCCCC'

        # 计算x位置 - 更均匀分布
        if len(taxonomic_levels) > 1:
            x_pos = 0.12 + (i * 0.76 / (len(taxonomic_levels) - 1))
        else:
            x_pos = 0.5

        # 突出显示Genus级别
        font_weight = "bold" if level == genus_level else "normal"
        font_size = 14 if level == genus_level else 12

        # 将图例移到图表最下方
        fig.add_annotation(
            x=x_pos,
            y=-0.09,  # 移到图表底部下方
            xref="paper",
            yref="paper",
            text=f"<b>{level}</b>",
            showarrow=False,
            font=dict(size=font_size, color="#1a1a1a", family=args.font_family, weight=font_weight),
            bgcolor=color,
            bordercolor="rgba(0,0,0,0.5)",
            borderwidth=2,
            borderpad=10,
            opacity=0.95,
            align="center"
        )

    # 在底部添加图例标题
    fig.add_annotation(
        text="<b>Taxonomic Levels Legend</b>",
        showarrow=False,
        xref="paper", yref="paper",
        x=0.5, y=-0.10,  # 放在图例下方
        font_size=14,
        font_color="#333333",
        font_weight="bold"
    )

    # 统计信息注解放到 Legend 下方 (y=-0.16)
    fig.add_annotation(
        text=f"Total sequences: {args.total_sequences:,} | Nodes: {num_nodes:,} | Links: {num_links:,}",
        showarrow=False,
        xref="paper", yref="paper",
        x=0.5, y=-0.16,
        font_size=11,
        font_color="#555555",
    )

    return fig

def print_color_info(color_info, args):
    """Print color configuration information."""
    print("\n" + "="*60)
    print("COLOR CONFIGURATION")
    print("="*60)
    
    print(f"Color scheme: {color_info['scheme']}")
    print(f"Palette: {args.palette}")
    print(f"Opacity: {args.opacity}")
    print(f"Brightness: {args.brightness}")
    print(f"Saturation: {args.saturation}")
    
    if args.reverse_palette:
        print("Palette direction: Reversed")
    
    if color_info['scheme'] == 'level_based' and 'level_colors' in color_info:
        print("\nLevel colors:")
        for level, color in color_info['level_colors'].items():
            print(f"  {level:15}: {color}")
    
    if color_info['scheme'] == 'random' and 'unique_taxa_colored' in color_info:
        print(f"\nUnique taxa colored: {color_info['unique_taxa_colored']:,}")
    
    if color_info['scheme'] == 'hierarchical' and 'unique_paths' in color_info:
        print(f"\nUnique taxonomic paths: {color_info['unique_paths']:,}")
    
    if color_info['scheme'] == 'uniform' and 'uniform_color' in color_info:
        print(f"\nUniform color: {color_info['uniform_color']}")
    
    if args.link_color:
        print(f"\nLink color: {args.link_color}")
    elif args.link_opacity is not None:
        print(f"\nLink opacity: {args.link_opacity}")

def validate_input_file(input_file):
    """Validate that input file exists and is readable."""
    if not os.path.exists(input_file):
        print(f"Error: Input file '{input_file}' does not exist", file=sys.stderr)
        sys.exit(1)
    
    if not os.access(input_file, os.R_OK):
        print(f"Error: Input file '{input_file}' is not readable", file=sys.stderr)
        sys.exit(1)

def load_data(input_file, verbose=False):
    """Load and validate the input data."""
    try:
        df = pd.read_csv(input_file, sep='\t')
        
        if verbose:
            print(f"Loaded {len(df)} records from {input_file}")
            print(f"Columns found: {', '.join(df.columns)}")
        
        return df
    except Exception as e:
        print(f"Error loading input file: {e}", file=sys.stderr)
        sys.exit(1)

def prepare_taxonomic_levels(df, levels_str, missing_label='Unknown', verbose=False):
    """Prepare taxonomic levels and handle missing data."""
    taxonomic_levels = [level.strip() for level in levels_str.split(',')]
    
    # Validate that all levels exist in the dataframe
    missing_levels = [level for level in taxonomic_levels if level not in df.columns]
    if missing_levels:
        print(f"Warning: The following taxonomic levels are not in the input file: {', '.join(missing_levels)}", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        taxonomic_levels = [level for level in taxonomic_levels if level in df.columns]
    
    if verbose:
        print(f"Using taxonomic levels: {', '.join(taxonomic_levels)}")
    
    # Handle missing values
    for level in taxonomic_levels:
        if level in df.columns:
            df[level] = df[level].fillna(f'{missing_label}_{level}')
        else:
            print(f"Error: Level '{level}' not found in input data", file=sys.stderr)
            sys.exit(1)
    
    return df, taxonomic_levels

def save_figure(fig, output_file, output_format):
    """Save figure in specified format."""
    try:
        if output_format == 'html':
            fig.write_html(output_file)
            print(f"Interactive plot saved to: {output_file}")
        else:
            try:
                import kaleido
                fig.write_image(output_file, engine='kaleido')
                print(f"Static plot saved to: {output_file}")
            except ImportError:
                print("Error: For PNG/JPG/PDF/SVG output, please install kaleido:", file=sys.stderr)
                print("  pip install kaleido", file=sys.stderr)
                sys.exit(1)
    except Exception as e:
        print(f"Error saving figure: {e}", file=sys.stderr)
        sys.exit(1)

def print_statistics(df, taxonomic_levels):
    """Print detailed statistics to console."""
    print("\n" + "="*60)
    print("TAXONOMIC CLASSIFICATION STATISTICS")
    print("="*60)
    
    print(f"\nTotal sequences: {df.shape[0]:,}")
    print(f"Taxonomic levels analyzed: {len(taxonomic_levels)}")
    
    for level in taxonomic_levels:
        unique_count = df[level].nunique()
        print(f"\n{level}:")
        print(f"  Unique taxa: {unique_count:,}")
        
        # Print most common taxa
        counts = df[level].value_counts()
        total = counts.sum()
        
        if unique_count <= 10:
            print(f"  Taxon distribution:")
        else:
            print(f"  Top 10 taxa:")
        
        for taxon, count in counts.head(10).items():
            percentage = (count / total) * 100
            print(f"    - {taxon}: {count:,} sequences ({percentage:.1f}%)")
        
        if unique_count > 10:
            others = counts.iloc[10:].sum()
            if others > 0:
                print(f"    - Others ({unique_count-10} taxa): {others:,} sequences ({(others/total)*100:.1f}%)")

def main():
    """Main function."""
    args = parse_arguments()
    
    # Handle list-palettes option
    if args.list_palettes:
        list_color_palettes()
    
    # Validate input file
    if not args.input:
        print("Error: Input file is required. Use -i/--input to specify.", file=sys.stderr)
        sys.exit(1)
    
    validate_input_file(args.input)
    
    # Validate color arguments
    validate_color_arguments(args)
    
    # Load data
    df = load_data(args.input, args.verbose)
    args.total_sequences = len(df)
    
    # Prepare taxonomic levels
    df, taxonomic_levels = prepare_taxonomic_levels(
        df, args.levels, args.missing_label, args.verbose
    )
    
    # 检查是否有有效的分类级别
    if not taxonomic_levels:
        print("Error: No valid taxonomic levels found.", file=sys.stderr)
        print(f"Requested levels: {args.levels}", file=sys.stderr)
        print(f"Available columns: {', '.join(df.columns)}", file=sys.stderr)
        sys.exit(1)
    
    # Generate sankey diagram with colors
    fig, num_nodes, num_links, color_info, genus_level = generate_sankey_diagram(df, taxonomic_levels, args)
    
    # Configure layout - 传入taxonomic_levels参数
    fig = configure_layout(fig, args, num_nodes, num_links, color_info, taxonomic_levels, genus_level)
    
    # Determine output filename
    if args.output is None:
        base_name = os.path.splitext(args.input)[0]
        palette_suffix = f"_{args.palette}" if args.palette != 'set3' else ""
        scheme_suffix = f"_{args.color_scheme}" if args.color_scheme != 'level_based' else ""
        minflow_suffix = f"_minflow{args.min_flow}" if args.min_flow != 1 else ""
        genus_minflow_suffix = f"_genus{args.min_genus_flow}" if args.min_genus_flow else ""
        thickness_suffix = f"_thick{args.node_thickness}" if args.node_thickness != 20 else ""
        args.output = f"{base_name}_sankey{palette_suffix}{scheme_suffix}{minflow_suffix}{genus_minflow_suffix}{thickness_suffix}.{args.format}"
    
    # Save figure
    save_figure(fig, args.output, args.format)
    
    # Print statistics if not disabled
    if not args.no_stats:
        print_color_info(color_info, args)
        print_statistics(df, taxonomic_levels)
    
    # Display plot if not disabled
    if not args.no_show and args.format == 'html':
        try:
            fig.show()
        except Exception as e:
            print(f"Note: Could not display plot: {e}")
            print(f"Plot saved to: {args.output}")

if __name__ == "__main__":
    main()
