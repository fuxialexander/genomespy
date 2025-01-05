# GenomeSpy

[![PyPI - Version](https://img.shields.io/pypi/v/genomespy.svg)](https://pypi.org/project/genomespy)
[![PyPI - Python Version](https://img.shields.io/pypi/pyversions/genomespy.svg)](https://pypi.org/project/genomespy)
[![Documentation Status](https://readthedocs.org/projects/genomespy/badge/?version=latest)](https://genomespy.readthedocs.io/en/latest/?badge=latest)

GenomeSpy is a Python wrapper for the JavaScript visualization grammar library [GenomeSpy](https://genomespy.app). It enables the creation of interactive, GPU-accelerated genomic visualizations through a declarative grammar.

The Python package provides a high-level interface to create GenomeSpy visualizations. Under the hood, it:

- Generates JSON specifications that GenomeSpy can interpret
- Creates a temporary HTTP server to serve local data files
- Handles rendering in Jupyter notebooks by creating temporary HTML files
- Supports exporting standalone HTML files for sharing (note: referenced data files must remain accessible)

The javascript (official) version of GenomeSpy is developed in The Systems Biology of Drug Resistance in Cancer group at the University of Helsinki by [Kari Lavikka](https://github.com/tuner). Original repo is at [genome-spy/genome-spy](https://github.com/genome-spy/genome-spy). 

> **Note**: This package has been primarily tested with the hg38 (GRCh38) human genome assembly. While it may work with other assemblies, full compatibility is not guaranteed. PRs are welcomed!

## Features

- 🚀 **GPU-accelerated rendering** for fluid interaction with large datasets
- 🧬 **Native support for genomic data formats** including BigWig, BigBed, and Indexed FASTA
- 📊 **Declarative visualization grammar** for mapping data to visual channels
- 🔍 **Interactive zooming and navigation** for exploring genomic regions
- 📱 **Responsive design** that adapts to different screen sizes
- 🔧 **Highly customizable** with extensive configuration options

## Installation

```console
pip install genomespy
```

## Quick Start

Here's a simple example to get you started:

```python
import genomespy as gs
import pandas as pd

# Create sample genomic data
data = pd.DataFrame({
    'chromosome': ['chr1', 'chr1', 'chr2'],
    'start': [1000, 2000, 1500],
    'end': [1500, 2500, 2000],
    'value': [1.5, 2.0, 1.8]
})

# Create visualization
vis = gs.GenomeSpy()

# Configure the visualization
vis.mark("rect", tooltip={"content": "data"})
vis.encode(
    x={"chrom": "chromosome", "pos": "start", "type": "locus"},
    x2={"chrom": "chromosome", "pos": "end"},
    y={"field": "value", "type": "quantitative"},
    color={"field": "value", "type": "quantitative"}
)

# Display the visualization
vis.show()
```

## Documentation

Comprehensive documentation is available at [genomespy.readthedocs.io](https://genomespy.readthedocs.io/), including:

- Detailed installation instructions
- API reference
- Usage examples
- Contributing guidelines

## Key Concepts

GenomeSpy uses a grammar of graphics approach where visualizations are built by:

1. **Marks**: Basic visual elements (rect, point, line, etc.)
2. **Encodings**: Mappings from data to visual properties
3. **Scales**: Functions that transform data values to visual values
4. **Transforms**: Data transformations and calculations
5. **Parameters**: Dynamic variables for interaction

## Development

To contribute to GenomeSpy:

```console
# Clone the repository
git clone https://github.com/fuxialexander/genomespy.git

# Install development dependencies
cd genomespy
pip install -e ".[dev]"

# Run tests
pytest
```

## License

`genomespy` is distributed under the terms of the [MIT](https://spdx.org/licenses/MIT.html) license.


## Support

- 📖 [Documentation](https://genomespy.readthedocs.io/)
- 🐛 [Issue Tracker](https://github.com/fuxialexander/genomespy/issues)
- 💬 [Discussions](https://github.com/fuxialexander/genomespy/discussions)
