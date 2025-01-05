import pandas as pd
import numpy as np
import pytest
from genomespy import GenomeSpy, igv

def test_igv_visualization():
    """Test IGV-style visualization with bigwig track"""
    tracks = {
        "ZBTB7A": {
            "url": "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX3161009.bw",
            "height": 40,
            "type": "bigwig"
        }
    }
    region = {"chrom": "chr7", "start": 66600000, "end": 66800000}
    
    plot = igv(tracks, region=region)
    assert isinstance(plot, GenomeSpy)
    assert plot.height > 0
    assert "vconcat" in plot.spec
    plot.close()

def test_basic_scatter_plot():
    """Test basic scatter plot functionality"""
    plot = GenomeSpy()
    
    # Test data and visualization setup
    df = pd.DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})
    plot.data(df) \
        .mark('point') \
        .encode(
            x={'field': 'x', 'type': 'quantitative'}, 
            y={'field': 'y', 'type': 'quantitative'}
        ) \
        .transform([{'type': 'filter', 'expr': 'datum.x > 1'}]) \
        .scale(x={'domain': [0, 5]}) \
        .view({'name': 'main_view'}) \
        .expression('expr1', 'datum.x * 2') \
        .parameter('param1', 10)
    
    assert isinstance(plot, GenomeSpy)
    assert "data" in plot.spec
    assert "values" in plot.spec["data"]
    plot.close()

def test_clustermap():
    """Test clustermap visualization"""
    # Create sample data
    np.random.seed(42)  # For reproducibility
    data = pd.DataFrame(
        np.random.rand(10, 10),
        index=[f"Row{i}" for i in range(10)],
        columns=[f"Col{i}" for i in range(10)]
    )
    
    plot = GenomeSpy()
    plot.clustermap(
        data, 
        x_label="Columns", 
        y_label="Rows",
        z_score=1, 
        method="ward", 
        metric="euclidean",
        row_cluster=True, 
        col_cluster=True,
        vmax=1, 
        vmin=0, 
        center=0.5, 
        cmap="viridis"
    )
    
    assert isinstance(plot, GenomeSpy)
    assert "data" in plot.spec
    assert "values" in plot.spec["data"]
    plot.close()

def test_invalid_colormap():
    """Test that invalid colormap raises ValueError"""
    data = pd.DataFrame(np.random.rand(5, 5))
    plot = GenomeSpy()
    
    with pytest.raises(ValueError, match="Invalid color map"):
        plot.clustermap(data, cmap="invalid_colormap")
    plot.close()

def test_cleanup():
    """Test proper cleanup of resources"""
    plot = GenomeSpy()
    plot.data(pd.DataFrame({'x': [1, 2], 'y': [3, 4]}))
    plot.show()
    plot.close()
    
    # Verify that temporary files are cleaned up
    import os
    temp_files = [f for f in os.listdir() if f.startswith('.genomespy_temp_')]
    assert len(temp_files) == 0 