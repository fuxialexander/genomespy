Examples
========

Here are some examples of using GenomeSpy for various visualization tasks.

Basic Visualization
-----------------

.. code-block:: python

    import genomespy as gs
    import pandas as pd

    # Create sample data
    data = pd.DataFrame({
        'chromosome': ['chr1', 'chr1', 'chr2'],
        'start': [1000, 2000, 1500],
        'end': [1500, 2500, 2000],
        'value': [1.5, 2.0, 1.8]
    })

    # Create visualization
    vis = gs.Visualization()
    
    # Add track
    vis.add_track(data)
    
    # Show the visualization
    vis.show()

Advanced Examples
---------------

(Add more examples based on your package's features) 