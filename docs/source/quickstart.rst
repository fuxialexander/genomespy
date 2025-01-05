Quickstart Guide
==============

This guide will help you get started with GenomeSpy quickly.

Basic Usage
----------

Here's a simple example of using GenomeSpy:

.. code-block:: python

    import pandas as pd

    from genomespy import GenomeSpy, igv

    tracks = {
        "ZBTB7A": {
            "url": "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX3161009.bw",
            "height": 40,
            "type": "bigwig"
        }}
    plot = igv(tracks, region={"chrom": "chr7", "start": 66600000, "end": 66800000}, server_port=18089)
    plot.show()

Note
----

If you are running this on a remote server, you need to port forward the server port to your local machine.
In vscode, you can do this by clicking the "Ports" button in the bottom panel and then clicking the "+" button to add a new port forward.