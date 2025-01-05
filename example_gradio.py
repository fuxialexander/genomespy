from genomespy import igv, GenomeSpy
import gradio as gr

gs = GenomeSpy(server_port=18090)

def create_genomespy_plot(gs):
    """Create a GenomeSpy plot and return the HTML content."""
    tracks = {
        "ZBTB7A": {
            "url": "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX3161009.bw",
            "height": 40,
            "type": "bigwig"
        }
    }
    plot = igv(tracks, region={"chrom": "chr7", "start": 66600000, "end": 66800000}, gs=gs)
    plot.cleanup()
    return plot.show_gradio()

# Create a Gradio interface
iface = gr.Interface(
    fn=lambda: create_genomespy_plot(gs),
    inputs=[],
    outputs=gr.HTML(),
    live=True
)

iface.launch()
