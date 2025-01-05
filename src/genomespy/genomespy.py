import json
import os
from pathlib import Path
import re
import shutil
import webbrowser
from http.server import HTTPServer, SimpleHTTPRequestHandler
from threading import Thread
from typing import Any, Dict, List, Optional, Union

import numpy as np
import pandas as pd
from IPython.display import IFrame
from scipy.cluster.hierarchy import dendrogram, linkage

# Constants
DEFAULT_TRACK_HEIGHT = 50
DEFAULT_REGION = {
    "chrom": "chr7",
    "start": 66600000,
    "end": 66800000
}

class RangeRequestHandler(SimpleHTTPRequestHandler):
    """HTTP handler that supports range requests for bigwig/bigbed files.

    This handler extends the SimpleHTTPRequestHandler to support HTTP range requests,
    which are necessary for serving large genomic data files like bigwig and bigbed.

    Attributes:
        BINARY_EXTENSIONS (list): List of file extensions considered as binary.
    """
    
    BINARY_EXTENSIONS = ['.bw', '.bigwig', '.bb', '.bigbed']
    
    def log_message(self, format, *args):
        pass
    
    def guess_type(self, path):
        """Guess the type of a file based on its extension.

        Args:
            path (str): The file path.

        Returns:
            str: The MIME type of the file.
        """
        base, ext = os.path.splitext(path)
        if ext in self.BINARY_EXTENSIONS:
            return 'application/octet-stream'
        return super().guess_type(path)
    
    def _handle_range_request(self, f, file_size):
        """Handle range request for a file.

        Args:
            f (file object): The file object to read from.
            file_size (int): The total size of the file.

        Returns:
            file object or None: The file object positioned at the start of the requested range,
            or None if the range is invalid.
        """
        range_header = self.headers.get('Range')
        if not range_header:
            return None
            
        range_match = re.match(r'bytes=(\d+)-(\d*)', range_header)
        if not range_match:
            return None
            
        start = int(range_match.group(1))
        end = int(range_match.group(2)) if range_match.group(2) else file_size - 1
        
        if start >= file_size:
            self.send_error(416, "Requested range not satisfiable")
            f.close()
            return None
        
        self.send_response(206)
        self.send_header("Content-Range", f"bytes {start}-{end}/{file_size}")
        self.send_header("Content-Length", str(end - start + 1))
        self.send_header("Content-Type", self.guess_type(self.path))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        
        f.seek(start)
        return f
    
    def send_head(self):
        """Common code for GET and HEAD commands.

        Returns:
            file object or None: The file object to be sent to the client, or None if an error occurs.
        """
        path = self.translate_path(self.path)
        
        try:
            f = open(path, 'rb')
            file_size = os.fstat(f.fileno())[6]
        except:
            self.send_error(404, "File not found")
            return None
        
        # Handle range request if present
        range_response = self._handle_range_request(f, file_size)
        if range_response is not None:
            return range_response
        
        # Normal request (no range)
        self.send_response(200)
        self.send_header("Content-Length", str(file_size))
        self.send_header("Content-Type", self.guess_type(path))
        self.send_header("Accept-Ranges", "bytes")
        self.end_headers()
        return f
    
    def _handle_client_disconnect(self, func, *args, **kwargs):
        """Wrapper to handle client disconnections gracefully.

        Args:
            func (callable): The function to execute.
            *args: Variable length argument list for the function.
            **kwargs: Arbitrary keyword arguments for the function.

        Returns:
            Any: The result of the function call, or None if a disconnection occurs.
        """
        try:
            return func(*args, **kwargs)
        except (ConnectionResetError, BrokenPipeError):
            pass
    
    def handle(self):
        """Handle multiple requests if necessary."""
        self._handle_client_disconnect(super().handle)
    
    def handle_one_request(self):
        """Handle a single HTTP request."""
        self._handle_client_disconnect(super().handle_one_request)
    
    def copyfile(self, source, outputfile):
        """Copy all data between two file objects.

        Args:
            source (file object): The source file object.
            outputfile (file object): The destination file object.
        """
        self._handle_client_disconnect(super().copyfile, source, outputfile)

class GenomeSpy:
    """A Python wrapper for GenomeSpy visualization library.

    Parameters
    ----------
    height : int, optional
        The height of the visualization in pixels, by default 600

    Attributes
    ----------
    height : int
        The height of the visualization in pixels
    spec : dict
        The GenomeSpy specification defining the visualization structure
    _server_port : int
        The port number of the local HTTP server
    _template : str
        The HTML template for rendering the visualization

    Notes
    -----
    GenomeSpy is a toolkit for interactive visualization of genomic and other data. It enables 
    tailored visualizations through a declarative grammar inspired by Vega-Lite, allowing mapping 
    of data to visual channels (position, color, etc.) and composing complex visualizations from 
    primitive graphical marks (points, rectangles, etc.).

    Key Features:
    - GPU-accelerated rendering for fluid interaction with large datasets
    - Support for specialized genomic file formats (BigWig, BigBed, Indexed FASTA)
    - Built-in genomic coordinate handling and transformations
    - Interactive zooming and navigation
    - Composable visualization grammar
    """
    
    def __init__(self, height: int = 600):
        """Initialize a GenomeSpy instance.

        Args:
            height (int, optional): The height of the visualization. Defaults to 600.
        """
        self.height = height
        self.spec = {
            "$schema": "https://unpkg.com/@genome-spy/core/dist/schema.json",
            "data": {},
            "mark": {},
            "encoding": {},
            "transform": [],  # Initialize transform as an empty list
            "scales": {},
            "views": [],
            "parameters": {},
            "expressions": {},
        }
        self._server_port = None
        self._template = self._load_template()
    
    @staticmethod
    def _load_template():
        """Load the HTML template for visualization.

        Returns:
            str: The HTML template as a string.
        """
        return """
        <!DOCTYPE html>
        <html>
        <head>
            <title>GenomeSpy</title>
            <meta charset="UTF-8">
            <link rel="stylesheet" type="text/css" 
                  href="https://cdn.jsdelivr.net/npm/@genome-spy/app@0.51.x/dist/style.css" />
            <style>
                .genome-spy-container {{
                    width: 100%;
                    height: {height}px;
                    margin: 0 auto;
                    padding: 20px;
                    box-sizing: border-box;
                }}
            </style>
        </head>
        <body>
            <div class="genome-spy-container" id="visualization-container"></div>
            <script type="text/javascript" 
                    src="https://cdn.jsdelivr.net/npm/@genome-spy/app@0.51.x/dist/index.js">
            </script>
            <script>
                document.addEventListener('DOMContentLoaded', function() {{
                    setTimeout(() => {{
                        const container = document.getElementById('visualization-container');
                        genomeSpyApp.embed(container, {spec}, {{
                            defaultOptions: {{
                                width: "container",
                                height: "container"
                            }}
                        }});
                    }}, 100);
                }});
            </script>
        </body>
        </html>
        """

    def _start_server(self):
        """Start the local HTTP server."""
        # copy the shared directory to the current working directory
        shared_path = Path(__file__).parent / 'shared'
        dest_shared = Path.cwd() / '.genomespy_shared'
        shutil.copytree(shared_path, dest_shared, dirs_exist_ok=True)
        def server_thread():
            httpd = HTTPServer(('127.0.0.1', 0), RangeRequestHandler)
            self._server_port = httpd.server_port
            httpd.serve_forever()
        
        thread = Thread(target=server_thread)
        thread.daemon = True
        thread.start()
        
        # Give the server time to start
        import time
        time.sleep(0.1)

    def load_spec(self, spec: Union[str, Dict[str, Any]], is_url: bool = False):
        """Load a GenomeSpy specification.

        GenomeSpy specifications define how data should be visualized, including data sources, transformations,
        and visual encodings. Specifications can be loaded from a JSON file or directly as a dictionary.

        Args:
            spec (Union[str, Dict[str, Any]]): Either a JSON string/dict containing the spec or a URL to a spec file.
            is_url (bool, optional): Whether the spec is a URL to a JSON file. Defaults to False.

        Returns:
            GenomeSpy: The current instance for method chaining.
        """
        if is_url:
            self.spec = spec
        else:
            if isinstance(spec, str):
                try:
                    with open(spec, 'r') as f:
                        self.spec = json.load(f)
                        # Convert local file paths to server URLs for data files
                        self._process_local_data_files(self.spec)
                except FileNotFoundError:
                    raise FileNotFoundError(f"Could not find the file: {spec}")
                except json.JSONDecodeError:
                    raise ValueError(f"Invalid JSON format in file: {spec}")
                except Exception as e:
                    raise Exception(f"Error loading spec from file {spec}: {str(e)}")
            else:
                self.spec = spec
                self._process_local_data_files(self.spec)
        return self

    def _process_local_data_files(self, spec_obj):
        """Recursively process the spec to convert local file paths to server URLs.

        Args:
            spec_obj (dict or list): The specification object to process.
        """
        if isinstance(spec_obj, dict):
            # Handle data section
            if "data" in spec_obj and isinstance(spec_obj["data"], dict):
                if "lazy" in spec_obj["data"]:
                    lazy_data = spec_obj["data"]["lazy"]
                    if "url" in lazy_data and not lazy_data["url"].startswith(("http://", "https://")):
                        # Convert local file path to server URL
                        file_path = lazy_data["url"]
                        if os.path.exists(file_path):
                            lazy_data["url"] = f"http://127.0.0.1:{self._server_port}/{file_path}"

            # Recursively process all dictionary values
            for key, value in spec_obj.items():
                self._process_local_data_files(value)
        elif isinstance(spec_obj, list):
            # Recursively process all list items
            for item in spec_obj:
                self._process_local_data_files(item)

    def save_html(self, filename: str):
        """Save the visualization as a standalone HTML file.

        Args:
            filename (str): Output HTML file path.
        """
        spec_json = json.dumps(self.spec) if isinstance(self.spec, dict) else f'"{self.spec}"'
        html_content = self._template.format(height=self.height, spec=spec_json)
        
        with open(filename, 'w', encoding='utf-8') as f:
            f.write(html_content)

    def show(self, filename: Optional[str] = None):
        """Display the visualization in a browser or Jupyter notebook.

        Parameters
        ----------
        filename : str, optional
            Optional filename to save the HTML file. If None, creates a temporary file.

        Notes
        -----
        When running in a Jupyter notebook, the visualization will be displayed inline.
        Otherwise, it will open in the default web browser.

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> # Configure visualization...
        >>> plot.show()  # Display inline in notebook
        >>> 
        >>> # Save to specific file
        >>> plot.show("visualization.html")
        """
        if filename is None:
            filename = f'.genomespy_temp_{os.getpid()}.html'
        
        try:
            from IPython.display import display
            
            # Start local HTTP server in a separate thread
            def start_server():
                httpd = HTTPServer(('127.0.0.1', 0), RangeRequestHandler)
                self._server_port = httpd.server_port
                httpd.serve_forever()
            
            server_thread = Thread(target=start_server)
            server_thread.daemon = True
            server_thread.start()
            
            # Give the server a moment to start
            import time
            time.sleep(0.1)
            
            # Now process the spec and save the HTML
            self.save_html(filename)
            
            return display(IFrame(
                src=f'http://127.0.0.1:{self._server_port}/{os.path.basename(filename)}',
                width='100%',
                height=self.height + 40
            ))
            
        except ImportError:
            self.save_html(filename)
            webbrowser.open(f'file://{os.path.abspath(filename)}')

    def _repr_html_(self):
        """Jupyter notebook representation.

        Returns:
            str: The HTML representation of the visualization.
        """
        spec_json = json.dumps(self.spec) if isinstance(self.spec, dict) else f'"{self.spec}"'
        return self._template.format(height=self.height, spec=spec_json)
    
    def close(self):
        """Close the server if it's running and cleanup temporary files.

        Notes
        -----
        This method should be called when you're done with the visualization to:
        - Stop the local HTTP server if running
        - Remove any temporary files created during visualization
        - Free up system resources

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> # Create visualization...
        >>> plot.show()
        >>> plot.close()  # Cleanup when done
        """
        if hasattr(self, '_server_port') and self._server_port:
            self._server_port = None
        
        # Cleanup temporary files
        current_pid = os.getpid()
        temp_file = f'.genomespy_temp_{current_pid}.html'
        if os.path.exists(temp_file):
            try:
                os.remove(temp_file)
                # additionally remove any previous temp files
                for file in os.listdir():
                    if file.startswith('.genomespy_temp_'):
                        os.remove(file)
                # also remove the shared directory
                if os.path.exists('.genomespy_shared'):
                    shutil.rmtree('.genomespy_shared')
            except OSError:
                pass  # Ignore errors during cleanup

    def data(self, data: Union[pd.DataFrame, np.ndarray, str], format: str = "json"):
        """Set the data for the visualization.

        Parameters
        ----------
        data : Union[pd.DataFrame, np.ndarray, str]
            The data to visualize. Can be:
            - pandas DataFrame: Converted to records format
            - numpy array: Converted to list format
            - str: URL or path to data file
        format : str, optional
            The format of the data file if using URL/path, by default "json"
            Options include:
            - "json": JSON data
            - "csv": Comma-separated values
            - "tsv": Tab-separated values
            - "bigwig": BigWig genomic data
            - "bigbed": BigBed genomic data
            - "fasta": FASTA sequence data
            - "gff3": GFF3 genomic features

        Returns
        -------
        GenomeSpy
            The current instance for method chaining

        Notes
        -----
        GenomeSpy utilizes a tabular data structure as its fundamental data model, similar to a 
        spreadsheet or database table. Each dataset consists of records containing named data fields.

        Data Sources:
        - Eager data: Fully loaded during initialization (CSV, TSV, JSON)
        - Lazy data: Loaded on-demand (BigWig, BigBed, Indexed FASTA)
        - Named data: Can be dynamically updated using the API

        Examples
        --------
        >>> import pandas as pd
        >>> from genomespy import GenomeSpy
        >>> 
        >>> # Using pandas DataFrame
        >>> df = pd.DataFrame({'x': [1, 2, 3], 'y': [4, 5, 6]})
        >>> plot = GenomeSpy()
        >>> plot.data(df)
        >>> 
        >>> # Using file path
        >>> plot.data("data.bigwig", format="bigwig")
        """
        if isinstance(data, pd.DataFrame):
            self.spec["data"] = {"values": data.to_dict(orient="records")}
        elif isinstance(data, np.ndarray):
            self.spec["data"] = {"values": data.tolist()}
        elif isinstance(data, str):
            self.spec["data"] = {"url": data, "format": {"type": format}}
        return self

    def transform(self, transform: List[Dict[str, Any]]):
        """Add transformations to the visualization specification.

        Parameters
        ----------
        transform : List[Dict[str, Any]]
            A list of transformation specifications. Each transformation is a dictionary 
            with at least a "type" field and transformation-specific parameters.

        Returns
        -------
        GenomeSpy
            The current instance for method chaining

        Notes
        -----
        Transformations allow data manipulation before visualization. GenomeSpy provides 
        specialized transformations for genomic data visualization and analysis tasks.

        Common Transformations:
        - formula: Calculate new fields using expressions
        - filter: Filter data based on conditions
        - flatten: Flatten nested data structures
        - coverage: Calculate coverage from interval data
        - pileup: Create piled-up layout for overlapping features
        - flattenSequence: Split sequences into individual bases
        - collect: Group and sort data
        - project: Select and rename fields

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> plot.transform([
        ...     {
        ...         "type": "formula",
        ...         "expr": "datum.end - datum.start",
        ...         "as": "length"
        ...     },
        ...     {
        ...         "type": "filter",
        ...         "expr": "datum.length > 1000"
        ...     }
        ... ])
        """
        self.spec["transform"].extend(transform)
        return self

    def mark(self, mark_type: str, **kwargs):
        """Set the mark type for the visualization.

        Parameters
        ----------
        mark_type : str
            The type of mark to use
        **kwargs : dict
            Additional mark properties to configure appearance and behavior

        Returns
        -------
        GenomeSpy
            The current instance for method chaining

        Notes
        -----
        Marks are the basic graphical elements used to represent data. GenomeSpy provides 
        various mark types suitable for genomic data visualization.

        Mark Types:
        - rect: Rectangles (good for intervals, exons)
        - point: Points (good for variants, peaks)
        - line: Lines (good for continuous data)
        - rule: Rules (good for boundaries)
        - text: Text labels
        - area: Filled areas

        Mark Properties:
        - size: Size of the mark
        - color: Color of the mark
        - opacity: Transparency
        - strokeWidth: Width of stroke
        - tooltip: Tooltip configuration
        - minWidth: Minimum width for visibility
        - minOpacity: Minimum opacity for visibility

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> plot.mark("rect",
        ...     size=5,
        ...     minWidth=0.5,
        ...     tooltip={"content": "data"}
        ... )
        """
        self.spec["mark"] = {"type": mark_type, **kwargs}
        return self

    def encode(self, **kwargs):
        """Set the encoding for the visualization.

        Encodings map data fields to visual properties. GenomeSpy supports various encoding types and
        provides special support for genomic coordinates.

        Parameters
        ----------
        **kwargs : dict
            Encoding specifications for different channels.
            Each specification should be a dictionary defining the encoding properties.

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        Supported Channels
        ----------------
        - x, y: Position encoding
        - x2, y2: Secondary position for intervals
        - color: Color encoding
        - opacity: Transparency
        - size: Size of marks
        - text: Text content
        - tooltip: Tooltip content
        - sample: Sample ID for multi-sample visualizations

        Data Types
        ---------
        - quantitative: Numerical values
        - nominal: Categorical values
        - ordinal: Ordered categories
        - locus: Genomic coordinates (requires chrom and pos fields)

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> plot.encode(
        ...     x={"chrom": "chr", "pos": "start", "type": "locus"},
        ...     y={"field": "value", "type": "quantitative"},
        ...     color={"field": "category", "type": "nominal"}
        ... )
        """
        self.spec["encoding"] = kwargs
        return self

    def scale(self, **kwargs):
        """Set the scales for the visualization.

        Scales are functions that map abstract data values (e.g., a type of mutation) to visual values 
        (e.g., colors). GenomeSpy implements most of Vega-Lite's scale types and adds specialized scales 
        for genomic data.

        Parameters
        ----------
        **kwargs : dict
            Scale specifications for different channels. Each specification can include:
            - type: The type of scale to use
            - domain: Input domain range
            - range: Output range values
            - nice: Whether to extend domain to nice round numbers
            - padding: Padding to add around domain
            - scheme: Color scheme for color scales

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        Supported Scale Types
        -------------------
        - linear: Linear mapping for quantitative data
        - pow: Power scale for quantitative data
        - sqrt: Square root scale for quantitative data
        - symlog: Symmetric log scale
        - log: Logarithmic scale
        - ordinal: Discrete mapping for categorical data
        - band: Special scale for discrete ranges
        - point: Position-based scale
        - quantize: Binning for continuous data
        - threshold: Threshold-based binning

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> plot.scale(
        ...     y={
        ...         "type": "linear",
        ...         "domain": [0, 1],
        ...         "range": [0, 100],
        ...         "nice": True
        ...     },
        ...     color={
        ...         "type": "ordinal",
        ...         "domain": ["A", "C", "G", "T"],
        ...         "range": ["red", "blue", "green", "yellow"]
        ...     }
        ... )
        """
        self.spec["scales"].update(kwargs)
        return self

    def view(self, view_spec: Dict[str, Any]):
        """Add a view to the visualization.

        Views in GenomeSpy allow for hierarchical composition of visualizations. Views can be concatenated,
        layered, or arranged in other ways. Each view inherits data and encoding from its parent but can
        override them with its own specifications.

        Parameters
        ----------
        view_spec : Dict[str, Any]
            The view specification defining the visualization properties, data, marks, 
            and encodings for this view.

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        View Properties
        --------------
        - data : Data source for the view
        - transform : Data transformations
        - mark : Visual marks to represent data
        - encoding : Visual encodings
        - height : View height
        - width : View width
        - name : Unique identifier for the view
        - title : View title
        - description : View description
        - padding : Space around the view
        - opacity : View opacity
        - configurableVisibility : Whether view can be toggled

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> plot.view({
        ...     "name": "genes",
        ...     "height": 120,
        ...     "data": {"url": "genes.bed"},
        ...     "mark": "rect",
        ...     "encoding": {
        ...         "x": {"chrom": "chr", "pos": "start", "type": "locus"},
        ...         "x2": {"chrom": "chr", "pos": "end"}
        ...     }
        ... })
        """
        self.spec["views"].append(view_spec)
        return self

    def import_view(self, url: str):
        """Import a view from a URL.

        This function allows importing external view specifications, enabling reuse and sharing of
        visualization components. Common uses include importing standard genomic tracks like:
        - Chromosome ideograms
        - Gene annotation tracks
        - Reference genome sequences

        Parameters
        ----------
        url : str
            The URL or path to the view specification to import.
            Can be absolute URL or relative to the base URL.

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        Built-in Views
        -------------
        The following views are available in the .genomespy_shared/ directory:
        - cytobands.json : Chromosome ideogram track
        - genes.json : Gene annotation track
        - hg38.json : Reference genome sequence

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> # Import chromosome ideogram
        >>> plot.import_view(".genomespy_shared/cytobands.json")
        >>> 
        >>> # Import gene annotations
        >>> plot.import_view(".genomespy_shared/genes.json")
        >>> 
        >>> # Import reference genome
        >>> plot.import_view(".genomespy_shared/hg38.json")
        """
        self.spec["views"].append({"import": {"url": url}})
        return self

    def expression(self, name: str, expr: str):
        """Add an expression to the visualization.

        Expressions in GenomeSpy allow for computing new data fields or modifying existing ones.
        They use a JavaScript-like syntax and can access the current data object using 'datum'.
        Expressions can be used in transforms, encodings, and other places where dynamic 
        computation is needed.

        Parameters
        ----------
        name : str
            The name of the expression to be referenced elsewhere in the specification.
        expr : str
            The expression string using GenomeSpy's expression syntax.
            Can access current data object via 'datum'.

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        Common Uses
        ----------
        - Computing derived values
        - Conditional logic
        - String manipulation
        - Mathematical calculations
        - Accessing parameters

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> # Calculate length of genomic interval
        >>> plot.expression("length", "datum.end - datum.start")
        >>> 
        >>> # Compute log ratio
        >>> plot.expression("logRatio", "log2(datum.value / datum.control)")
        >>> 
        >>> # Create conditional label
        >>> plot.expression(
        ...     "label", 
        ...     "datum.score > 0.05 ? 'High impact' : 'Low impact'"
        ... )
        """
        self.spec["expressions"][name] = expr
        return self

    def parameter(self, name: str, value: Any):
        """Add a parameter to the visualization.

        Parameters enable dynamic behaviors and interactions in GenomeSpy visualizations.
        They can be used for interactive selections, conditional encoding, data filtering,
        and parameterizing imported specifications.

        Parameters
        ----------
        name : str
            The name of the parameter to be referenced in expressions and conditions.
        value : Any
            The parameter value or configuration. Can be a simple value
            or a parameter definition object.

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        Parameter Types
        --------------
        - Selection parameters : Enable interactive data selection
        - Value parameters : Store single values
        - Range parameters : Store numeric ranges
        - Vector parameters : Store arrays of values

        Common Uses
        ----------
        - Interactive filtering
        - Conditional encoding
        - Dynamic thresholds
        - Coordinated selections
        - View parameterization

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> # Selection parameter for interactive highlighting
        >>> plot.parameter("highlight", {
        ...     "select": {"type": "point", "on": "pointerover"}
        ... })
        >>> 
        >>> # Value parameter for filtering
        >>> plot.parameter("threshold", 0.05)
        >>> 
        >>> # Use in encoding
        >>> plot.encode(
        ...     opacity={
        ...         "condition": {"param": "highlight", "value": 1.0},
        ...         "value": 0.3
        ...     }
        ... )
        """
        self.spec["parameters"][name] = value
        return self

    def to_json(self):
        """Convert the specification to a JSON string.

        This function serializes the current GenomeSpy specification into a JSON string, which can be used for
        saving or sharing the visualization configuration.

        Returns
        -------
        str
            The JSON string representation of the specification.

        Examples
        --------
        >>> plot = GenomeSpy()
        >>> plot.encode(x={"field": "value", "type": "quantitative"})
        >>> json_spec = plot.to_json()
        """
        return json.dumps(self.spec, indent=2)

    def heatmap(self, data: pd.DataFrame, x_label: str = "x", y_label: str = "y"):
        """Create a heatmap from a pandas DataFrame.

        Heatmaps are a common way to visualize matrix-like data, where values are represented by colors. This
        function prepares the data and sets up the GenomeSpy specification for rendering a heatmap.

        Parameters
        ----------
        data : pd.DataFrame
            A pandas DataFrame containing the data for the heatmap.
        x_label : str, optional
            The label for the x-axis. Defaults to "x".
        y_label : str, optional
            The label for the y-axis. Defaults to "y".

        Returns
        -------
        GenomeSpy
            The current instance for method chaining.

        Examples
        --------
        >>> import pandas as pd
        >>> plot = GenomeSpy()
        >>> data = pd.DataFrame({
        ...     'A': [1, 2, 3],
        ...     'B': [4, 5, 6],
        ...     'C': [7, 8, 9]
        ... })
        >>> plot.heatmap(data, x_label="Samples", y_label="Features")
        """
        # Ensure the DataFrame has a name for the index
        if data.index.name is None:
            data.index.name = "index"

        # Melt the DataFrame to long format
        values = (
            data.reset_index().melt(id_vars=data.index.name).to_dict(orient="records")
        )

        # Define the GenomeSpy specification
        self.spec = {
            "$schema": "https://unpkg.com/@genome-spy/core/dist/schema.json",
            "data": {"values": values},
            "mark": {
                "type": "rect",
                "encoding": {
                    "x": {
                        "field": "variable",
                        "type": "nominal",
                        "axis": {"title": x_label},
                    },
                    "y": {
                        "field": data.index.name,
                        "type": "nominal",
                        "axis": {"title": y_label},
                    },
                    "color": {
                        "field": "value",
                        "type": "quantitative",
                        "scale": {
                            "scheme": "viridis",
                        },
                    },
                },
            },
        }
        return self

    def clustermap(
        self,
        data: pd.DataFrame,
        x_label: str = "x",
        y_label: str = "y",
        method: str = "ward",
        metric: str = "euclidean",
        z_score: Optional[int] = None,
        standard_scale: Optional[int] = None,
        row_cluster: bool = True,
        col_cluster: bool = True,
        vmax: Optional[float] = None,
        vmin: Optional[float] = None,
        center: Optional[float] = None,
        cmap: str = "viridis",
    ):
        """Create a clustermap from a pandas DataFrame.

        A clustermap combines a heatmap with hierarchical clustering dendrograms on both axes.
        The clustering helps reveal patterns and relationships in the data by grouping similar
        rows and columns together.

        Parameters
        ----------
        data : pd.DataFrame
            Input data matrix to be clustered and visualized
        x_label : str, optional
            Label for x-axis, by default "x"
        y_label : str, optional
            Label for y-axis, by default "y" 
        method : str, optional
            Linkage method for hierarchical clustering, by default "ward"
        metric : str, optional
            Distance metric for clustering, by default "euclidean"
        z_score : int, optional
            Standardize the data along rows (0) or columns (1), by default None
        standard_scale : int, optional
            Scale data along rows (0) or columns (1), by default None
        row_cluster : bool, optional
            Whether to cluster rows, by default True
        col_cluster : bool, optional
            Whether to cluster columns, by default True
        vmax : float, optional
            Maximum value for color scaling, by default None
        vmin : float, optional
            Minimum value for color scaling, by default None
        center : float, optional
            Center value for diverging colormaps, by default None
        cmap : str, optional
            Colormap name, either "viridis" or "blues", by default "viridis"

        Returns
        -------
        GenomeSpy
            The current instance for method chaining

        Examples
        --------
        >>> import pandas as pd
        >>> from genomespy import GenomeSpy
        >>> 
        >>> # Create sample data
        >>> data = pd.DataFrame({
        ...     'A': [1, 2, 3],
        ...     'B': [2, 4, 6],
        ...     'C': [3, 6, 9]
        ... })
        >>> 
        >>> # Create and display clustermap
        >>> plot = GenomeSpy()
        >>> plot.clustermap(
        ...     data,
        ...     x_label="Samples",
        ...     y_label="Features",
        ...     z_score=1,
        ...     method="ward"
        ... )
        """
        # Ensure the DataFrame has a name for the index
        if data.index.name is None:
            data.index.name = "index"

        if cmap not in ["viridis", "blues"]:
            raise ValueError("Invalid color map. Please use 'viridis' or 'blues'.")

        # Apply z-score normalization
        if z_score is not None:
            if z_score == 0:
                data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=1)
            elif z_score == 1:
                data = data.apply(lambda x: (x - x.mean()) / x.std(), axis=0)

        # Apply standard scaling
        if standard_scale is not None:
            if standard_scale == 0:
                data = data.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=1)
            elif standard_scale == 1:
                data = data.apply(lambda x: (x - x.min()) / (x.max() - x.min()), axis=0)

        # Perform hierarchical clustering on rows
        if row_cluster:
            Z_rows = linkage(data, method=method, metric=metric)
            row_dendro = dendrogram(
                Z_rows, labels=data.index, orientation="left", no_plot=True
            )
            data = data.iloc[row_dendro["leaves"], :]

        # Perform hierarchical clustering on columns
        if col_cluster:
            Z_cols = linkage(data.T, method=method, metric=metric)
            col_dendro = dendrogram(
                Z_cols, labels=data.columns, orientation="top", no_plot=True
            )
            data = data.iloc[:, col_dendro["leaves"]]

        # Define the GenomeSpy specification for the heatmap
        color_scale = {
            "scheme": cmap,
        }

        if vmin is not None:
            color_scale["domainMin"] = vmin
        if vmax is not None:
            color_scale["domainMax"] = vmax
        if center is not None:
            color_scale["domainMid"] = center

        self.spec = {
            "$schema": "https://unpkg.com/@genome-spy/core/dist/schema.json",
            "data": {
                "values": data.reset_index()
                .melt(id_vars=data.index.name)
                .to_dict(orient="records")
            },
            "mark": "rect",
            "encoding": {
                "x": {
                    "field": "variable",
                    "type": "nominal",
                    "axis": {"title": x_label},
                },
                "y": {
                    "field": data.index.name,
                    "type": "nominal",
                    "axis": {"title": y_label},
                },
                "color": {
                    "field": "value",
                    "type": "quantitative",
                    "scale": color_scale,
                },
            },
        }

        return self

    def dendrogram(
        self,
        data: pd.DataFrame,
        method: str = "ward",
        metric: str = "euclidean"
    ):
        """Create a dendrogram using GenomeSpy.

        Dendrograms are tree-like diagrams used to visualize the arrangement of clusters produced by hierarchical
        clustering.

        Parameters
        ----------
        data : pd.DataFrame
            Input data matrix for clustering
        method : str, optional
            Linkage method for clustering, by default "ward"
        metric : str, optional
            Distance metric for clustering, by default "euclidean"

        Returns
        -------
        GenomeSpy
            The current instance for method chaining

        Examples
        --------
        >>> import pandas as pd
        >>> plot = GenomeSpy()
        >>> data = pd.DataFrame({
        ...     'A': [1, 2, 3],
        ...     'B': [4, 5, 6]
        ... })
        >>> plot.dendrogram(data, method="ward", metric="euclidean")
        """
        # Perform hierarchical clustering
        from scipy.cluster.hierarchy import dendrogram, linkage

        Z = linkage(data, method=method, metric=metric)
        dendro = dendrogram(Z, no_plot=True)

        # Prepare data for GenomeSpy
        icoord = np.array(dendro["icoord"])
        dcoord = np.array(dendro["dcoord"])
        lines = []

        for xs, ys in zip(icoord, dcoord):
            for i in range(3):
                lines.append({"x": xs[i], "x2": xs[i + 1], "y": ys[i], "y2": ys[i + 1]})

        # Define the GenomeSpy specification
        self.spec = {
            "$schema": "https://unpkg.com/@genome-spy/core/dist/schema.json",
            "data": {"values": lines},
            "mark": {"type": "rule", "strokeDash": [6, 3]},
            "encoding": {
                "x": {"field": "x", "type": "quantitative"},
                "x2": {"field": "x2", "type": "quantitative"},
                "y": {"field": "y", "type": "quantitative"},
                "y2": {"field": "y2", "type": "quantitative"},
                "color": {"field": "y", "type": "nominal"},
            },
        }
        return self

# Additional helper functions and classes can be added here as needed. 

def _get_track_height(track_spec):
    """Helper function to get track height.

    Parameters
    ----------
    track_spec : dict
        The track specification

    Returns
    -------
    int
        The height of the track in pixels

    Notes
    -----
    Handles both numeric and dictionary height values. For tracks with step-based
    height (like Gencode), uses a fixed viewport height.
    """
    if "viewportHeight" in track_spec:
        return track_spec["viewportHeight"]
    
    height = track_spec.get("height", DEFAULT_TRACK_HEIGHT)
    if isinstance(height, dict):
        return 280  # Fixed height for Gencode track
    return height

def create_track_spec(track_name: str, track_config: Dict[str, Any], region: Dict[str, Any]) -> Dict[str, Any]:
    """Create a track specification for GenomeSpy.

    Parameters
    ----------
    track_name : str
        The name of the track
    track_config : Dict[str, Any]
        The configuration for the track
    region : Dict[str, Any]
        The genomic region for the track

    Returns
    -------
    Dict[str, Any]
        The complete track specification

    Examples
    --------
    >>> region = {"chrom": "chr1", "start": 1000, "end": 2000}
    >>> config = {
    ...     "type": "bigwig",
    ...     "url": "data.bw",
    ...     "height": 100
    ... }
    >>> spec = create_track_spec("Coverage", config, region)
    """
    height = track_config.get('height', DEFAULT_TRACK_HEIGHT)
    
    track_spec = {
        "height": height,
        "name": track_name,
        "view": {"stroke": "lightgray"},
        "data": {
            "lazy": {
                "type": track_config.get('type', 'bigwig'),
                "pixelsPerBin": 1
            }
        },
        "encoding": {
            "x": {
                "chrom": "chrom",
                "pos": "start",
                "type": "locus",
                "scale": {
                    "domain": [
                        {"chrom": region["chrom"], "pos": region["start"]},
                        {"chrom": region["chrom"], "pos": region["end"]}
                    ]
                }
            },
            "x2": {
                "chrom": "chrom",
                "pos": "end"
            },
            "y": {
                "field": "score",
                "type": "quantitative",
                "scale": {"nice": True},
                "axis": {
                    "title": track_name,
                    "grid": True,
                    "gridDash": [2, 2],
                    "maxExtent": 35
                }
            }
        },
        "mark": {
            "type": "rect",
            "minWidth": 0.5,
            "minOpacity": 1,
            "tooltip": None
        }
    }
    
    # Set the URL/path
    if 'url' in track_config:
        track_spec["data"]["lazy"]["url"] = track_config["url"]
    elif 'path' in track_config:
        track_spec["data"]["lazy"]["url"] = track_config["path"]
    else:
        raise ValueError(f"Track {track_name} must have either 'url' or 'path' specified")
    
    return track_spec

def create_base_spec(region: Dict[str, Any]) -> Dict[str, Any]:
    """Create the base specification for GenomeSpy visualization.

    Parameters
    ----------
    region : Dict[str, Any]
        The genomic region for the visualization

    Returns
    -------
    Dict[str, Any]
        The base specification including schema and default tracks

    Examples
    --------
    >>> region = {"chrom": "chr1", "start": 1000, "end": 2000}
    >>> base_spec = create_base_spec(region)
    """
    return {
        "$schema": "https://unpkg.com/@genome-spy/core/dist/schema.json",
        "genome": {"name": "hg38"},
        "resolve": {"axis": {"x": "shared"}},
        "vconcat": [
            {"import": {"url": ".genomespy_shared/cytobands.json"}},
            {"import": {"url": ".genomespy_shared/genes.json"}},
            {"import": {"url": ".genomespy_shared/hg38.json"}}
        ]
    }

def igv(file_dict: Dict[str, Dict[str, Any]], region: Optional[Dict[str, Any]] = None) -> GenomeSpy:
    """Create a GenomeSpy visualization with custom tracks in IGV style.

    This function creates a genome browser visualization similar to IGV (Integrative Genomics Viewer),
    with support for various genomic data formats and customizable tracks.

    Parameters
    ----------
    file_dict : Dict[str, Dict[str, Any]]
        A dictionary mapping track names to their configurations.
        Each track configuration should specify:
        - url or path : Path to the data file
        - type : Data format (e.g., "bigwig", "bigbed")
        - height : Track height in pixels
    region : Optional[Dict[str, Any]], optional
        The genomic region to display, by default None.
        Should contain:
        - chrom : Chromosome name
        - start : Start position
        - end : End position

    Returns
    -------
    GenomeSpy
        The configured GenomeSpy instance ready for display

    Examples
    --------
    >>> from genomespy import igv
    >>> # Configure tracks
    >>> tracks = {
    ...     "ZBTB7A": {
    ...         "url": "https://chip-atlas.dbcls.jp/data/hg38/eachData/bw/SRX3161009.bw",
    ...         "height": 40,
    ...         "type": "bigwig"
    ...     }
    ... }
    >>> # Create visualization
    >>> plot = igv(
    ...     tracks,
    ...     region={"chrom": "chr7", "start": 66600000, "end": 66800000}
    ... )
    >>> plot.show()
    """
    region = region or DEFAULT_REGION
    gs = GenomeSpy()
    
    # Start server
    gs._start_server()
    
    # Create base specification
    spec = create_base_spec(region)
    
    # Add custom tracks
    for track_name, track_config in file_dict.items():
        track_spec = create_track_spec(track_name, track_config, region)
        spec["vconcat"].append(track_spec)
    
    # Add standard tracks (cCRE, Gencode)
    spec["vconcat"].extend([
        create_ccre_track(region),
        create_gencode_track(region)
    ])
    
    # Update visualization height
    total_height = sum(_get_track_height(track) for track in spec["vconcat"]) + 100
    gs.height = total_height
    gs.spec = spec
    
    return gs

def create_ccre_track(region: Dict[str, Any]) -> Dict[str, Any]:
    """Create the cCRE track specification.

    Parameters
    ----------
    region (Dict[str, Any]): The genomic region for the track.

    Returns
    -------
    Dict[str, Any]
        The cCRE track specification.
    """
    return {
        "view": {"stroke": "lightgray"},
        "height": DEFAULT_TRACK_HEIGHT,
        "name": "ENCODE cCRE",
        "data": {
            "lazy": {
                "type": "bigbed",
                "url": "https://data.genomespy.app/sample-data/encodeCcreCombined.hg38.bb"
            }
        },
        "encoding": {
            "x": {
                "chrom": "chrom",
                "pos": "chromStart",
                "type": "locus",
                "scale": {
                    "domain": [
                        {"chrom": region["chrom"], "pos": region["start"]},
                        {"chrom": region["chrom"], "pos": region["end"]}
                    ]
                }
            },
            "x2": {
                "chrom": "chrom",
                "pos": "chromEnd"
            },
            "color": {
                "field": "ucscLabel",
                "type": "nominal",
                "scale": {
                    "domain": ["prom", "enhP", "enhD", "K4m3", "CTCF"],
                    "range": ["#FF0000", "#FFA700", "#FFCD00", "#FFAAAA", "#00B0F0"]
                }
            }
        },
        "mark": "rect"
    }

def create_gencode_track(region: Dict[str, Any]) -> Dict[str, Any]:
    """Create the Gencode track specification.

    Parameters
    ----------
    region (Dict[str, Any]): The genomic region for the track.

    Returns
    -------
    Dict[str, Any]
        The Gencode track specification.
    """
    return {
        "height": {"step": 28},  # Height per row
        "name": "Gencode v43",
        "viewportHeight": 280,  # Total viewport height
        "data": {
            "lazy": {
                "type": "gff3",
                "url": "https://data.genomespy.app/sample-data/gencode.v43.annotation.sorted.gff3.gz",
                "windowSize": 2000000,
                "debounceDomainChange": 300
            }
        },
        "transform": [
            {"type": "flatten"},
            {
                "type": "formula",
                "expr": "datum.attributes.gene_name[0]",
                "as": "gene_name"
            },
            {
                "type": "flatten",
                "fields": ["child_features"]
            },
            {
                "type": "flatten",
                "fields": ["child_features"],
                "as": ["child_feature"]
            },
            {
                "type": "project",
                "fields": [
                    "gene_name",
                    "child_feature.type",
                    "child_feature.strand",
                    "child_feature.seq_id",
                    "child_feature.start",
                    "child_feature.end",
                    "child_feature.attributes.gene_type",
                    "child_feature.attributes.transcript_type",
                    "child_feature.attributes.gene_id",
                    "child_feature.attributes.transcript_id",
                    "child_feature.attributes.transcript_name",
                    "child_feature.attributes.tag",
                    "source",
                    "child_feature.child_features"
                ],
                "as": [
                    "gene_name",
                    "type",
                    "strand",
                    "seq_id",
                    "start",
                    "end",
                    "gene_type",
                    "transcript_type",
                    "gene_id",
                    "transcript_id",
                    "transcript_name",
                    "tag",
                    "source",
                    "_child_features"
                ]
            },
            {
                "type": "collect",
                "sort": {
                    "field": ["seq_id", "start", "transcript_id"]
                }
            },
            {
                "type": "pileup",  # Add pileup transform to stack genes
                "start": "start",
                "end": "end",
                "as": "_lane"
            }
        ],
        "encoding": create_gencode_encoding(region),
        "layer": create_gencode_layers()
    }

def create_gencode_encoding(region: Dict[str, Any]) -> Dict[str, Any]:
    """Create the encoding specification for the Gencode track.

    Args:
        region (Dict[str, Any]): The genomic region for the track.

    Returns:
        Dict[str, Any]: The encoding specification.
    """
    return {
        "x": {
            "chrom": "seq_id",
            "pos": "start",
            "offset": 1,
            "type": "locus",
            "axis": {
                "orient": "top",
                "chromGrid": True,
                "chromGridColor": "lightgray",
                "grid": True,
                "chromGridDash": [3, 3],
                "gridDash": [1, 5],
                "gridColor": "#e0e0e0"
            },
            "scale": {
                "domain": [
                    {"chrom": region["chrom"], "pos": region["start"]},
                    {"chrom": region["chrom"], "pos": region["end"]}
                ]
            }
        },
        "x2": {
            "chrom": "seq_id",
            "pos": "end"
        },
        "y": {
            "field": "_lane",
            "type": "index",
            "axis": None,
            "scale": {
                "zoom": False,
                "reverse": True,
                "domain": [0, 50],
                "padding": 0.5
            }
        }
    }

def create_gencode_layers() -> list:
    """Create the layer specifications for the Gencode track.

    Returns
    -------
    list
        The list of layer specifications.
    """
    return [
        {
            "name": "gencode-transcript",
            "layer": [
                {
                    "name": "gencode-tooltip-trap",
                    "title": "GENCODE transcript",
                    "mark": {
                        "type": "rule",
                        "color": "#b0b0b0",
                        "opacity": 0,
                        "size": 7
                    }
                },
                {
                    "name": "gencode-transcript-body",
                    "mark": {
                        "type": "rule",
                        "color": "#b0b0b0",
                        "tooltip": None
                    }
                }
            ]
        },
        create_gencode_exons_layer(),
        create_gencode_labels_layer()
    ]

def create_gencode_exons_layer() -> Dict[str, Any]:
    """Create the exons layer specification for the Gencode track.

    Returns
    -------
    Dict[str, Any]
        The exons layer specification.
    """
    return {
        "name": "gencode-exons",
        "transform": [
            {
                "type": "flatten",
                "fields": ["_child_features"]
            },
            {
                "type": "flatten",
                "fields": ["_child_features"],
                "as": ["child_feature"]
            },
            {
                "type": "project",
                "fields": [
                    "gene_name",
                    "_lane",
                    "child_feature.type",
                    "child_feature.seq_id",
                    "child_feature.start",
                    "child_feature.end",
                    "child_feature.attributes.exon_number",
                    "child_feature.attributes.exon_id"
                ],
                "as": [
                    "gene_name",
                    "_lane",
                    "type",
                    "seq_id",
                    "start",
                    "end",
                    "exon_number",
                    "exon_id"
                ]
            }
        ],
        "layer": [
            create_exon_layer(),
            create_feature_layer(),
            create_utr_label_layer()
        ]
    }

def create_exon_layer() -> Dict[str, Any]:
    """Create the exon sublayer specification.

    Returns
    -------
    Dict[str, Any]
        The exon sublayer specification.
    """
    return {
        "title": "GENCODE exon",
        "transform": [
            {"type": "filter", "expr": "datum.type == 'exon'"}
        ],
        "mark": {
            "type": "rect",
            "minWidth": 0.5,
            "minOpacity": 0.5,
            "stroke": "#505050",
            "fill": "#fafafa",
            "strokeWidth": 1.0
        }
    }

def create_feature_layer() -> Dict[str, Any]:
    """Create the feature sublayer specification.

    Returns
    -------
    Dict[str, Any]
        The feature sublayer specification.
    """
    return {
        "title": "GENCODE feature",
        "transform": [
            {
                "type": "filter",
                "expr": "datum.type != 'exon' && datum.type != 'start_codon' && datum.type != 'stop_codon'"
            }
        ],
        "mark": {
            "type": "rect",
            "minWidth": 0.5,
            "minOpacity": 0,
            "strokeWidth": 1.0,
            "strokeOpacity": 0.0,
            "stroke": "gray"
        },
        "encoding": {
            "fill": {
                "field": "type",
                "type": "nominal",
                "scale": {
                    "domain": ["five_prime_UTR", "CDS", "three_prime_UTR"],
                    "range": ["#83bcb6", "#ffbf79", "#d6a5c9"]
                }
            }
        }
    }

def create_utr_label_layer() -> Dict[str, Any]:
    """Create the UTR label sublayer specification.

    Returns
    -------
    Dict[str, Any]
        The UTR label sublayer specification.
    """
    return {
        "transform": [
            {
                "type": "filter",
                "expr": "datum.type == 'three_prime_UTR' || datum.type == 'five_prime_UTR'"
            },
            {
                "type": "formula",
                "expr": "datum.type == 'three_prime_UTR' ? \"3'\" : \"5'\"",
                "as": "label"
            }
        ],
        "mark": {
            "type": "text",
            "color": "black",
            "size": 11,
            "opacity": 0.7,
            "paddingX": 2,
            "paddingY": 1.5,
            "tooltip": None
        },
        "encoding": {
            "text": {
                "field": "label"
            }
        }
    }

def create_gencode_labels_layer() -> Dict[str, Any]:
    """Create the labels layer specification for the Gencode track.

    Returns
    -------
    Dict[str, Any]
        The labels layer specification.
    """
    return {
        "name": "gencode-transcript-labels",
        "transform": [
            {
                "type": "formula",
                "expr": "(datum.strand == '-' ? '< ' : '') + datum.transcript_name + ' - ' + datum.transcript_id + (datum.strand == '+' ? ' >' : '')",
                "as": "label"
            }
        ],
        "mark": {
            "type": "text",
            "size": 10,
            "yOffset": 12,
            "tooltip": None,
            "color": "#505050"
        },
        "encoding": {
            "text": {
                "field": "label"
            }
        }
    }

