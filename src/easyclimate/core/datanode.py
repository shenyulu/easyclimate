"""
Hierarchical data structure that dynamically manages attributes and nested nodes
"""

import html as html_escape
from xarray.core.formatting_html import _load_static_files
from typing import Literal
import xarray as xr
import inspect
import json
from pathlib import Path
from typing import Union
import zarr

__all__ = ["DataNode", "open_datanode"]


class DataNode:
    """
    A hierarchical data structure that dynamically manages attributes and nested nodes.

    The :py:class:`DataNode <DataNode>` class provides a flexible way to organize and access data in a tree-like
    structure. It supports automatic creation of nested nodes, path-style access,
    and rich HTML representation in Jupyter environments.

    Parameters
    ----------
    name : :py:class:`str <str>`, optional
        The name of the root node (default: ``"root"``).
    """

    def __init__(self, name="root"):
        """
        Initialize a new DataNode instance.

        Parameters
        ----------
        name : :py:class:`str <str>`, optional
            The name of the node (default: "root").
        """
        self._attributes = {}  # Storing Dynamic Properties
        self.name = name  # Setting names for each node

    def __getattr__(self, key):
        """
        Dynamically access or create node attributes.

        Automatically creates nested :py:class:`DataNode <DataNode>` for non-existent attributes,
        while filtering out special IPython/Jupyter attribute requests.

        Parameters
        ----------
        key : :py:class:`str <str>`
            The attribute name to access.

        Returns
        -------
        Any
            The requested attribute or a new :py:class:`DataNode <DataNode>` if attribute doesn't exist.

        Raises
        ------
        AttributeError
            If the attribute is a special IPython/Jupyter attribute.
        """
        # Filter out all IPython requests for special attributes
        if (
            key.startswith("_repr_")
            or key.startswith("_ipython_")
            or key.startswith("to_pandas")
            or key.startswith("toPandas")
            or key.startswith("__dataframe__")
        ):
            raise AttributeError(key)
        # Automatically creates a nested DynamicClass if the accessed property does not exist.
        if key not in self._attributes:
            self._attributes[key] = DataNode(name=key)
        return self._attributes[key]

    def __setattr__(self, key, value):
        """
        Set node attributes while protecting internal attributes.

        Parameters
        ----------
        key : :py:class:`str <str>`
            The attribute name to set.
        value : Any
            The value to assign to the attribute.
        """
        if key in [
            "_attributes",
            "name",
        ]:  # Avoid overwriting internal storage or names
            super().__setattr__(key, value)
        else:
            self._attributes[key] = value

    def __getitem__(self, key):
        """
        Access attributes using path-style notation (e.g., "path/to/attribute").

        Parameters
        ----------
        key : :py:class:`str <str>`
            The attribute path to access, with '/' separators for nested nodes.

        Returns
        -------
        Any
            The value at the specified path.

        Raises
        ------
        KeyError
            If any part of the path doesn't exist.
        """
        # Support for path style access
        if "/" in key:
            parts = key.split("/")
            node = self
            for part in parts:
                if part not in node._attributes:
                    raise KeyError(f"'{key}' not found")
                node = node._attributes[part]
            return node
        return self._attributes[key]

    def __contains__(self, key):
        return key in self._attributes

    def __setitem__(self, key, value):
        """
        Set attributes using path-style notation, creating intermediate nodes as needed.

        Parameters
        ----------
        key : :py:class:`str <str>`
            The attribute path to set, with '/' separators for nested nodes.
        value : Any
            The value to assign at the specified path.
        """
        # Support for path style settings
        if "/" in key:
            parts = key.split("/")
            node = self
            for part in parts[
                :-1
            ]:  # Iterate over intermediate nodes and automatically create
                if part not in node._attributes:
                    node._attributes[part] = DataNode(name=part)
                node = node._attributes[part]
            node._attributes[parts[-1]] = value
        else:
            self._attributes[key] = value

    def _repr_html_(self):
        """
        Generate HTML representation for Jupyter notebooks.

        Returns
        -------
        :py:class:`str <str>`
            HTML string representing the node and its contents.
        """
        # Returns a string directly instead of an HTML object
        return self._format_html()

    def _format_html(self):
        """
        Generate complete HTML representation including styles and scripts.

        Returns
        -------
        :py:class:`str <str>`
            Complete HTML document as a string.
        """
        html = []
        icons_svg, css_style = _load_static_files()

        # Add css
        html.append(
            f"""
        <style>
            .datanode-container {{
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
                font-size: 14px;
                line-height: 1.4;
                margin: 10px;
            }}
            .xarray-html-repr {{
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .node-header {{
                display: flex;
                align-items: center;
                cursor: pointer;
                padding: 2px 0;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .node-header.root {{
                font-weight: bold;
                font-size: 1.1em;
                margin-bottom: 5px;
            }}
            .toggle {{
                margin-right: 5px;
                color: #666;
                font-size: 10px;
                width: 10px;
                display: inline-block;
            }}
            .toggle-placeholder {{
                margin-right: 15px;
                visibility: hidden;
            }}
            .node-name {{
                font-weight: bold;
                color: #0366d6;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .node-children {{
                margin-left: 20px;
                border-left: 1px dotted #ddd;
                padding-left: 10px;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .node-children.collapsed {{
                display: none;
            }}
            .node-attribute {{
                display: flex;
                flex-direction: column;
                margin: 2px 0;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .attr-header {{
                display: flex;
                align-items: center;
                cursor: pointer;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .attr-name {{
                color: #d63384;
                margin-right: 5px;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .attr-value {{
                color: #333;
                margin-right: 10px;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .attr-type {{
                color: #6f42c1;
                font-style: italic;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
            .xarray-html-repr {{
                margin: 10px 0;
                margin-left: 15px;
                font-family: "Helvetica Neue", Helvetica, Arial, "PingFang SC", "Hiragino Sans GB", "Heiti SC", "Microsoft YaHei", "WenQuanYi Micro Hei", sans-serif;
            }}
        </style>
        """
        )

        # Add header information
        html.append(
            f"""
        <div class="datanode-container">
        {icons_svg}<style>{css_style}</style>


        <div class="xr-header">
            <div class="xr-obj-type">{html_escape.escape('easyclimate.'+type(self).__name__)}</div>
            <div class="xr-array-name">{html_escape.escape("'" + self.name + "'")}</div>
        </div>
        """
        )

        html.append(
            """
        <script>
        function easyclimateToggleNode(nodeId) {
            const header = document.getElementById(nodeId + '-header');
            const children = document.getElementById(nodeId + '-children');
            const toggle = header.querySelector('.toggle');

            if (children.classList.contains('collapsed')) {
                children.classList.remove('collapsed');
                toggle.textContent = '▼';
            } else {
                children.classList.add('collapsed');
                toggle.textContent = '▶';
            }
        }

        function easyclimateToggleAttr(attrId) {
            const content = document.getElementById(attrId + '-content');
            const toggle = document.getElementById(attrId + '-toggle');

            if (content.style.display === 'none') {
                content.style.display = 'block';
                toggle.textContent = '▼';
            } else {
                content.style.display = 'none';
                toggle.textContent = '▶';
            }
        }
        </script>
        """
        )

        html.append(self._format_node_html(self))
        html.append("</div>")
        return "".join(html)

    def _format_node_html(self, node, level=0, parent_id=None):
        """
        Recursively generate HTML for a node and its children.

        Parameters
        ----------
        node : :py:class:`DataNode <DataNode>`
            The node to format.
        level : :py:class:`int <int>`, optional
            Current nesting level (default: 0).
        parent_id : :py:class:`str <str>`, optional
            ID of parent node for DOM construction (default: None).

        Returns
        -------
        :py:class:`str <str>`
            HTML string representing the node.
        """
        node_id = (
            f"node-{id(node)}"
            if parent_id is None
            else f"{parent_id}-{html_escape.escape(node.name)}"
        )
        is_root = level == 0

        html = []

        # Node title line
        html.append(
            f"<div class='node-header{' root' if is_root else ''}' id='{node_id}-header'>"
        )

        # Add collapse/expand buttons (if there are child nodes)
        has_children = any(
            not k.startswith("_ipython_") for k in node._attributes.keys()
        )
        if has_children:
            html.append(
                f"<span class='toggle' onclick='easyclimateToggleNode(\"{node_id}\")'>▶</span>"
            )
        else:
            html.append("<span class='toggle-placeholder'></span>")

        # Node name
        html.append(f"<span class='node-name'>{html_escape.escape(node.name)}</span>")
        html.append("</div>")  # 关闭 node-header

        # Child node container
        html.append(
            f"<div class='node-children{' collapsed' if not is_root else ''}' id='{node_id}-children'>"
        )

        # Iterate over properties, filtering out IPython specific properties
        for key, value in sorted(node._attributes.items()):
            if key.startswith("_ipython_"):
                continue

            attr_id = f"{node_id}-{html_escape.escape(key)}"

            if isinstance(value, DataNode):
                html.append(self._format_node_html(value, level + 1, node_id))
            else:
                # All properties have collapsing functionality
                html.append(f"<div class='node-attribute'>")
                html.append(
                    f"<div class='attr-header' onclick='easyclimateToggleAttr(\"{attr_id}\")'>"
                )
                html.append(f"<span id='{attr_id}-toggle' class='toggle'>▶</span>")
                html.append(
                    f"<span class='attr-name'>{html_escape.escape(key)}:</span>"
                )

                # Show full type information
                if isinstance(value, (xr.DataArray, xr.Dataset, xr.DataTree)):
                    type_name = f"xarray.{type(value).__name__}"
                else:
                    type_name = type(value).__name__

                type_name = (
                    f"xarray.{type(value).__name__}"
                    if isinstance(value, xr.DataArray)
                    or isinstance(value, xr.Dataset)
                    or isinstance(value, xr.DataTree)
                    else type(value).__name__
                )
                html.append(f"<span class='attr-type'>{type_name}</span>")
                html.append("</div>")  # Close attr-header

                # Attribute Content Area
                html.append(f"<div id='{attr_id}-content' style='display:none;'>")
                if isinstance(value, (xr.DataArray, xr.Dataset, xr.DataTree)):
                    # Using xarray native HTML representation
                    html_repr = value._repr_html_()
                    html.append(f"<div class='xarray-html-repr'>{html_repr}</div>")
                else:
                    html.append(
                        f"<div class='attr-value'>{html_escape.escape(str(self._format_value(value)))}</div>"
                    )
                html.append("</div>")  # Close content
                html.append("</div>")  # Close node-attribute

        html.append("</div>")  # Close node-children

        return "".join(html)

    def _format_value(self, value):
        """
        Format values for display, truncating long sequences.

        Parameters
        ----------
        value : Any
            The value to format.

        Returns
        -------
        :py:class:`str <str>`
            Formatted string representation of the value.
        """
        if isinstance(value, (list, tuple)) and len(value) > 3:
            return f"[{value[0]}, {value[1]}, ..., {value[-1]}] (length: {len(value)})"
        elif isinstance(value, dict):
            return (
                "{"
                + ", ".join(f"{k}: ..." for k in list(value.keys())[:3])
                + "}"
                + (f" (length: {len(value)})" if len(value) > 3 else "")
            )
        return str(value)

    def format_tree(self, level=0, html=False, is_last_list=None):
        """
        Generate a tree-structured representation of the node.

        Parameters
        ----------
        level : :py:class:`int <int>`, optional
            Current indentation level (default: 0).
        html : :py:class:`bool <bool>`, optional
            Whether to generate HTML output (default: False).
        is_last_list : :py:class:`list <list>` of :py:class:`bool <bool>`, optional
            Track position in hierarchy for proper indentation (default: None).

        Returns
        -------
        :py:class:`str <str>`
            Formatted tree representation.
        """
        if html:
            return self._format_html()

        if is_last_list is None:
            is_last_list = []

        # Header line
        if level == 0:
            lines = [f"<easyclimate.DataNode '{self.name}'>"]
            lines.append(f"{self.name}: /")
        else:
            lines = []

        # Get sorted attributes
        sorted_attrs = sorted(self._attributes.items())

        for i, (key, value) in enumerate(sorted_attrs):
            if key.startswith("_ipython_"):
                continue

            current_is_last = i == len(sorted_attrs) - 1

            # Build the prefix for this line
            prefix = ""
            for is_last in is_last_list:
                prefix += "    " if is_last else "│   "

            connector = "└── " if current_is_last else "├── "

            if isinstance(value, DataNode):
                lines.append(f"{prefix}{connector}{key}: /")
                # For DataNodes, we need to handle both attributes and xarray objects
                child_lines = value.format_tree(
                    level + 1, html, is_last_list + [current_is_last]
                ).split("\n")
                # Skip the header line and add the rest
                if len(child_lines) > 1:
                    lines.extend(child_lines[0:])
            elif isinstance(value, (xr.DataArray, xr.Dataset)):
                dims = ", ".join(f"{k}: {v}" for k, v in value.sizes.items())
                lines.append(
                    f"{prefix}{connector}{key}: <xarray.{type(value).__name__}>"
                )
                lines.append(
                    f"{prefix}{'    ' if current_is_last else '│   '}Dimensions:  ({dims})"
                )

                # Handle coordinates
                if hasattr(value, "coords") and value.coords:
                    lines.append(
                        f"{prefix}{'    ' if current_is_last else '│   '}Coordinates:"
                    )
                    for coord in sorted(value.coords):
                        coord_dims = ", ".join(
                            f"{k}: {v}" for k, v in value.coords[coord].sizes.items()
                        )
                        lines.append(
                            f"{prefix}{'    ' if current_is_last else '│   '}  * {coord.ljust(10)} ({coord_dims}): {value.coords[coord].dtype}"
                        )

                # Handle data variables for Datasets
                if hasattr(value, "data_vars") and value.data_vars:
                    lines.append(
                        f"{prefix}{'    ' if current_is_last else '│   '}Data variables:"
                    )
                    for var in sorted(value.data_vars):
                        var_dims = ", ".join(
                            f"{k}: {v}" for k, v in value.data_vars[var].sizes.items()
                        )
                        lines.append(
                            f"{prefix}{'    ' if current_is_last else '│   '}    {var.ljust(10)} ({var_dims}): {value.data_vars[var].dtype}"
                        )
            else:
                lines.append(f"{prefix}{connector}{key}: {value}")

        return "\n".join(lines)

    def __repr__(self):
        return self.format_tree()

    def _repr_pretty_(self, p, cycle):
        """Support the IPython of pretty printing"""
        if cycle:
            p.text(f"{self.name}(...)")
        else:
            p.text(self.format_tree())

    # Explicitly prohibited _repr_mimebundle_
    def _repr_mimebundle_(self, include=None, exclude=None):
        return None

    # Or, more radically, disable all possible representations
    def __dir__(self):
        # Filter out all _repr_* methods
        return [
            attr
            for attr in super().__dir__()
            if not attr.startswith("_repr_")
            and not attr.startswith("_ipython_")
            and not attr.startswith("to_pandas")
        ]

    def to_zarr(self, filepath: Union[str, Path], zarr_format: Literal[2, 3] = 2):
        """
        Save the DataNode and its contents to a Zarr storage format.

        Parameters
        ----------
        filepath : Union[str, Path]
            Directory path to save the data.
        zarr_format : Literal[2, 3], optional
            Zarr storage format version (default: 2).
        """
        filepath = Path(filepath)
        filepath.mkdir(parents=True, exist_ok=True)

        metadata = {"name": self.name, "attributes": {}}

        for key, value in self._attributes.items():
            if isinstance(value, DataNode):
                child_path = filepath / key
                value.to_zarr(child_path)
                metadata["attributes"][key] = {
                    "__type__": "DataNode",
                    "path": str(child_path).replace(
                        "\\", "/"
                    ),  # Convert to forward slashes,
                }
            elif isinstance(value, (xr.DataArray, xr.Dataset, xr.DataTree)):
                zarr_path = filepath / f"{key}.zarr"
                # Disable consolidated metadata
                value.to_zarr(
                    zarr_path, consolidated=False, mode="w", zarr_format=zarr_format
                )
                metadata["attributes"][key] = {
                    "__type__": "xarray",
                    "path": str(zarr_path).replace(
                        "\\", "/"
                    ),  # Convert to forward slashes,
                }
            else:
                metadata["attributes"][key] = value

        with open(filepath / "metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

    @classmethod
    def load(cls, filepath: Union[str, Path]):
        """
        Load a DataNode from a Zarr storage directory.

        Parameters
        ----------
        filepath : Union[str, Path]
            Directory path containing the saved DataNode.

        Returns
        -------
        :py:class:`DataNode <DataNode>`
            The reconstructed DataNode with all its attributes.
        """
        filepath = Path(filepath)

        # Load metadata
        with open(filepath / "metadata.json", "r", encoding="utf-8") as f:
            metadata = json.load(f)

        node = cls(name=metadata["name"])

        for key, value in metadata["attributes"].items():
            if isinstance(value, dict) and value.get("__type__") == "DataNode":
                # Recursively load child node - convert path to Path object
                child_path = Path(value["path"].replace("\\", "/"))
                node._attributes[key] = cls.load(child_path)
            elif isinstance(value, dict) and value.get("__type__") == "xarray":
                # Load xarray data - convert path to Path object
                zarr_path = Path(value["path"].replace("\\", "/"))
                node._attributes[key] = xr.open_zarr(zarr_path, consolidated=False)
            else:
                # Load normal data
                node._attributes[key] = value

        return node


def open_datanode(filepath: str) -> DataNode:
    """
    Load a DataNode object from native location.

    This function provides a convenient way to load a DataNode that was previously saved
    using the ``DataNode.to_zarr()`` method.

    Parameters
    ----------
    filepath : :py:class:`str <str>`
        The path to the directory containing the saved DataNode data.
        This should be the same path used with DataNode.to_zarr().

    Returns
    -------
    :py:class:`DataNode <DataNode>`
        The loaded DataNode object with all its attributes and nested structure.

    Examples
    --------
    >>> node = open_datanode("path/to/saved_node")
    >>> node.some_attribute  # Access attributes as usual

    .. minigallery::
        :add-heading: Example(s) related to the function

        ./dynamic_docs/plot_multieof.py
    """
    return DataNode.load(filepath)
