from IPython.display import HTML
import html as html_escape
from xarray.core.formatting_html import _load_static_files
import xarray as xr
import inspect
import json
from pathlib import Path
from typing import Union
import zarr


class DataNode:
    def __init__(self, name="root"):
        self._attributes = {}  # 存储动态属性
        self.name = name  # 为每个节点设置名称

    def __getattr__(self, key):
        # 过滤掉所有 IPython 的特殊属性请求
        if (
            key.startswith("_repr_")
            or key.startswith("_ipython_")
            or key.startswith("to_pandas")
            or key.startswith("toPandas")
            or key.startswith("__dataframe__")
        ):
            raise AttributeError(key)
        # 如果访问的属性不存在，自动创建嵌套的 DynamicClass
        if key not in self._attributes:
            self._attributes[key] = DataNode(name=key)
        return self._attributes[key]

    def __setattr__(self, key, value):
        if key in ["_attributes", "name"]:  # 避免覆盖内部存储或名称
            super().__setattr__(key, value)
        else:
            self._attributes[key] = value

    def __getitem__(self, key):
        # 支持路径风格访问
        if "/" in key:
            parts = key.split("/")
            node = self
            for part in parts:
                if part not in node._attributes:
                    raise KeyError(f"'{key}' not found")
                node = node._attributes[part]
            return node
        return self._attributes[key]

    def __setitem__(self, key, value):
        # 支持路径风格设置
        if "/" in key:
            parts = key.split("/")
            node = self
            for part in parts[:-1]:  # 遍历中间节点，自动创建
                if part not in node._attributes:
                    node._attributes[part] = DataNode(name=part)
                node = node._attributes[part]
            node._attributes[parts[-1]] = value
        else:
            self._attributes[key] = value

    def _repr_html_(self):
        """生成类似 xarray 的 HTML 表示"""
        # 直接返回字符串而不是 HTML 对象
        return self._format_html()

    def _format_html(self):
        """生成完整的 HTML 内容"""
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

        # 添加头部信息
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
        """递归生成节点的 HTML 表示"""
        node_id = (
            f"node-{id(node)}"
            if parent_id is None
            else f"{parent_id}-{html_escape.escape(node.name)}"
        )
        is_root = level == 0

        html = []

        # 节点标题行
        html.append(
            f"<div class='node-header{' root' if is_root else ''}' id='{node_id}-header'>"
        )

        # 添加折叠/展开按钮（如果有子节点）
        has_children = any(
            not k.startswith("_ipython_") for k in node._attributes.keys()
        )
        if has_children:
            html.append(
                f"<span class='toggle' onclick='easyclimateToggleNode(\"{node_id}\")'>▶</span>"
            )
        else:
            html.append("<span class='toggle-placeholder'></span>")

        # 节点名称
        html.append(f"<span class='node-name'>{html_escape.escape(node.name)}</span>")
        html.append("</div>")  # 关闭 node-header

        # 子节点容器
        html.append(
            f"<div class='node-children{' collapsed' if not is_root else ''}' id='{node_id}-children'>"
        )

        # 遍历属性，过滤掉 IPython 的特殊属性
        for key, value in sorted(node._attributes.items()):
            if key.startswith("_ipython_"):
                continue

            attr_id = f"{node_id}-{html_escape.escape(key)}"

            if isinstance(value, DataNode):
                html.append(self._format_node_html(value, level + 1, node_id))
            else:
                # 所有属性都有折叠功能
                html.append(f"<div class='node-attribute'>")
                html.append(
                    f"<div class='attr-header' onclick='easyclimateToggleAttr(\"{attr_id}\")'>"
                )
                html.append(f"<span id='{attr_id}-toggle' class='toggle'>▶</span>")
                html.append(
                    f"<span class='attr-name'>{html_escape.escape(key)}:</span>"
                )

                # 显示完整的类型信息
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
                html.append("</div>")  # 关闭 attr-header

                # 属性内容区域
                html.append(f"<div id='{attr_id}-content' style='display:none;'>")
                if isinstance(value, (xr.DataArray, xr.Dataset, xr.DataTree)):
                    # 使用xarray原生的HTML表示
                    html_repr = value._repr_html_()
                    html.append(f"<div class='xarray-html-repr'>{html_repr}</div>")
                else:
                    html.append(
                        f"<div class='attr-value'>{html_escape.escape(str(self._format_value(value)))}</div>"
                    )
                html.append("</div>")  # 关闭 content
                html.append("</div>")  # 关闭 node-attribute

        html.append("</div>")  # 关闭 node-children

        return "".join(html)

    def _format_value(self, value):
        """格式化值以便显示"""
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

    def format_tree(self, level=0, html=False):
        """保留原有的树形结构输出方法"""
        if html:
            return self._format_html()
        # 原有的纯文本格式输出
        indent = "  " * level
        lines = [f"{indent}- {self.name}"]
        for key, value in sorted(self._attributes.items()):
            if key.startswith("_ipython_"):
                continue
            if isinstance(value, DataNode):
                lines.append(value.format_tree(level + 1, html))
            elif isinstance(value, xr.DataArray):
                lines.append(
                    f"{indent}  {key}: <xarray.DataArray> (shape: {value.shape}, dtype: {value.dtype})"
                )
            elif inspect.ismethod(value) or inspect.isfunction(value):
                # 跳过方法和函数
                continue
            else:
                lines.append(f"{indent}  {key}: {value}")
        return "\n".join(lines)

    def __repr__(self):
        return self.format_tree()

    def _repr_pretty_(self, p, cycle):
        """支持 IPython 的 pretty printing"""
        if cycle:
            p.text(f"{self.name}(...)")
        else:
            p.text(self.format_tree())

    # 明确禁用 _repr_mimebundle_
    def _repr_mimebundle_(self, include=None, exclude=None):
        return None

    # 或者更彻底地，禁用所有可能的表示方法
    def __dir__(self):
        # 过滤掉所有 _repr_* 方法
        return [
            attr
            for attr in super().__dir__()
            if not attr.startswith("_repr_")
            and not attr.startswith("_ipython_")
            and not attr.startswith("to_pandas")
        ]

    def to_zarr(self, filepath: Union[str, Path]):
        filepath = Path(filepath)
        filepath.mkdir(parents=True, exist_ok=True)

        metadata = {"name": self.name, "attributes": {}}

        for key, value in self._attributes.items():
            if isinstance(value, DataNode):
                child_path = filepath / key
                value.to_zarr(child_path)
                metadata["attributes"][key] = {
                    "__type__": "DataNode",
                    "path": str(child_path),
                }
            elif isinstance(value, (xr.DataArray, xr.Dataset, xr.DataTree)):
                zarr_path = filepath / f"{key}.zarr"
                # Disable consolidated metadata
                value.to_zarr(zarr_path, consolidated=False, mode="w")
                metadata["attributes"][key] = {
                    "__type__": "xarray",
                    "path": str(zarr_path),
                }
            else:
                metadata["attributes"][key] = value

        with open(filepath / "metadata.json", "w", encoding="utf-8") as f:
            json.dump(metadata, f, indent=2, ensure_ascii=False)

    @classmethod
    def load(cls, filepath: Union[str, Path]):
        """
        Load DataNode from file

        Parameters
        -------------------------
            filepath: File path

        Returns
        -------------------------
            Loaded DataNode object
        """
        filepath = Path(filepath)

        # Load metadata
        with open(filepath / "metadata.json", "r", encoding="utf-8") as f:
            metadata = json.load(f)

        node = cls(name=metadata["name"])

        for key, value in metadata["attributes"].items():
            if isinstance(value, dict) and value.get("__type__") == "DataNode":
                # Recursively load child node
                child_path = Path(value["path"])
                node._attributes[key] = cls.load(child_path)
            elif isinstance(value, dict) and value.get("__type__") == "xarray":
                # Load xarray data
                zarr_path = Path(value["path"])
                node._attributes[key] = xr.open_zarr(zarr_path, consolidated=False)
            else:
                # Load normal data
                node._attributes[key] = value

        return node
