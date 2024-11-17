class DataNode:
    def __init__(self, name="root"):
        self._attributes = {}  # 存储动态属性
        self.name = name  # 为每个节点设置名称

    def __getattr__(self, key):
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

    def format_tree(self, level=0, html=False):
        # 递归生成树状结构
        indent = "  " * level
        html_indent = f"<div style='margin-left:{level * 20}px;'>"
        lines = []

        if not html:
            lines.append(f"{indent}- {self.name}")
        else:
            lines.append(f"{html_indent}- {self.name}</div><br>")

        for key, value in self._attributes.items():
            if isinstance(value, DataNode):
                lines.append(value.format_tree(level + 1, html))
            else:
                if not html:
                    lines.append(f"{indent}  {key}: {value}")
                else:
                    lines.append(
                        f"<div style='margin-left:{(level + 1) * 20}px;'>{key}: {value}</div><br>"
                    )
        return "\n".join(lines) if not html else "".join(lines)

    def __repr__(self):
        return self.format_tree()
