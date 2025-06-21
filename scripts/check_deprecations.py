#!/usr/bin/env python3
"""
Deprecation Checker with Rich Color Output

Usage:
    python check_deprecations.py [--next-version X.Y.Z] [--json FILE] [--html FILE]

Example:
    # Basic inspection
    python scripts/check_deprecations.py

    # Specify the next version and save the JSON report
    python scripts/check_deprecations.py --next-version 2.0.0 --json deprecations.json

Sample:
    @deprecated(
        version="2025.4.0",
        removal_version="2025.7.0",
        replacement="easyclimate.field.typhoon.calc_potential_intensity_Bister_Emanuel_2002",
    )
    def calc_potential_intensity_Bister_Emanuel_2002(
        sst_data: xr.DataArray,
        sst_data_units: Literal["celsius", "kelvin", "fahrenheit"],
        surface_pressure_data: xr.DataArray,
        surface_pressure_data_units: Literal["hPa", "Pa", "mbar"],
        temperature_data: xr.DataArray,
        temperature_data_units: Literal["celsius", "kelvin", "fahrenheit"],
        specific_humidity_data: xr.DataArray,
        specific_humidity_data_units: str,
        vertical_dim: str,
        vertical_dim_units: str,
        CKCD: float = 0.9,
        ascent_flag: bool = False,
        diss_flag: bool = True,
        V_reduc: float = 0.8,
        ptop: float = 50,
        miss_handle: bool = True,
    ) -> xr.Dataset:
"""

import ast
import importlib.util
import json
import sys
import io
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Optional

from packaging import version
from rich.console import Console
from rich.panel import Panel
from rich.syntax import Syntax
from rich.table import Table
from rich.text import Text
from rich.style import Style

# Windows encoding compatibility processing
if sys.platform == "win32":
    sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace")


class DeprecationChecker:
    def __init__(self, package_name: str = "easyclimate"):
        # Replace Unicode characters with Windows-compatible alternative text
        if sys.platform == "win32":
            self.scan_icon = "[SEARCH]"
        else:
            self.scan_icon = "🔍"

        self.package_name = package_name
        self.current_version = self._get_package_version()
        self.next_version = None
        self.deprecations = defaultdict(list)
        self.console = Console()

    def _get_package_version(self) -> str:
        """Get package version from version.py"""
        try:
            spec = importlib.util.find_spec(f"{self.package_name}.version")
            if spec and spec.loader:
                version_module = importlib.util.module_from_spec(spec)
                spec.loader.exec_module(version_module)
                return version_module.__version__
        except (ImportError, AttributeError):
            pass

        version_file = Path(f"src/{self.package_name}/version.py")
        if version_file.exists():
            with open(version_file, "r", encoding="utf-8") as f:
                content = f.read()
            tree = ast.parse(content)
            for node in ast.walk(tree):
                if isinstance(node, ast.Assign) and len(node.targets) == 1:
                    if getattr(node.targets[0], "id", None) == "__version__":
                        # if isinstance(node.value, ast.Str):
                        # return node.value.s
                        if isinstance(node.value, ast.Constant):
                            return node.value.value
                        elif isinstance(node.value, ast.Constant):
                            return str(node.value.value)
        return "0.0.0"

    def set_next_version(self, next_version: str):
        """Set the next planned version"""
        self.next_version = next_version

    def check_file(self, filepath: Path):
        """Check a single Python file for deprecations"""
        with open(filepath, "r", encoding="utf-8") as f:
            content = f.read()

        try:
            tree = ast.parse(content)
        except SyntaxError as e:
            self.console.print(f"[yellow]⚠️  Syntax error in {filepath}: {e}[/]")
            return

        for node in ast.walk(tree):
            if not isinstance(node, ast.FunctionDef):
                continue

            for decorator in node.decorator_list:
                if not isinstance(decorator, ast.Call):
                    continue

                if not isinstance(decorator.func, ast.Name):
                    continue

                if decorator.func.id == "deprecated":
                    self._process_deprecated_node(node, decorator, filepath)

    def _process_deprecated_node(
        self, node: ast.FunctionDef, decorator: ast.Call, filepath: Path
    ):
        """Process a single @deprecated decorator"""
        dep_version = None
        removal_version = None
        replacement = None

        for keyword in decorator.keywords:
            if keyword.arg == "version":
                dep_version = self._get_ast_value(keyword.value)
            elif keyword.arg == "removal_version":
                removal_version = self._get_ast_value(keyword.value)
            elif keyword.arg == "replacement":
                replacement = self._get_ast_value(keyword.value)

        if not dep_version or not removal_version:
            return

        status = self._determine_status(dep_version, removal_version)
        code_line = self._get_code_line(filepath, node.lineno)

        self.deprecations[status].append(
            {
                "name": node.name,
                "file": str(filepath),
                "line": node.lineno,
                "code": code_line,
                "version": dep_version,
                "removal_version": removal_version,
                "replacement": replacement,
                "status": status,
            }
        )

    def _get_ast_value(self, node: ast.AST) -> Optional[str]:
        """Extract value from AST node"""
        # if isinstance(node, ast.Str):
        #     return node.s
        if isinstance(node, ast.Constant):
            return node.value
        elif isinstance(node, ast.Constant):
            return str(node.value)
        elif isinstance(node, ast.Name):
            return node.id
        return None

    def _get_code_line(self, filepath: Path, lineno: int) -> str:
        """Get specific line of code from file"""
        with open(filepath, "r", encoding="utf-8") as f:
            for i, line in enumerate(f, 1):
                if i == lineno:
                    return line.strip()
        return ""

    def _determine_status(self, dep_version: str, removal_version: str) -> str:
        """Determine deprecation status"""
        current_v = version.parse(self.current_version)
        dep_v = version.parse(dep_version)
        removal_v = version.parse(removal_version)

        if self.next_version and version.parse(self.next_version) >= removal_v:
            return "to_be_removed"
        elif current_v >= removal_v:
            return "removed_but_still_used"
        elif current_v >= dep_v:
            return "deprecated"
        else:
            return "future_deprecation"

    def check_project(self, project_dir: str = "src"):
        """Check entire project directory"""
        project_path = Path(project_dir)
        self.console.print(
            f"[bold]{self.scan_icon} Scanning {project_path} for deprecations...[/]"
        )

        with self.console.status("[bold green]Checking files...", spinner="dots"):
            for py_file in project_path.rglob("*.py"):
                if py_file.name == "version.py":
                    continue
                self.check_file(py_file)

    def print_report(self):
        """Print beautiful rich formatted report"""
        # Header Panel
        header = Text.assemble(
            ("Deprecation Report for ", Style(bold=True)),
            (f"{self.package_name} ", Style(bold=True, color="blue")),
            (f"v{self.current_version}", Style(bold=True, color="cyan")),
        )
        if self.next_version:
            header.append(f"\nNext version: v{self.next_version}", style="bold yellow")

        self.console.print(
            Panel.fit(
                header,
                border_style="blue",
                padding=(1, 2),
                title="📝 Deprecation Status",
                subtitle=f"Generated {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}",
            )
        )

        # Status Summary
        summary_table = Table.grid(padding=(0, 2))
        summary_table.add_column(style="bold")
        summary_table.add_column(style="bold")

        status_counts = {
            status: len(items) for status, items in self.deprecations.items()
        }
        total = sum(status_counts.values())

        if total == 0:
            self.console.print("[bold green]✅ No deprecations found![/]")
            return

        summary_table.add_row("Total deprecations found:", f"[bold]{total}[/]")
        self._add_status_row(
            summary_table, "removed_but_still_used", "🚨 Should Be Removed", "red"
        )
        self._add_status_row(
            summary_table,
            "to_be_removed",
            "⚠️ To Be Removed Next Version",
            "dark_orange",
        )
        self._add_status_row(
            summary_table, "deprecated", "🔶 Currently Deprecated", "yellow"
        )
        self._add_status_row(
            summary_table, "future_deprecation", "🔹 Future Deprecations", "blue"
        )

        self.console.print(
            Panel.fit(summary_table, title="📊 Summary", border_style="green")
        )

        # Detailed Tables
        status_order = [
            "removed_but_still_used",
            "to_be_removed",
            "deprecated",
            "future_deprecation",
        ]

        for status in status_order:
            items = self.deprecations.get(status, [])
            if not items:
                continue

            color = {
                "removed_but_still_used": "red",
                "to_be_removed": "dark_orange",
                "deprecated": "yellow",
                "future_deprecation": "blue",
            }.get(status, "white")

            table = Table(
                title=f"[bold {color}]{self._get_status_title(status)} ({len(items)})[/]",
                border_style=color,
                show_header=True,
                header_style=f"bold {color}",
                expand=True,
            )
            table.add_column("Function/Method", style="bold")
            table.add_column("Location")
            table.add_column("Deprecated Since")
            table.add_column("Removal Version")
            table.add_column("Replacement")

            for item in items:
                table.add_row(
                    item["name"],
                    f"{Path(item['file']).name}:{item['line']}",
                    f"v{item['version']}",
                    f"v{item['removal_version']}",
                    item["replacement"] or "[dim]-[/]",
                )

            self.console.print(table)

            # Show code example
            sample = items[0]
            self.console.print(f"[bold]Example from {Path(sample['file']).name}:[/]")
            with open(sample["file"], "r", encoding="utf-8") as f:
                lines = f.readlines()
                start_line = max(0, sample["line"] - 2)
                end_line = min(len(lines), sample["line"] + 2)
                code = "".join(lines[start_line:end_line])
                self.console.print(
                    Syntax(
                        code,
                        "python",
                        line_numbers=True,
                        line_range=(start_line + 1, end_line + 1),
                        highlight_lines={sample["line"]},
                        word_wrap=True,
                    )
                )
            self.console.print()

    def _add_status_row(self, table, status_key, label, color):
        """Add a status row to summary table"""
        count = len(self.deprecations.get(status_key, []))
        if count > 0:
            table.add_row(f"[{color}]{label}:[/]", f"[bold {color}]{count}[/]")

    def _get_status_title(self, status):
        """Get human-readable status title"""
        return {
            "removed_but_still_used": "Should Be Removed (But Still in Code)",
            "to_be_removed": "To Be Removed in Next Version",
            "deprecated": "Currently Deprecated",
            "future_deprecation": "Future Deprecations",
        }.get(status, status.replace("_", " ").title())


def main():
    import argparse

    parser = argparse.ArgumentParser(
        description="Check project deprecations with rich output"
    )
    parser.add_argument("--next-version", help="Next planned version number")
    parser.add_argument("--json", help="Save report to JSON file")
    args = parser.parse_args()

    checker = DeprecationChecker()
    if args.next_version:
        checker.set_next_version(args.next_version)

    checker.check_project()
    checker.print_report()

    if args.json:
        report_data = {
            "timestamp": datetime.now().isoformat(),
            "package": checker.package_name,
            "current_version": checker.current_version,
            "next_version": checker.next_version,
            "deprecations": [
                {**item, "status": status}
                for status, items in checker.deprecations.items()
                for item in items
            ],
        }
        with open(args.json, "w", encoding="utf-8") as f:
            json.dump(report_data, f, indent=2)
        checker.console.print(f"\n[green]✓ JSON report saved to {args.json}[/]")

    if any(
        status in checker.deprecations
        for status in ["removed_but_still_used", "to_be_removed"]
    ):
        sys.exit(1)


if __name__ == "__main__":
    main()
