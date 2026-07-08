"""Update the tutorial dataset TOML registry from a local data repository."""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
DEFAULT_TUTORIAL_REPO = Path(r"..\easyclimate-tutorial")
DEFAULT_CONFIG = REPO_ROOT / "src" / "easyclimate" / "core" / "tutorial_data.toml"
DATA_SUFFIXES = {".nc", ".csv", ".CSV", ".grib"}


def parse_simple_toml(path: Path) -> dict[str, dict[str, object]]:
    if not path.exists():
        return {}

    config: dict[str, dict[str, object]] = {}
    section: str | None = None
    for raw_line in path.read_text(encoding="utf-8").splitlines():
        line = raw_line.strip()
        if not line or line.startswith("#"):
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1].strip()
            config.setdefault(section, {})
            continue
        if section is None or "=" not in line:
            continue
        key, value = line.split("=", 1)
        key = key.strip().strip('"')
        value = value.strip()
        if value.startswith('"') and value.endswith('"'):
            parsed_value: object = value[1:-1]
        else:
            parsed_value = int(value)
        config[section][key] = parsed_value
    return config


def quote(value: str) -> str:
    return '"' + value.replace("\\", "\\\\").replace('"', '\\"') + '"'


def format_value(value: object) -> str:
    if isinstance(value, int):
        return str(value)
    return quote(str(value))


def render_section(name: str, values: dict[str, object]) -> list[str]:
    lines = [f"[{name}]"]
    for key in sorted(values):
        lines.append(f"{quote(key)} = {format_value(values[key])}")
    lines.append("")
    return lines


def parse_md5_file(path: Path) -> str | None:
    if not path.exists():
        return None
    match = re.search(r"\b[0-9a-fA-F]{32}\b", path.read_text(encoding="utf-8"))
    if match is None:
        return None
    return match.group(0).lower()


def calculate_md5(path: Path) -> str:
    digest = hashlib.md5()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def infer_format(path: Path) -> int | str:
    if path.suffix in {".csv", ".CSV"}:
        return "csv"
    if path.suffix == ".grib":
        return "grib"
    return 4


def update_config(tutorial_repo: Path, config_path: Path) -> None:
    if not tutorial_repo.exists():
        raise FileNotFoundError(f"Tutorial repository not found: {tutorial_repo}")

    config = parse_simple_toml(config_path)
    file_formats = dict(config.get("file_formats", {}))
    tutorial_hashes = dict(config.get("tutorial_hashes", {}))
    external_urls = dict(config.get("external_urls", {}))
    external_hashes = dict(config.get("external_hashes", {}))

    for data_file in sorted(tutorial_repo.iterdir()):
        if not data_file.is_file() or data_file.suffix not in DATA_SUFFIXES:
            continue

        dataset_name = data_file.stem
        file_formats[dataset_name] = infer_format(data_file)

        md5 = parse_md5_file(data_file.with_suffix(".md5"))
        if md5 is None:
            md5 = calculate_md5(data_file)
        tutorial_hashes[data_file.name] = f"md5:{md5}"

    lines = [
        "# Tutorial dataset registry.",
        "#",
        "# Update this file with:",
        "#   python scripts/update_tutorial_config.py",
        "",
    ]
    lines.extend(render_section("file_formats", file_formats))
    lines.extend(render_section("tutorial_hashes", tutorial_hashes))
    lines.extend(render_section("external_urls", external_urls))
    lines.extend(render_section("external_hashes", external_hashes))

    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text("\n".join(lines), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--tutorial-repo",
        type=Path,
        default=DEFAULT_TUTORIAL_REPO,
        help="Path to the local easyclimate-tutorial repository.",
    )
    parser.add_argument(
        "--config",
        type=Path,
        default=DEFAULT_CONFIG,
        help="Path to the TOML file to update.",
    )
    args = parser.parse_args()

    update_config(args.tutorial_repo, args.config)
    print(f"Updated {args.config} from {args.tutorial_repo}")


if __name__ == "__main__":
    main()
