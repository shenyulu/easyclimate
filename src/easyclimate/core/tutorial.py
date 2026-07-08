"""
Tutorials in the `easyclimate` documentation
"""

# idea borrowed from xarray

from __future__ import annotations

import os
import pooch
import requests
import pathlib
import re
import warnings
from typing import TYPE_CHECKING
from glob import glob
import xarray as xr
from pathlib import Path
import pandas as pd
from rich.progress import (
    Progress,
    BarColumn,
    DownloadColumn,
    TransferSpeedColumn,
    TimeRemainingColumn,
    TextColumn,
)

if TYPE_CHECKING:
    from xarray.backends.api import T_Engine

__all__ = ["open_tutorial_dataset"]

_default_cache_dir_name = "easylimate_tutorial_data"
base_url = "https://github.com/shenyulu/easyclimate-tutorial"
version = "main"
download_retries = 3
_known_hash_cache = {}
_tutorial_config_path = Path(__file__).with_name("tutorial_data.toml")


def _construct_cache_dir(path):
    if isinstance(path, os.PathLike):
        path = os.fspath(path)
    elif path is None:
        path = pooch.os_cache(_default_cache_dir_name)
    return path


def _parse_simple_toml(path: Path) -> dict[str, dict[str, object]]:
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


def _load_tutorial_config(path: Path) -> dict[str, dict[str, object]]:
    try:
        import tomllib
    except ModuleNotFoundError:
        return _parse_simple_toml(path)

    with path.open("rb") as stream:
        return tomllib.load(stream)


_tutorial_config = _load_tutorial_config(_tutorial_config_path)
file_formats = _tutorial_config.get("file_formats", {})
tutorial_hashes = _tutorial_config.get("tutorial_hashes", {})
external_urls = _tutorial_config.get("external_urls", {})
external_hashes = _tutorial_config.get("external_hashes", {})


def _check_netcdf_engine_installed(name):
    version = file_formats.get(name)
    if version == 3:
        try:
            import scipy  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either scipy or "
                    "netCDF4 to be installed."
                )
    if version == 4:
        try:
            import h5netcdf  # noqa
        except ImportError:
            try:
                import netCDF4  # noqa
            except ImportError:
                raise ImportError(
                    f"opening tutorial dataset {name} requires either h5netcdf "
                    "or netCDF4 to be installed."
                )
    if version == "csv":
        try:
            import pandas  # noqa
        except ImportError:
            raise ImportError(
                f"opening tutorial dataset {name} requires pandas " " to be installed."
            )


class RichDownloader:
    def __init__(self, headers=None):
        self.progress = Progress(
            TextColumn("[bold blue]{task.fields[filename]}", justify="right"),
            BarColumn(bar_width=None),
            "[progress.percentage]{task.percentage:>3.1f}%",
            "•",
            DownloadColumn(),
            "•",
            TransferSpeedColumn(),
            "•",
            TimeRemainingColumn(),
        )
        self.task_id = None
        # pass user-agent to allow rtd downloading from github
        self.headers = headers or {"User-Agent": "easyclimate"}  # Default User-Agent

    def __call__(self, url, output_file, pooch_instance):
        # Extract filename from URL
        filename = os.path.basename(url)

        def download_with_progress():
            with self.progress:
                self.task_id = self.progress.add_task(
                    "download", filename=filename, total=None
                )
                try:
                    # Create a new requests session
                    with requests.Session() as session:
                        response = session.get(url, stream=True)
                        response.raise_for_status()
                        total_size = int(response.headers.get("content-length", 0))

                        # Update the task with the total size if available
                        if total_size > 0:
                            self.progress.update(self.task_id, total=total_size)

                        # Write the file in chunks and update progress
                        downloaded = 0
                        with open(output_file, "wb") as f:
                            for chunk in response.iter_content(chunk_size=8192):
                                if chunk:
                                    f.write(chunk)
                                    downloaded += len(chunk)
                                    self.progress.update(
                                        self.task_id, advance=len(chunk)
                                    )

                        # Ensure the progress bar reaches 100% if the download is complete
                        if total_size > 0 and downloaded >= total_size:
                            self.progress.update(
                                self.task_id, completed=total_size, refresh=True
                            )
                finally:
                    self.task_id = None

        # Call the download function
        download_with_progress()
        return output_file


def _parse_md5_content(content: str) -> str:
    match = re.search(r"\b[0-9a-fA-F]{32}\b", content)
    if match is None:
        raise ValueError("MD5 file does not contain a valid checksum.")
    return f"md5:{match.group(0).lower()}"


def _get_known_hash(path: pathlib.Path, url: str, cache_dir) -> str | None:
    if url in external_hashes:
        return external_hashes[url]
    if path.name in tutorial_hashes:
        return tutorial_hashes[path.name]

    md5_name = path.with_suffix(".md5").name
    cache_key = (version, md5_name)
    if cache_key in _known_hash_cache:
        return _known_hash_cache[cache_key]

    md5_cache_path = pathlib.Path(cache_dir) / "md5" / md5_name
    if md5_cache_path.exists():
        known_hash = _parse_md5_content(md5_cache_path.read_text())
        _known_hash_cache[cache_key] = known_hash
        return known_hash

    md5_url = f"{base_url}/raw/{version}/{path.with_suffix('.md5').name}"
    warnings.warn(
        f"Tutorial dataset checksum for {path.name!r} is not recorded in "
        f"{_tutorial_config_path.name}; fetching {md5_url!r}. "
        "Run `python scripts/update_tutorial_config.py` to update the local "
        "tutorial dataset registry.",
        UserWarning,
        stacklevel=2,
    )
    response = requests.get(md5_url, headers={"User-Agent": "easyclimate"}, timeout=30)
    if response.status_code == 404:
        return None
    response.raise_for_status()
    known_hash = _parse_md5_content(response.text)

    try:
        md5_cache_path.parent.mkdir(parents=True, exist_ok=True)
        md5_cache_path.write_text(response.text)
    except OSError:
        pass
    _known_hash_cache[cache_key] = known_hash
    return known_hash


def _remove_cached_download(cache_dir, filename: str):
    for cached_file in pathlib.Path(cache_dir).glob(f"*-{filename}"):
        cached_file.unlink(missing_ok=True)


def _retrieve_with_retries(
    *,
    url: str,
    known_hash: str | None,
    path,
    progressbar: bool,
    downloader,
    filename: str,
):
    last_error = None
    for attempt in range(download_retries):
        try:
            return pooch.retrieve(
                url=url,
                known_hash=known_hash,
                path=path,
                progressbar=progressbar,
                downloader=downloader,
            )
        except (ValueError, requests.RequestException) as exc:
            last_error = exc
            _remove_cached_download(path, filename)
            if attempt == download_retries - 1:
                break
    raise RuntimeError(
        f"Failed to download a valid copy of {filename} after "
        f"{download_retries} attempts."
    ) from last_error


def open_tutorial_dataset(
    name: str,
    cache: bool = True,
    cache_dir: None | str | os.PathLike = None,
    progressbar: bool = False,
    *,
    engine: T_Engine = None,
    **kws,
) -> xr.Dataset:
    """
    Open a dataset from the online repository (requires internet).

    If a local copy is found then always use that to avoid network traffic.

    Available datasets:

    * ``"air_202201_mon_mean"``: 2m air temperature of the NCEP reanalysis subset
    * ``"hgt_202201_mon_mean"``: Geopotential height of the NCEP reanalysis subset
    * ``"precip_202201_mon_mean"``: Precipitation of the NCEP reanalysis subset
    * ``"pressfc_202201_mon_mean"``: Mean sea surface pressure of the NCEP reanalysis subset
    * ``"shum_202201_mon_mean"``: Absolute humidity of the NCEP reanalysis subset
    * ``"uwnd_202201_mon_mean"``: Zonal wind of the NCEP reanalysis subset
    * ``"vwnd_202201_mon_mean"``: Meridional wind of the NCEP reanalysis subset
    * ``"omega_202201_mon_mean"``: Vertical velocity of the NCEP reanalysis subset
    * ``"mini_HadISST_ice"``: Hadley Centre Sea Ice and Sea Surface Temperature data set (HadISST) subset
    * ``"PressQFF_202007271200_872"``: Observational data from European stations (from https://github.com/EXCITED-CO2/xarray-regrid)
    * ``"pr_wtr_eatm_2022"``: Precipitable water of the NCEP reanalysis subset in the 2022
    * ``"sst_mnmean_oisst"``: NOAA Optimum Interpolation (OI) SST V2 (from https://psl.noaa.gov/data/gridded/data.noaa.oisst.v2.html)

    Parameters
    ----------
    name : :py:class:`str <str>`
        Name of the file containing the dataset.
        e.g. 'air_202201_mon_mean'
    cache_dir : path-like, optional
        The directory in which to search for and write cached data.
    cache : dim: :py:class:`bool <bool>`, optional
        If True, then cache data locally for use on subsequent calls
    progressbar: :py:class:`bool <bool>`, default `False`.
        If True, will print a progress bar of the download to standard error (stderr).
    **kws : :py:class:`dict <dict>`, optional
        Passed to xarray.open_dataset

    Returns
    -------
    :py:class:`xarray.Dataset<xarray.Dataset>`

    Reference
    --------------
    - Kalnay et al.,The NCEP/NCAR 40-year reanalysis project, Bull. Amer. Meteor. Soc., 77, 437-470, 1996
    - Rayner, N. A.; Parker, D. E.; Horton, E. B.; Folland, C. K.; Alexander, L. V.; Rowell, D. P.; Kent, E. C.; Kaplan, A. (2003) Global analyses of sea surface temperature, sea ice, and night marine air temperature since the late nineteenth century J. Geophys. Res.Vol. 108, No. D14, 4407 10.1029/2002JD002670  (pdf ~9Mb)

    .. seealso::
        - :py:func:`xarray.tutorial.load_dataset<xarray.tutorial.load_dataset>`
        - :py:func:`xarray.open_dataset<xarray.open_dataset>`
        - :py:func:`xarray.load_dataset<xarray.load_dataset>`
    """
    logger = pooch.get_logger()
    logger.setLevel("WARNING")

    # pass user-agent to allow rtd downloading from github
    # downloader = pooch.HTTPDownloader(headers={"User-Agent": "easyclimate"})
    downloader = RichDownloader()

    cache_dir = _construct_cache_dir(cache_dir)
    known_hash = None
    if name in external_urls:
        url = external_urls[name]
        path = pathlib.Path(url)
        known_hash = external_hashes.get(name) or external_hashes.get(url)
    else:
        path = pathlib.Path(name)
        if not path.suffix:
            # process the name
            default_extension = ".nc"
            if engine is None:
                _check_netcdf_engine_installed(name)
            path = path.with_suffix(default_extension)
        elif path.suffix == ".grib":
            if engine is None:
                engine = "cfgrib"
                try:
                    import cfgrib  # noqa
                except ImportError as e:
                    raise ImportError(
                        "Reading this tutorial dataset requires the cfgrib package."
                    ) from e

        url = f"{base_url}/raw/{version}/{path.name}"
        known_hash = _get_known_hash(path, url, cache_dir)

    # retrieve the file
    filepath = _retrieve_with_retries(
        url=url,
        known_hash=known_hash,
        path=cache_dir,
        progressbar=progressbar,
        downloader=downloader,
        filename=path.name,
    )

    if Path(filepath).suffix == ".nc" or Path(filepath).suffix == ".grib":
        # print(f"The sample data is located in {filepath} with {engine}")
        ds = xr.open_dataset(filepath, engine=engine, **kws)
        if not cache:
            ds = ds.load()
            pathlib.Path(filepath).unlink()
    elif Path(filepath).suffix == ".csv" or Path(filepath).suffix == ".CSV":
        ds = pd.read_csv(filepath)

    return ds
