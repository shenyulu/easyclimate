from pathlib import Path
import shutil
import uuid

import pytest

from easyclimate.core import tutorial


def _local_tmpdir():
    path = Path(".pytest_tmp") / f"tutorial-{uuid.uuid4().hex}"
    path.mkdir(parents=True, exist_ok=False)
    return path


def test_parse_md5_content_accepts_common_formats():
    assert (
        tutorial._parse_md5_content(
            "d771ab1a6e5ced52b3082e4b9afbf947  vwnd_202201_mon_mean.nc"
        )
        == "md5:d771ab1a6e5ced52b3082e4b9afbf947"
    )


def test_retrieve_with_retries_redownloads_bad_cache(monkeypatch):
    filename = "vwnd_202201_mon_mean.nc"
    cache_dir = _local_tmpdir()
    calls = []

    def fake_retrieve(**kwargs):
        calls.append(kwargs)
        cached_file = cache_dir / f"cached-{filename}"
        cached_file.write_text("bad download")
        if len(calls) < 3:
            raise ValueError("MD5 hash mismatch")
        return str(cached_file)

    monkeypatch.setattr(tutorial.pooch, "retrieve", fake_retrieve)

    try:
        result = tutorial._retrieve_with_retries(
            url="https://example.test/vwnd_202201_mon_mean.nc",
            known_hash="md5:d771ab1a6e5ced52b3082e4b9afbf947",
            path=cache_dir,
            progressbar=False,
            downloader=None,
            filename=filename,
        )

        assert Path(result).exists()
        assert len(calls) == 3
    finally:
        shutil.rmtree(cache_dir, ignore_errors=True)


def test_retrieve_with_retries_fails_after_three_bad_downloads(monkeypatch):
    filename = "vwnd_202201_mon_mean.nc"
    cache_dir = _local_tmpdir()
    calls = []

    def fake_retrieve(**kwargs):
        calls.append(kwargs)
        (cache_dir / f"cached-{filename}").write_text("bad download")
        raise ValueError("MD5 hash mismatch")

    monkeypatch.setattr(tutorial.pooch, "retrieve", fake_retrieve)

    try:
        with pytest.raises(RuntimeError, match="after 3 attempts"):
            tutorial._retrieve_with_retries(
                url="https://example.test/vwnd_202201_mon_mean.nc",
                known_hash="md5:d771ab1a6e5ced52b3082e4b9afbf947",
                path=cache_dir,
                progressbar=False,
                downloader=None,
                filename=filename,
            )

        assert not list(cache_dir.glob(f"*-{filename}"))
        assert len(calls) == 3
    finally:
        shutil.rmtree(cache_dir, ignore_errors=True)
