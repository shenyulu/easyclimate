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


def test_get_known_hash_uses_cached_md5_without_http(monkeypatch):
    cache_dir = _local_tmpdir()
    md5_dir = cache_dir / "md5"
    md5_dir.mkdir()
    (md5_dir / "vwnd_202201_mon_mean.md5").write_text(
        "d771ab1a6e5ced52b3082e4b9afbf947  vwnd_202201_mon_mean.nc"
    )

    def fail_get(*args, **kwargs):
        raise AssertionError("md5 cache should avoid HTTP")

    monkeypatch.setattr(tutorial.requests, "get", fail_get)
    tutorial._known_hash_cache.clear()

    try:
        known_hash = tutorial._get_known_hash(
            Path("vwnd_202201_mon_mean.nc"),
            "https://example.test/vwnd_202201_mon_mean.nc",
            cache_dir,
        )

        assert known_hash == "md5:d771ab1a6e5ced52b3082e4b9afbf947"
    finally:
        tutorial._known_hash_cache.clear()
        shutil.rmtree(cache_dir, ignore_errors=True)


def test_get_known_hash_uses_builtin_hash_before_http(monkeypatch):
    cache_dir = _local_tmpdir()

    def fail_get(*args, **kwargs):
        raise AssertionError("built-in hash should avoid HTTP")

    monkeypatch.setattr(tutorial.requests, "get", fail_get)
    tutorial._known_hash_cache.clear()

    try:
        known_hash = tutorial._get_known_hash(
            Path("uwnd_vwnd_hgt_equtorial_2021_2024.nc"),
            "https://example.test/uwnd_vwnd_hgt_equtorial_2021_2024.nc",
            cache_dir,
        )

        assert known_hash == "md5:e313295bbb8fb9631f65f6f122513ea1"
    finally:
        tutorial._known_hash_cache.clear()
        shutil.rmtree(cache_dir, ignore_errors=True)


def test_get_known_hash_writes_remote_md5_to_cache(monkeypatch):
    cache_dir = _local_tmpdir()
    calls = []

    class Response:
        status_code = 200
        text = "d771ab1a6e5ced52b3082e4b9afbf947  custom_dataset.nc"

        def raise_for_status(self):
            return None

    def fake_get(*args, **kwargs):
        calls.append((args, kwargs))
        return Response()

    monkeypatch.setattr(tutorial.requests, "get", fake_get)
    tutorial._known_hash_cache.clear()

    try:
        with pytest.warns(UserWarning, match="update_tutorial_config.py"):
            known_hash = tutorial._get_known_hash(
                Path("custom_dataset.nc"),
                "https://example.test/custom_dataset.nc",
                cache_dir,
            )
        second_known_hash = tutorial._get_known_hash(
            Path("custom_dataset.nc"),
            "https://example.test/custom_dataset.nc",
            cache_dir,
        )

        assert known_hash == "md5:d771ab1a6e5ced52b3082e4b9afbf947"
        assert second_known_hash == known_hash
        assert (cache_dir / "md5" / "custom_dataset.md5").exists()
        assert len(calls) == 1
    finally:
        tutorial._known_hash_cache.clear()
        shutil.rmtree(cache_dir, ignore_errors=True)
