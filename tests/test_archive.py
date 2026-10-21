import zipfile

from fus import archive


def _zip(tmp_path, names):
    path = tmp_path / "US_Data.zip"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_STORED) as z:
        for n in names:
            z.writestr(n, b"x" * 10)
    return path


def test_recordings_lists_targets_only(tmp_path):
    path = _zip(tmp_path, [
        "US_Data/061126_FUS_AAV_Controller/Mouse_Cntr_01_Target2_Repeat.mat",
        "US_Data/070126_FUS_AAV_Controller/Mouse_Cntr_17_Target5Repeat.mat",
        "US_Data/070126_FUS_AAV_Controller/Mouse_Cntr_13_Baseline.mat",
        "US_Data/070126_FUS_AAV_Controller/Reduced_Mouse_Cntr_13_Target1.mat",
        "US_Data/061626_FUS_AAV_Controller/Mouse_Cntr_07_Target1.mat",
    ])
    recs = archive.recordings(path)
    assert [(r.day, r.mouse, r.recording) for r in recs] == [
        ("20260611", 1, "Target2_Repeat"), ("20260616", 7, "Target1"), ("20260701", 17, "Target5Repeat")]


def test_local_path_temporary_copy_is_removed(tmp_path, monkeypatch):
    monkeypatch.setattr(archive, "LOCAL", tmp_path / "none")
    path = _zip(tmp_path, ["US_Data/070726_FUS_AAV_Controller/Mouse_Cntr_24_Target3.mat"])
    (m,) = archive.recordings(path)
    with archive.local_path(m, path) as p:
        assert p.read_bytes() == b"x" * 10
    assert not p.exists()
