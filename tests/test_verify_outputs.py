"""Exercise the verifier's PNG tolerance and exact content rejection end to end."""
import importlib.util
from pathlib import Path
import shutil

from PIL import Image, PngImagePlugin
import pytest

spec = importlib.util.spec_from_file_location(
    "verify_outputs", Path(__file__).resolve().parents[1] / "scripts/verify_outputs.py"
)
verify = importlib.util.module_from_spec(spec)
spec.loader.exec_module(verify)


@pytest.fixture
def outputs(tmp_path, monkeypatch):
    expected = tmp_path / "results"
    expected.mkdir()
    for name in verify.OUTPUTS[:-1]:
        (expected / name).write_bytes(b"exact text\n")
    Image.new("RGBA", (2, 3), (10, 20, 30, 255)).save(expected / "summary.png")
    monkeypatch.setattr(verify, "ROOT", tmp_path)

    def regenerate(_input, _resource, destination):
        for path in expected.iterdir():
            shutil.copyfile(path, destination / path.name)
        return destination

    return regenerate


@pytest.mark.parametrize("level", [0, 9])
def test_png_container_changes_match(outputs, monkeypatch, capsys, level):
    def regenerate(*args):
        destination = outputs(*args)
        path = destination / "summary.png"
        metadata = PngImagePlugin.PngInfo()
        metadata.add_text("Software", "different encoder")
        with Image.open(path) as image:
            image.load()
            image.save(path, pnginfo=metadata, compress_level=level)
        assert path.read_bytes() != (verify.ROOT / "results/summary.png").read_bytes()

    monkeypatch.setattr(verify, "run", regenerate)
    verify.main()
    assert "Verified 5/5 reproducible outputs" in capsys.readouterr().out


@pytest.mark.parametrize("change", ["pixel", "mode", "size"])
def test_png_content_changes_fail(outputs, monkeypatch, change):
    def regenerate(*args):
        path = outputs(*args) / "summary.png"
        with Image.open(path) as original:
            image = original.copy()
        if change == "pixel":
            image.putpixel((0, 0), (10, 20, 31, 255))
        elif change == "mode":
            image = image.convert("RGB")
        else:
            # Identical pixel bytes, different dimensions must still fail.
            image = Image.frombytes(image.mode, (3, 2), image.tobytes())
        image.save(path)

    monkeypatch.setattr(verify, "run", regenerate)
    with pytest.raises(SystemExit, match="^Reproducibility mismatch: summary.png$"):
        verify.main()


@pytest.mark.parametrize("name", verify.OUTPUTS[:-1])
def test_text_byte_changes_fail(outputs, monkeypatch, name):
    def regenerate(*args):
        path = outputs(*args) / name
        path.write_bytes(path.read_bytes() + b"\n")

    monkeypatch.setattr(verify, "run", regenerate)
    with pytest.raises(SystemExit, match=f"^Reproducibility mismatch: {name}$"):
        verify.main()
