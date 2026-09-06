"""Check the versioned release's input identity and file inventory."""

import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent


def verify(root=HERE):
    manifest = json.loads((root / "manifest.json").read_text())
    files = {p.relative_to(root).as_posix(): p for p in root.rglob("*") if p.is_file()
             and not {"build", "__pycache__", ".pytest_cache"}.intersection(p.relative_to(root).parts)
             and p.name != "manifest.json"}
    declared = {record["file"]: record for record in manifest["files"]}
    if set(files) != set(declared):
        raise ValueError(f"File inventory mismatch: {sorted(set(files) ^ set(declared))}")
    for name, path in files.items():
        payload = path.read_bytes()
        if len(payload) != declared[name]["bytes"] or hashlib.sha256(payload).hexdigest() != declared[name]["sha256"]:
            raise ValueError(f"Changed release input: {name}")
        if len(payload) > 1_000_000 or b"\x00" in payload:
            raise ValueError(f"Unexpected binary or oversized input: {name}")
        payload.decode("utf-8")
    return {"status": "PASS", "files": len(files),
            "bytes": sum(path.stat().st_size for path in files.values())}


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
