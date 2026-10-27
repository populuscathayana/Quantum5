"""Offline helpers; no cluster access or job execution."""
import hashlib
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def save_json(path,data):
    with Path(path).open('x') as stream:json.dump(data,stream,indent=2,ensure_ascii=False)
