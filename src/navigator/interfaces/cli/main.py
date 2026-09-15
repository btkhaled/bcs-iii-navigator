"""CLI predict."""
from __future__ import annotations

import json
import sys


def main():
    from ...application.predict import predict
    from ...engine.report import summarize

    if len(sys.argv) < 4:
        print("usage: predict <compound.json> <formulation.json> <config_dir>")
        raise SystemExit(2)
    r = predict(sys.argv[1], sys.argv[2], sys.argv[3])
    print(json.dumps(summarize(r), indent=2))


if __name__ == "__main__":
    main()
