# Champion SKILL.md artifacts

These are frozen outputs from clawpathy-autoresearch runs. Each is the best-scoring SKILL.md the loop found for its task. Dropped here as a public reference so others can download a single champion without running the tuner.

| File | Task | Final judge score (lower=better) |
|------|------|------|
| `trubetskoy_scz_finemap_0.235.md` | Trubetskoy et al. 2022 PGC3 SCZ fine-mapping (FINEMAP v1.4 + 1000G EUR LD) | 0.235 |
| `yengo_height_ldsc_h2.md` | Yengo et al. 2022 height LDSC h² reproduction | ~0.12 (target 0.15 beaten) |

## Fetch one without cloning the repo

```python
from pathlib import Path
import urllib.request

def download_champion(name: str = "trubetskoy_scz_finemap_0.235",
                      dest: str | Path = ".",
                      ref: str = "main") -> Path:
    """Download a champion SKILL.md from the ClawBio repo."""
    url = (f"https://raw.githubusercontent.com/ClawBio/ClawBio/{ref}"
           f"/skills/clawpathy_autoresearch/examples/champions/{name}.md")
    out = Path(dest) / f"{name}.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url) as r:
        out.write_bytes(r.read())
    return out
```

Shell:
```bash
curl -sL -o trubetskoy.md \
  https://raw.githubusercontent.com/ClawBio/ClawBio/main/skills/clawpathy_autoresearch/examples/champions/trubetskoy_scz_finemap_0.235.md
```
