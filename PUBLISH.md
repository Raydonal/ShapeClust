# Publishing this repository

Everything below assumes you are at the root of this directory.

## 1. Create the repository and push

```bash
git init
git add .
git commit -m "ShapeClust 0.2.0: package, manuscript, and reproduction scripts"
git branch -M main
git remote add origin https://github.com/raydonal/ShapeClust.git
git push -u origin main
```

## 2. Turn on the website

On GitHub: **Settings → Pages → Source: Deploy from a branch → `gh-pages` / root**.

The `pkgdown` workflow creates the `gh-pages` branch on the first push to `main`,
so the branch will only appear in the dropdown after that first Actions run
finishes. The site then lives at <https://raydonal.github.io/ShapeClust/>.

## 3. What the workflows do

| Workflow | Trigger | What it does |
|---|---|---|
| `R-CMD-check.yaml` | push / PR to `main` | `R CMD check` on macOS, Windows and Ubuntu, over R release, devel and oldrel |
| `pkgdown.yaml` | push to `main` | builds the documentation site and deploys it to `gh-pages` |

Both are standard `r-lib/actions` recipes and need no secrets: they use the
`GITHUB_TOKEN` that Actions provides automatically.

## 4. Two notes on `R CMD check`

The package currently passes with **0 errors, 0 warnings and 3 notes**. Two of
the three disappear once the repository is public:

- *CRAN incoming feasibility* — flags the GitHub URLs in `DESCRIPTION` as
  unreachable, because they 404 until the repository exists.
- *future file timestamps* — an artifact of building inside a container.
- *compilation flags* — reports the flags of the build machine, not of the
  package.

## 5. Before submitting to CRAN

- Add `cran-comments.md` describing the submission.
- Check on the winbuilder and macbuilder services.
- Bump the version in `DESCRIPTION` and add a `NEWS.md` entry.
- The vignette must build in under ten minutes on CRAN's machines; the
  reproduction scripts under `inst/` are *not* run at build time, so they do not
  affect this.

## 6. Archiving for the paper

Once the manuscript is accepted, mint a DOI by connecting the repository to
Zenodo (Zenodo → GitHub → toggle the repository on, then cut a release). Then
replace the placeholder sentence in the *Code and data availability* section of
`paper/paper.tex`, which currently reads that a permanent archived version will
be deposited on acceptance.
