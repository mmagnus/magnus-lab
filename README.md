# Magnus Lab — website

Website of the Laboratory of RNA Design and Therapeutics, Centre of New Technologies,
University of Warsaw (PI: dr Marcin Magnus).

Static site: one HTML file, no build step, no JavaScript framework. All styling is inline
in `index.html`; the only external request is the Google Fonts stylesheet.

## Files

```
index.html                       the whole site
assets/rna-diffusion.gif         hero animation (ligand-conditioned RNA diffusion)
assets/marcin-magnus.jpg         PI photo
assets/jadwiga-meissner.jpg      postdoc photo
assets/logo-*.png|jpg            UW, CeNT, NCN, Harvard, UWA, RNA Club, Do Science!
tools/make_rna_diffusion_gif.py  script that generates the hero animation
tools/bases.py                   atom coordinates used by that script
.nojekyll                        tells GitHub Pages to serve the files as they are
```

## Publishing on GitHub Pages

1. Push this directory to a repository, e.g. `mmagnus/magnus-lab`.
2. Settings → Pages → Source: *Deploy from a branch*, branch `main`, folder `/ (root)`.
3. The site appears at `https://<user>.github.io/<repo>/`.

For a custom domain (e.g. `magnuslab.uw.edu.pl`), add a `CNAME` file containing the domain
and point a DNS CNAME record at `<user>.github.io`.

## Editing

Everything is plain HTML with readable section ids: `#research`, `#news`, `#team`,
`#collaborators`, `#publications`, `#teaching`, `#sig`, `#join`, `#contact`. Colours, fonts
and spacing are CSS custom properties on `:root` near the top of the file — change
`--accent` to restyle the whole page.

The hero animation is generated from real coordinates of PDB entry
[1DDY](https://www.rcsb.org/structure/1DDY). To regenerate it:

```bash
pip install pillow
python tools/make_rna_diffusion_gif.py
```

## Notes

Two items are still marked as to be completed on the page: the author list for the
OpenRNAFold paper, and the venue/date details of positions that have not opened yet.
The Günter Mayer collaborator card is present in the source but commented out.
