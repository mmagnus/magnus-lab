# Magnus Lab — website

Website of the Laboratory of RNA Design and Therapeutics, Centre of New Technologies,
University of Warsaw (PI: dr Marcin Magnus). Live at https://magnus-lab.org

Static site: one HTML file, no build step, no framework. All styling is inline in
`index.html`; the only external request is the Google Fonts stylesheet.

## Files

```
index.html                       the whole site
CNAME                            custom domain for GitHub Pages
favicon.png                      lab mark
assets/logo-lab.png              lab mark, header
assets/rna-diffusion.gif         hero animation (ligand-conditioned RNA diffusion)
assets/marcin-magnus.jpg         PI photo
assets/jadwiga-meissner.jpg      postdoc photo
assets/logo-*.png|jpg            UW, CeNT, NCN, Harvard, UWA, RNA Club, Do Science!
tools/make_rna_diffusion_gif.py  script that generates the hero animation
tools/bases.py                   atom coordinates used by that script
.nojekyll                        serve the files as they are
```

## Publishing

GitHub Pages, branch `main`, folder `/ (root)`. The domain `magnus-lab.org` points at
GitHub's four A records (185.199.108–111.153); `www` is a CNAME to `mmagnus.github.io`.

## Editing

Plain HTML with readable section ids: `#research`, `#news`, `#team`, `#collaborators`,
`#publications`, `#teaching`, `#sig`, `#join`, `#contact`. Colours, fonts and spacing are
CSS custom properties on `:root` near the top of the file — change `--accent` to restyle
the whole page.

The hero animation is generated from real coordinates of PDB entry
[1DDY](https://www.rcsb.org/structure/1DDY):

```bash
pip install pillow
python tools/make_rna_diffusion_gif.py
```

## Still to complete

The author list for the OpenRNAFold paper. The Günter Mayer collaborator card is in the
source but commented out.
