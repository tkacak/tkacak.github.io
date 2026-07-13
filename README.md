# tkacak.github.io

Personal academic website of Tugay Kaçak — <https://tkacak.github.io/>

Built with **R Markdown websites** (`rmarkdown::render_site()`), the same setup used by
several psychometrics researchers. Source lives in the `.Rmd` files; the rendered `.html`
files are committed to the repo root so GitHub Pages can serve them.

## Structure

| File | Page |
|:--|:--|
| `_site.yml` | Site config: navbar, theme, CSS |
| `index.Rmd` | Home |
| `research.Rmd` | Research lines and education |
| `publications.Rmd` | Journal articles, chapters, work in progress |
| `talks.Rmd` | Conference papers and posters |
| `software.Rmd` | Software, materials, teaching |
| `style.css` | Custom styling on top of the `cosmo` Bootstrap theme |
| `header.html` | SEO meta tags and schema.org markup |
| `images/` | Profile photo and other images |

## How to build

Once, in R:

```r
install.packages("rmarkdown")
```

Then, from the project folder:

```r
rmarkdown::render_site()
```

This writes `index.html`, `research.html`, `publications.html`, `talks.html` and
`software.html` into the repo root (`output_dir: "."`). Commit **both** the `.Rmd` sources
and the generated `.html` files, then push:

```bash
git add -A
git commit -m "Update site"
git push
```

GitHub Pages (Settings → Pages → branch `main`, folder `/root`) serves the result.

## Adding a publication

Open `publications.Rmd`, add one line at the top of the right year block, bump the
superscript number, re-render, push. That's it.

## License

Content © Tugay Kaçak. Code released under the MIT License.
