# Configuration file of the nf-core/scalpel documentation (Sphinx, Read the Docs)
import re
from datetime import date

project = "nf-core/scalpel"
author = "Franz AKE"
copyright = f"{date.today().year}, {author}"
release = "1.0.0dev"

extensions = [
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
]
source_suffix = {".md": "markdown"}
master_doc = "index"
#- files of the repository docs/ that are not pages of the site
exclude_patterns = ["_build", "README.md", "vignettes/**", "Thumbs.db", ".DS_Store", "build_tutorials.sh"]

myst_enable_extensions = ["colon_fence", "deflist", "html_image", "html_admonition", "attrs_inline", "attrs_block"]
myst_heading_anchors = 3

# HTML output
# ++++++++++++
html_theme = "pydata_sphinx_theme"
html_title = "SCALPEL"
html_logo = "images/nf-core-scalpel_logo_light.png"
html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_show_sphinx = False
html_show_sourcelink = False
html_context = {
    "github_user": "plasslab",
    "github_repo": "nf-core-scalpel",
    "github_version": "dev",
    "doc_path": "docs",
    "default_mode": "light",
}
html_theme_options = {
    "logo": {"image_light": "images/nf-core-scalpel_logo_light.png", "image_dark": "images/nf-core-scalpel_logo_dark.png"},
    "icon_links": [
        {"name": "GitHub", "url": "https://github.com/plasslab/nf-core-scalpel", "icon": "fa-brands fa-github"},
        {"name": "scalpelR", "url": "https://github.com/Franzx7/scalpelR", "icon": "fa-brands fa-r-project"},
        {"name": "Paper", "url": "https://doi.org/10.1038/s41467-025-61118-0", "icon": "fa-solid fa-file-lines"},
    ],
    "navbar_align": "left",
    "header_links_before_dropdown": 8,
    "navbar_end": ["theme-switcher", "navbar-icon-links"],
    "navigation_depth": 3,
    "show_toc_level": 2,
    "secondary_sidebar_items": ["page-toc", "edit-this-page"],
    "article_header_end": [],
    "use_edit_page_button": True,
    "footer_start": ["copyright"],
    "footer_end": [],
}
html_sidebars = {"index": [], "changelog": []}


# GitHub alerts (> [!NOTE]) of the pipeline docs -> MyST admonitions
# +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
_ALERT = re.compile(r"^> \[!(NOTE|TIP|IMPORTANT|WARNING|CAUTION)\]\n((?:>.*\n?)*)", re.MULTILINE)


def _github_alerts(app, docname, source):
    def repl(m):
        body = "\n".join(re.sub(r"^> ?", "", line) for line in m.group(2).splitlines())
        return f":::{{{m.group(1).lower()}}}\n{body}\n:::\n"

    text = _ALERT.sub(repl, source[0])
    #- nf-core website notices & title prefix of the pipeline docs...
    text = re.sub(r"^## :warning: Please read this documentation on the nf-core website.*\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^> _Documentation of pipeline parameters is generated automatically.*\n", "", text, flags=re.MULTILINE)
    text = re.sub(r"^# nf-core/scalpel: ", "# ", text, count=1, flags=re.MULTILINE)
    #- links between the tutorials (standalone HTML pages) -> pages of the site
    text = re.sub(r"\]\(((?:chromium|dropseq)(?:_downstream)?)\.html\)", r"](\1.md)", text)
    #- no Pygments lexer for csv / tsv: plain text
    text = re.sub(r"^```(csv|tsv)\b", "```text", text, flags=re.MULTILINE)
    source[0] = text


def setup(app):
    app.connect("source-read", _github_alerts)
