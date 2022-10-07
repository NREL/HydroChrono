(TeX-add-style-hook
 "report"
 (lambda ()
   (TeX-add-to-alist 'LaTeX-provided-class-options
                     '(("article" "11pt")))
   (TeX-add-to-alist 'LaTeX-provided-package-options
                     '(("inputenc" "utf8") ("fontenc" "T1") ("ulem" "normalem") ("helvet" "scaled")))
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "path")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "url")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "nolinkurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperbaseurl")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperimage")
   (add-to-list 'LaTeX-verbatim-macros-with-braces-local "hyperref")
   (add-to-list 'LaTeX-verbatim-macros-with-delims-local "path")
   (TeX-run-style-hooks
    "latex2e"
    "article"
    "art11"
    "inputenc"
    "fontenc"
    "graphicx"
    "grffile"
    "longtable"
    "wrapfig"
    "rotating"
    "ulem"
    "amsmath"
    "textcomp"
    "amssymb"
    "capt-of"
    "hyperref"
    "mathpazo"
    "helvet"
    "courier")
   (LaTeX-add-labels
    "sec:org5a1202f"
    "fig:org8863164"
    "tab:orgbce3e80"
    "sec:org495e0d0"
    "sec:org7227b48"
    "sec:org89ffc0f"
    "sec:orgdfa82f3"
    "sec:org731c339"
    "sec:orgc32ec33"
    "fig:org99e9384"
    "sec:orge230bac"
    "fig:org6aa40c0"
    "sec:org8c1320b"
    "fig:org54dce10"
    "sec:org0e95c37"
    "sec:orgcf6146f"
    "sec:orge5c01d0"
    "sec:org6f61678"
    "sec:org39f9777"
    "sec:orgd3c8058"
    "sec:org1b17e99"
    "sec:org388de5b"
    "sec:org04f3989"))
 :latex)

