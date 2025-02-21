project = "Swimming Particle"
author = "Naoki Hori"
copyright = f"2025, {author}"

html_theme = "alabaster"
html_theme_options = {
    "fixed_sidebar": "false",
    "github_banner": "false",
    "github_button": "true",
    "github_count": "true",
    "github_repo": "SwimmingParticle",
    "github_type": "star",
    "github_user": "NaokiHori",
    "navigation_with_keys": "true",
    "nosidebar": "false",
}

mathjax_path = "https://cdn.jsdelivr.net/npm/mathjax@2/MathJax.js?config=TeX-AMS-MML_HTMLorMML"

macros = dict()
macros["vr"] = ["r", 0]
macros["vt"] = ["\\theta", 0]
macros["ur"] = ["u_{\\vr}", 0]
macros["ut"] = ["u_{\\vt}", 0]
macros["nr"] = ["N_{\\vr}", 0]
macros["nt"] = ["N_{\\vt}", 0]
macros["vat"] = ["\\left. #1 \\right|_{#2}", 2]
macros["pder"] = ["\\frac{\\partial #1}{\\partial #2}", 2]
macros["dder"] = ["\\frac{\\delta #1}{\\delta #2}", 2]

# discrete advective terms in scalar transport
macros["ave"] = ["\\overline{#1}^{#2}", 2]
macros["dif"] = ["{\\delta_{#2} {#1}}", 2]
macros["gcs"] = ["{\\xi^{#1}}", 1]
macros["sfact"] = ["h_{\\gcs{#1}}", 1]
macros["vel"] = ["{u_{#1}}", 1]
macros["scalar"] = ["c", 0]
macros["dscalaradv"] = [
        "\\frac{1}{J}"
        "\\ave{"
        "  \\frac{J}{\\sfact{#1}}"
        "  \\vel{#1}"
        "  \\dif{#2}{\\gcs{#1}}"
        "}{\\gcs{#1}}"
        , 2
]
# discrete diffusive terms in scalar transport
macros["dscalardif"] = [
        "\\frac{1}{J}"
        "\\dif{"
        "}{\\gcs{#1}}"
        "\\left("
        "  \\frac{J}{\\sfact{#1}}"
        "  \\frac{1}{\\sfact{#1}}"
        "  \\dif{#2}{\\gcs{#1}}"
        "\\right)"
        , 2
]

mathjax3_config = {"TeX": {"Macros": macros}}
