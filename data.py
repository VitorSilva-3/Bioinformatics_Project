

ec_manual = {
    "beta-fructofuranosidase":             "3.2.1.26",
    "beta-galactosidase":                  "3.2.1.23",
    "beta-d-galactosidase":                "3.2.1.23",
    "lactase-phlorizin hydrolase":         "3.2.1.62",
    "arabinanase/levansucrase/invertase":  "3.2.1.26",
    "fructokinase":                        "2.7.1.4",
    "maltase":                             "3.2.1.20",
    "beta-amylase":                        "3.2.1.2",
    "alpha-amylase":                       "3.2.1.1",
    "cellulase":                           "3.2.1.4"
}

sugar_manual = {
    "beta-fructofuranosidase":             "sacarose",
    "beta-galactosidase":                  "lactose",
    "beta-d-galactosidase":                "lactose",
    "lactase-phlorizin hydrolase":         "lactose",
    "arabinanase/levansucrase/invertase":  "sacarose",
    "fructokinase":                        "fructose",
    "maltase":                             "maltose",
    "beta-amylase":                        "starch",
    "alpha-amylase":                       "starch",
    "cellulase":                           "cellulose"
}

queries = [
    '("beta-galactosidase"[All Fields] OR "lactase"[All Fields])',
    '("beta-fructofuranosidase"[All Fields] OR "invertase"[All Fields] OR "sucrase"[All Fields])',
    '("fructokinase"[All Fields])',
    '("maltase"[All Fields])',
    '("beta-amylase"[All Fields])',
    '("alpha-amylase"[All Fields])',
    '("cellulase"[All Fields])'
    ]

taxa = ["Chlorophyta", "Rhodophyta", "Glaucophyta", "Bacillariophyta", "Haptophyta"]

enzymes = {
    "beta-galactosidase":         "3.2.1.23",
    "lactase":                    "3.2.1.23",
    "beta-fructofuranosidase":    "3.2.1.26",
    "invertase":                  "3.2.1.26",
    "sucrase":                    "3.2.1.26",
    "fructokinase":               "2.7.1.4",
    "maltase":                    "3.2.1.20",
    "beta-amylase":               "3.2.1.2",
    "alpha-amylase":              "3.2.1.1",
    "cellulase":                  "3.2.1.4"
}