# Ontology processing invariants

This note preserves implementation history that is useful for preventing
scientific regressions but is too detailed for production docstrings.

## Parental-background filter

The parental test asks whether a domain–child association is enriched within
the proteins annotated to the parent. The term index must therefore implement
the True Path Rule: a protein annotated to a child also belongs to every
ancestor's background. Domain membership remains direct because domains have no
hierarchy in this stage.

An earlier implementation indexed only direct term annotations. Parents without
direct annotations then appeared to have empty backgrounds, causing 54,951
human t0 associations to be rejected without a test. Only about 14% of
associations survived, which invalidated the affected True Path ablation.
Any ablation result produced with the unpropagated background must not be cited
or compared as scientific evidence; rerun it with the corrected index. Keep the
propagated-background regression tests when changing this code.

The maps are inverted once per filtering run. Re-scanning the proteome for each
domain–child–parent test produces the same contingency cells but is prohibitive
at production scale (roughly 10^10 Python membership checks for a representative
human run).

## Failure policy

An undersized background and an invalid contingency table are separate,
explicit data-dependent rejection categories. Each is counted and the
association is conservatively rejected. Unexpected exceptions are programming
failures and must propagate instead of silently removing results.

Graph metrics are calculated during ``OntologyProcessor`` initialization.
Failures abort initialization rather than leaving a processor whose graph
summary is partial or potentially misleading.
