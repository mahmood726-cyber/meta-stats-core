# sentinel-findings.md

*Written by Sentinel — WARN-tier findings.*

## [WARN] P1-claim-grounding
- **Location:** `README.md:48`
- **Detail:** document states 1 quantitative effect claim(s) (e.g. 'p<0.1') but carries no resolvable source identifier (DOI/PMID/NCT/URL) anywhere - the claims are ungrounded.
- **Fix hint:** add the DOI/PMID/NCT of the source each effect estimate comes from, or run the Overmind claim-grounding witness (overmind ground) to bind each claim to a corpus record
- **Source:** F:\e156\docs\assurance-standard.md#1-citation-verification
- **When:** 2026-06-13T08:04:20.808365+00:00
