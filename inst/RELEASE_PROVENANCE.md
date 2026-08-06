# Release provenance

`PACKAGE_SOURCE_MANIFEST.csv` records portable SHA-256 hashes for the R, C++,
test, vignette, authors-method script, validation driver, help, workflow, and
release-documentation files in the pre-build repository source snapshot. It
also records repository build controls that `R CMD build` may exclude. Each
rebuttal extdata directory has a separate package-subset manifest.

`R CMD build` canonicalizes package metadata and generates installed vignette
artifacts, so the manifest is not presented as a byte inventory of the tarball.
The exact built source tarball is accompanied by an external `.sha256` file
because a tarball cannot contain its own final hash.

The current workspace is not an initialized Git repository, so this snapshot
does not yet have a repository URL, commit, or immutable tag. Before a public
GitHub release, the authors should select the license, replace the `OWNER`
installation placeholders, add `URL` and `BugReports` to `DESCRIPTION`, create
an immutable tag, and attach the source tarball plus its checksum.

## Historical rebuttal source boundary

The seven-method tensor-regression archive records the SHA-256 of the exact
historical runner and solver. The current active solver was edited after that
freeze, and an exact copy of the historical solver is not present in this
workspace. Therefore SegreCP makes two narrower, testable claims:

1. the rebuttal table and figures are deterministically reconstructed from the
   paired frozen raw records and their manifests; and
2. the public RGD/RGN implementation is independently checked against direct
   numerical oracles, while the archived formal version-0.1.0
   decomposition-initialization rerun provides separate historical evidence
   under its recorded protocol.

It does not claim that the complete seven-method raw archive can be regenerated
bit-for-bit from a silently substituted later solver file. An immutable public
revision containing the exact historical source would be needed for that
stronger claim.

## Version-0.1.1 code-scope split

The installed package contains only the authors' RGD/RGN implementation and
initialization layer. The current comparator entry point lives under the
repository-only `reproducibility/` directory, which is excluded from package
builds. The mixed initialization runner used before this split is preserved
byte-for-byte inside
`reproducibility/archive/SegreCP_0.1.0_pre_scope_split.tar.gz`; its tar member
hash matches the `fa2d3ef0...` value already recorded by the frozen
`decomp_pkg/source_manifest.csv`. The frozen bundle was not rewritten to make
it appear that a later runner generated the original results.
