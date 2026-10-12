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
GitHub release, the authors should select the license, add a verified
repository address to the installation guide plus `URL` and `BugReports` to
`DESCRIPTION`, create
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

## Version-0.1.2 real-data boundary

The repository-only `reproducibility/real_data/` workflow contains source
locations, preprocessing code, and input checks, but no raw or processed trade
data. It is excluded from the installed package. The Markdown vignette documents
how to rebuild the selected inputs in a private directory and run the authors'
methods from the package. A fresh 0.1.2 RGD/RGN replay on the four frozen
processed inputs matched all 50 corresponding historical result rows to
numerical precision. This does not establish a fresh seven-method refit or
source-to-result reproduction: the archived table used a separate legacy
solver, whose active file has changed since its frozen source hash. The
historical sample screen also used 2019 source-availability information, so it
must not be described as a strictly 2018-information-only prospective design.

## Version-0.1.3 source-replay boundary

The missing CEPII Gravity V202211 archive and RDS were restored from the
official CEPII download and matched the frozen size and SHA-256 references.
An isolated offline rebuild checked all 82 official-source files and
reconstructed the four selected processed inputs byte-for-byte, including
the response tensor, six-feature design, folds, and preprocessing registry.
The source checkout contains only the scripts, URLs, and hash references;
neither the restored data nor derived values enter the R package or public
source tree. The independent 50-fit authors-method replay used inputs with
the same hashes. These two checks validate the frozen-source-to-RGD/RGN
result chain in stages, but not current availability of every upstream API
response or a new fit of the five comparison methods. The source-checkout
workflow is not included in the built R package tarball.

## Version-0.1.4 public interface and simulation boundary

The reader-facing README and getting-started vignette use only the public
RGD/RGN, CPCA/HOSVD, decomposition, and regression APIs. The private real-data
reader accepts a data directory and keeps the frozen historical file label in
source-only mapping code. No URL has been invented for a GitHub repository
that does not yet exist, and the license remains a publication decision for
the copyright holders.

Source-only simulation runners now create new author-method fits under the
rebuttal decomposition and tensor-regression designs. Both full authors-only
grids have been run. All 240 tensor-regression generated-data fingerprints
match the historical archive, and selected cells agree with the current
active manuscript solver, but noisy method metrics do not pointwise match the
frozen regression table. Its historical solver has a different recorded hash.
The installed package still
contains the frozen replication-level records needed to reconstruct the
published comparisons, but not comparison optimizers or the fresh-run scripts.
