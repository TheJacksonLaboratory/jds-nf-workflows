# Release Versioning

This repository uses [Semantic Versioning (SemVer)](https://semver.org/) for releases after `1.0.0`.

A release version has the form `MAJOR.MINOR.PATCH`, for example `1.2.3`:

- **MAJOR** changes introduce incompatibilities with the documented workflow interface or expected outputs.
- **MINOR** changes add functionality while preserving backward compatibility.
- **PATCH** changes correct backward-compatible defects or make internal improvements.

The release version applies to the repository as a whole. Published releases should use matching Git tags such as `v1.0.0`, `v1.1.0`, and `v1.1.1`.

## Compatibility Contract

For this repository, the public API includes:

- Workflow names and launch commands, such as `--workflow rnaseq`.
- Workflow parameters, including parameter names, required values, defaults, and accepted input formats.
- Input and output file formats, output filenames, directory organization, and sample identifiers.
- Module and subworkflow interfaces used by repository consumers or downstream workflows.
- Supported Nextflow versions and documented runtime requirements.
- The intended interpretation of analysis results and published reference or database resources.

The size of a code change does not determine its SemVer level. A large internal refactor can be a minor or patch release when it is transparent to users. A one-line parameter rename is a major release if it breaks existing launch commands.

## Version Increments

| Change | Version increment | Repository examples |
| --- | --- | --- |
| Breaking change | **MAJOR** | Remove or rename a parameter; change a required input; alter an output filename or schema; change a channel or module contract; drop a supported Nextflow version; change a default in a way that invalidates existing runs |
| Backward-compatible feature | **MINOR** | Add a workflow; add an optional parameter; add an optional output; support a new input format or genome build; add a new caller while preserving existing behavior |
| Backward-compatible fix | **PATCH** | Correct an error without changing the documented interface; fix metadata or logging; improve error handling; adjust resources; update a container while preserving expected tool behavior |
| Documentation or test-only change | Usually no release, or **PATCH** if published | Update the README or wiki, improve test coverage, or correct release-note wording |

After incrementing one component, reset the components to its right:

- `1.2.3` -> `1.2.4` for a patch release.
- `1.2.3` -> `1.3.0` for a minor release.
- `1.2.3` -> `2.0.0` for a major release.

## Repository-Specific Guidance

### Major releases

Use a major release when users must change commands, inputs, downstream processing, or result interpretation. Examples include:

- Renaming or removing a workflow parameter.
- Making an optional input required, or removing support for an existing input format.
- Changing output names, formats, or directory structure in a way that breaks downstream scripts.
- Changing the sample identifier or metadata contract.
- Changing module or subworkflow inputs and outputs used by external consumers.
- Dropping a supported Nextflow version or changing a runtime requirement that existing users cannot satisfy without an environment change.
- Changing defaults so that existing launch commands produce materially incompatible analyses and migration is required.

The `1.0.0` release is a logical major-version boundary because the Nextflow DSL2 and strict-syntax refactor changed the implementation across workflows, modules, and scripts. A future internal refactor should not automatically require another major release if the documented interface remains compatible.

### Minor releases

Use a minor release for additive capabilities that existing commands continue to support. Examples include:

- Adding a new workflow such as a new simulation or conversion workflow.
- Adding an optional parameter, optional output, or optional analysis mode.
- Supporting an additional genome build, sequencing technology, or input format without removing existing support.
- Adding a new caller or report while preserving existing outputs and defaults.
- Adding a new configuration profile that does not change existing profiles.

For example, adding an optional HLA-typing mode while leaving existing workflow defaults and outputs intact would be a minor release.

### Patch releases

Use a patch release for fixes and implementation changes that preserve the documented interface. Examples include:

- Fixing a workflow failure or incorrect command while retaining the same parameters and output contract.
- Correcting metadata, logging, or report formatting.
- Increasing memory or wall-clock allocations to allow existing analyses to complete.
- Improving error handling or validation.
- Updating a container image when the replacement is expected to provide equivalent behavior and output compatibility.

A patch release can still require a prominent release-note warning. For example, fixing a variant-calling bug may change results for affected inputs even though users do not need to change their launch command.

## Special Cases

### Container updates

Use a patch release when a container replacement is an implementation update with equivalent expected behavior. Use a minor release when it adds a new capability or changes supported behavior. Use a major release when it changes output compatibility or requires users to change how they invoke or consume a workflow.

### Reference and database updates

Reference and database changes can alter scientific results without changing the command-line interface. Use a minor release for a deliberate refresh of a reference or database resource, and identify the affected workflows and resource versions in the release notes. Use a major release when the update also changes output schemas, invalidates downstream integrations, or requires migration.

### Bug fixes that change results

A fix that restores the intended documented behavior is normally a patch release, even if corrected results differ from results produced by an earlier version. Release notes should identify the affected workflows, describe the previous behavior, and state whether prior analyses should be reviewed or rerun.

### Resource changes

Increasing memory, CPU, or wall-clock limits is normally a patch release because it does not change the workflow interface. If resource changes alter scheduling assumptions, required infrastructure, or supported deployment environments, consider a minor or major release based on the compatibility impact.

### Pre-releases

Use SemVer pre-release identifiers for release candidates or other versions that are not yet considered stable, for example:

- `1.1.0-alpha.1`
- `1.1.0-rc.1`

Pre-releases sort before the corresponding final release and should not be treated as production-stable versions.

## Release Process

Before publishing a release:

1. Review changes to workflows, parameters, inputs, outputs, defaults, containers, references, and supported Nextflow versions.
2. Classify the highest-impact change as major, minor, or patch.
3. Run the relevant `nf-test` tests and `nextflow lint`.
4. Update `ReleaseNotes.md` with affected workflows, user-visible changes, and migration or reproducibility notes.
5. Create an annotated Git tag whose name matches the release version, for example:

   ```bash
   git tag -a v1.0.1 -m "Release 1.0.1"
   git push origin v1.0.1
   ```

6. Publish the corresponding GitHub and Zenodo release using the same version.

When classification is uncertain, choose the higher compatibility level and document the reasoning. It is preferable to communicate a potentially breaking change clearly than to make users discover it through a failed run or changed downstream results.
