# Changelog

All notable changes to this project should be documented in this file.

The format is inspired by Keep a Changelog and uses semantic version sections.

## [Unreleased]

### Added
- Bootstrap confidence intervals for calibrated probabilities in CLI outputs.
- Calibration stability utilities (LOO + stratified/repeated CV summaries).
- Reference performance reporting helpers with subgroup breakdowns.
- Structured metadata file for reference labels:
  `/home/runner/work/ile-predict/ile-predict/data/reference_labels.metadata.json`.

### Changed
- Lookup/CLI records now expose action-oriented status (`bundled-lookup`, `computed-from-structure`).
- Domain-edge outputs now surface explicit low-confidence policy in CLI.

### Testing
- Added CLI behavior tests and new calibration helper tests.
