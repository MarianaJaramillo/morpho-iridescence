# Contributing

Thank you for helping improve the Morpho iridescence datasets and analysis. This repository supports scientific results, so changes to source data, cleaning rules, parameters, and reported statistics must be traceable and reviewable.

## Before making a change

1. Read the root [README](README.md), [methods](docs/methods.md), [data dictionary](docs/data_dictionary.md), and the relevant dataset README.
2. Identify whether the change affects the Gruson pipeline, the Angular Atlas pipeline, or both. They are distinct workflows.
3. Do not overwrite or silently regenerate published outputs. Review the diff and retain a copy of existing results before running scripts that write into the data or figures directories.
4. Check data-owner permissions before adding, removing, or redistributing acquisition files.

## Data and naming conventions

- Use English for documentation, comments, generated column names, plot labels, and newly introduced internal identifiers.
- Preserve original instrument filenames, column headers, species-name spellings, and other source identifiers where they are required for parsing or provenance. Add or update an explicit mapping to a canonical English identifier rather than renaming source material in place.
- Do not silently change a scientific threshold, exclusion, correction, model, or measurement interpretation. Explain the reason, cite the decision or source, and update the methods and parameter documentation.
- Keep raw source files immutable. Record corrections and transformations in the existing configuration or cleaning logs.
- Never commit credentials, private participant information, or data whose redistribution rights have not been confirmed.

## Validation

For changes to the Gruson pipeline, run the regression tests from the repository root:

```bash
python -m pytest gruson/tests
```

For Atlas changes, run the affected processing step(s) only when their input files and output destinations are understood. Review generated tables, logs, and viewer output for unintended changes. Step 1 needs upstream workbooks that are not included in this repository; see [atlas/README.md](atlas/README.md#reproduce-the-processing-steps).

Include the commands and environment used in the contribution description. If validation could not be run, state why. A passing regression test does not by itself validate every Atlas output or establish scientific correctness.

## Change descriptions

A useful contribution description includes:

- the affected dataset, pipeline stage, and files;
- the motivation and any relevant literature or data-owner decision;
- whether source data, parameters, or scientific interpretations changed;
- validation commands and a concise summary of their results;
- any remaining limitations or required follow-up.
