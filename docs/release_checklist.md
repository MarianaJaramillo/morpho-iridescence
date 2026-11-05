# External release checklist

Use this checklist before sharing this repository with an external research team or making it public. Several items require decisions or evidence from the project owners; they cannot be inferred safely from the source code or publications.

## Ownership, permissions, and citation — required

- [x] Identify the responsible authors, data custodians, and contact email for questions and corrections. (Mariana Jaramillo Acero, Universidad de Antioquia.)
- [ ] Confirm who owns the code and who owns each dataset and acquisition file (institution, co-authors, collaborators).
- [ ] Obtain explicit permission to share every included workbook, spectrum, table, figure, notebook output, and generated viewer.
- [x] Select and add a code license. (MIT, see `LICENSE`.)
- [x] Select and add a data license. (CC BY 4.0, see `DATA_LICENSE.md`.)
- [x] Add a project citation record (`CITATION.cff`). Review the author list and add co-authors if any.
- [ ] Add an archival DOI (for example through Zenodo) once the first release exists, and add it to `CITATION.cff`.
- [ ] Confirm that third-party software, fonts, images, and source materials are used and redistributed under compatible terms.

## Data provenance and curation — required

- [ ] Record the source, version/date, custodian, and checksum for every original acquisition archive.
- [ ] For Atlas step 1, provide stable source locations or access instructions for `Spectra.zip`, the remeasurement files, and the two correction workbooks listed in [atlas/README.md](../atlas/README.md#reproduce-the-processing-steps), subject to permission.
- [ ] Document any replacements, manual corrections, exclusions, and duplicate-resolution decisions, including who approved them and when.
- [ ] Check that no workbook metadata, comments, hidden sheets, or embedded objects contain sensitive or non-shareable information.
- [ ] Confirm that public names, taxonomy, specimen identifiers, and provenance statements are accurate and appropriate.

## Reproducibility and validation — required

- [ ] Record the operating system, Python version, R version, and exact resolved package versions used to generate the release outputs. (Python packages: see `requirements-lock.txt`; R is not recorded.)
- [x] Add a reproducible Python environment specification (`requirements.txt`, `requirements-lock.txt`). No lock file exists for R.
- [x] Run the Gruson regression suite and the Atlas integrity tests (14 passed, Python 3.11.5).
- [x] Make the Gruson pipeline output deterministic: row/column order no longer depends on `PYTHONHASHSEED`; `gruson/data/processed/` was regenerated and a test guards it.
- [ ] Validate the Atlas outputs against the included inputs; clearly distinguish the reproducible steps 2–7 from step 1, whose original workbooks are absent here.
- [ ] Regenerate the offline viewer from a clean output directory and confirm it opens without a network connection.
- [ ] Confirm that notebooks execute from a clean kernel in the documented order, or document and justify any notebook that requires unavailable files or manual state.
- [ ] Check generated outputs, figure provenance, CSV encodings, and data counts against the documented summaries.
- [x] Add automated tests for Atlas input/output integrity (`atlas/tests/`) and a GitHub Actions workflow that runs all tests.

## GitHub repository hygiene — required

- [ ] Verify the actual Git working tree and inspect all staged, unstaged, and ignored files before the first push.
- [x] Complete `.gitignore` (Python, R, editor and OS files).
- [x] Decide where the 31 MB offline viewer lives: inside Git (`atlas/viewer/`).
- [x] Notebook outputs are kept so readers can see results without running them.
- [ ] Add issue and contribution guidance, a code of conduct if appropriate for the collaboration, and a security/contact policy if the repository is public.
- [ ] Ensure all repository documentation and newly introduced labels are in English. Preserve source-language filenames and headers only where needed for fidelity, provenance, or compatibility, and explain the mapping to English internal names.
- [ ] Create a versioned release only after permissions, licenses, provenance, and validation are complete.

## Current known limitations

- The original Atlas acquisition workbooks and stable upstream source locations/checksums are not included in this distribution (decision pending).
- No lock file exists for R packages.
- Notebooks and Atlas steps 1 and 3–7 have not been re-run in a clean environment as part of this checklist.
