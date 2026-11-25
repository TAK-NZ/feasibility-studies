# feasibility-studies

Short, time-boxed studies that answer one question before TAK-NZ builds something: **can we do this, and is it worth the effort?**

Each study gets its own folder. It states the question and the pass/fail criteria up front (before looking at results), collects real data, runs the test, and ends with a written recommendation: go, no-go, or go with conditions. The code is kept so anyone can repeat the test or run it again on new data. Nothing here is production code.

## Studies

| Folder | Question | Verdict | Date |
|---|---|---|---|
| [drift-hindcast](drift-hindcast) | Is an OpenDrift ocean drift and search-area service for New Zealand SAR worth building on free data only (Copernicus Marine currents, GFS and ECMWF winds)? | **No-go as specified.** The model beat the hand-calculation baseline by 23.6% at 24 h (target: 30%). The 95% contour held 92% of observed endpoints, but only after a tuned layer outside OpenDrift. See the report. | 2026-11 |

## Adding a study

1. Create a folder named after the topic, in lowercase with hyphens.
2. Put these in it:
   - `README.md`: how to reproduce the study, step by step, with the exact commands.
   - `REPORT.md`: opens with the question, the pass/fail criteria that were agreed in advance, and the verdict. It then gives the numbers, the limitations, and what the test cannot show.
   - `LICENSE`: the licence that applies to that folder (see below).
   - `scripts/`, `results/` (small tables and figures that the report quotes), and `requirements.txt`.
3. Fix the pass/fail criteria before running the test on the data that decides them. Where a method was tuned, keep a tuning set apart from the test set and say so in the report.
4. Add a row to the table above.

## Rules for all studies

- **Licences are per folder.** The repository has no root licence on purpose. A study that uses a GPL-2.0-only library such as OpenDrift must not copy code from AGPL-3.0 TAK-NZ repositories, and the other way round. Each folder's `LICENSE` says what applies.
- **Downloaded data stays out of git.** `data/` is ignored. Each study's README must say how to download the data again and roughly how large it is.
- **No credentials in the repository.** Use the tool's own login command (for example `copernicusmarine login`, which stores the file in your home directory) or a git-ignored `.env`.
- **Report limitations honestly.** Say what the test cannot show as plainly as what it shows.
