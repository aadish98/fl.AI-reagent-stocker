# 07. Allada Lab Sheet Snapshot

## MVP Decision: Upload, Not Google Auth

For the first version, the GUI does **not** integrate Google sign-in or
Google Drive/Sheets APIs. Instead, the user downloads the shared Allada Lab
Order Sheet as `.xlsx` and uploads that snapshot to the GUI.

Rationale:

- Avoids Google verification overhead and broad-scope permissions.
- Single-file uploads are explicit and auditable.
- Read-only by design: the GUI never writes back to the snapshot or the
  shared sheet.
- Easy to swap for a Google Drive Picker integration later, scoped to
  `drive.file` only.

The screen-planning flow stays **disabled** until a valid `.xlsx` is uploaded
and validated.

## Reference Fixture

[`Maddy_Playground/RNAi plan script TEST/Allada Lab Common Stocks 2026 (1).xlsx`](../Maddy_Playground/RNAi%20plan%20script%20TEST/Allada%20Lab%20Common%20Stocks%202026%20(1).xlsx).

That workbook today contains 27 worksheets with very heterogeneous layouts.
The parser must tolerate this fixture and any reasonable evolution of it.

## Documentation Without Leaking Lab Data

The expected snapshot format must be documented in a public docs page (e.g.,
`docs/allada_lab_sheet_snapshot.md`). The fixture itself stays out of git
(see [`09-data-privacy.md`](09-data-privacy.md)). Documentation describes
**supported layouts**, not private screen contents.

## Parser Requirements

The parser must be robust to:

- Multiple worksheets, some unrelated to stock inventory.
- New sheets being added by lab members.
- Repeated side-by-side box blocks within one sheet (e.g., "BL Box 1" through
  "BL Box 4" on the same tab, each with its own stock-ID column).
- Stock identifiers under varied headers, including:
  - `Stock Number`, `Stock #`, `Stk #`
  - `Bloomington Stock Number`, `Bloomington`
  - `Stock`, `Stock ID`, `Line ID`
  - `VDRC ID`, `VDRC`
- Stock IDs found below the first row (headers may sit on row 2 with row 1
  containing block titles).
- Numeric stock IDs (e.g., `5026.0` rounded to `5026`).
- Vienna IDs with optional leading `v` (`v60100` and `60100` should match).

What the parser must do per worksheet:

1. Detect zero or more likely stock-ID columns by scanning the first ~8 rows.
2. Index every recognizable stock ID under each detected column.
3. Preserve `(sheet name, block title, row number, column letter,
   raw value)` for every match so the GUI can show a useful "found in"
   location.
4. Skip only sheets that have no recognizable stock identifiers, with a
   visible warning rather than silent omission.
5. Surface ambiguous header detections to the user for review.

## Behavior Improvements Over Maddy's PowerShell Scripts

Today's PowerShell `Build-LabStockIndex` only picks **one** stock-ID column
per sheet, which misses many inventory boxes (e.g., `(GW) AD + DBD`,
`AD & DBD (GW) (archive)`, `CR (NU)`, `YK`). The GUI parser must:

- Scan every sheet.
- Detect headers beyond row 1.
- Handle multiple repeated blocks per sheet.
- Index multiple stock-ID columns per sheet.
- Never silently discard duplicate matches.
- Export ambiguous/skipped sheets for review.

## Status Card In The GUI

When a snapshot is uploaded, the GUI must display:

- Uploaded filename
- Detected sheets / tabs (with skip reasons where applicable)
- Number of stock IDs indexed total and per sheet
- Whether lab-stock matching is active
- Ambiguous header detections that may need review

## Output Effects

When a screen-planning run uses an uploaded snapshot:

- Candidate rows already present in the lab are flagged
  (e.g., `In lab stocks - sheet tab: <tab> (Allada Lab Common Stocks workbook)`).
- Row sorting can prioritize already-in-lab rows or, per user preference,
  exclude them from the order form.
- The status, tab name, and local row/block location are included as columns
  in `RNAi screening plan.xlsx` and `RNAi Order Form.xlsx`.
- The schema flowchart includes a "Lab snapshot match" checkpoint count.

## Hard Don'ts

- Don't write back to the uploaded workbook.
- Don't auto-merge multiple uploaded snapshots silently. If the user wants to
  combine snapshots, ask explicitly.
- Don't fail the run when one sheet is unrecognized; warn instead.
- Don't move the uploaded `.xlsx` outside the run folder unless the user
  explicitly saves it as reference data.
