# 10. Deployment Choice

## Decision For The MVP

Use **Streamlit**, run **locally**.

That keeps the iteration loop tight and avoids deployment / auth complexity
while the workflow is still being shaped by lab users.

The canonical plan describes the MVP as a local Streamlit app at
[`app/streamlit_app.py`](../app/streamlit_app.py), launched via
`streamlit run app/streamlit_app.py` and consuming the existing pipeline
APIs directly:

- `StockFindingPipeline.run` for Stage 1
- `StockSplittingPipeline.run` for Stage 2 (and optional Stage 3)
- `load_split_config` for validation of edited config

## Streamlit Trade-Offs (Honest Assessment)

Streamlit is great for:

- Rapid validation of the phenotype-first workflow.
- Fast iteration on the keyword-hit/no-hit result browser shape.
- Throwaway UI for in-lab user testing.

Streamlit is fragile for:

- Senior research scientists who do not want to deal with terminals,
  Python environments, dependency errors, or multiple browser tabs.
- Polished long-form workflows with multi-step state.
- Hosted multi-user access with auth.

These are explicitly accepted trade-offs for the first version.

## When To Graduate

Move beyond Streamlit when:

- The phenotype-first workflow and result-browser shape are validated by lab
  users.
- The RNAi planning module is mature.
- Senior users start asking for shared / hosted access.

## Successor Options

1. **Hosted FastAPI + React/Next.js**
   - Cleanest UX. Open URL → sign in → upload gene list → enter keywords →
     get phenotype hits and RNAi outputs.
   - Allowlist Allada/lab accounts at the auth layer.
   - Same Python pipeline runs server-side.
   - Best long-term default for non-coding scientists.
2. **Tauri or Electron desktop app**
   - One-click launcher; no terminal, no Python environment for end users.
   - Bundled Python backend.
   - Useful when offline use is required.
3. **NiceGUI (Python-native)**
   - Closer to Streamlit's ergonomics with more layout/flow control.
   - Does not solve deployment as cleanly as a hosted web app.

## Open Questions Deferred

- Hosted deployment, authentication, and multi-user job queues are explicitly
  out of scope for this first Streamlit version.
- Google Drive/Sheets integration is deferred until the upload-snapshot
  workflow is validated with lab users (see
  [`07-allada-snapshot.md`](07-allada-snapshot.md)).
- Polished React/FastAPI rebuild can come later if Streamlit proves the
  workflow.
