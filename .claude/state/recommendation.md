**Start #115 now (with #116 and #10 in the same PR).** #34 is decided (keep), so #115 is
the head of the only serial chain to 0.99.8: #115 → #117 → the cut, due ~9 Oct. That leaves
about a week for two PRs. #118 and #120 depend on nothing: run them in parallel sessions.

Two standing signals:

- **Large PRs get no automatic review right now.** `.claude/settings.json` defines
  `PostToolUse` twice, so `review-after-pr.py` is dropped (#148). Run `/pr-review <N>` by
  hand on Wave 2 PRs until #148 is fixed.
- **The parking lot must be empty before the 0.99.8 cut (#143).** Triage it in the next
  `/mesa-sweep`.
