# Working Instructions

## Issue Workflow

Follow these steps exactly when working on a GitHub issue.

### 1. Starting an issue

Before writing any code:

```bash
# Mark the issue in progress
gh issue edit <N> --add-label "in-progress"
gh issue comment <N> --body "Starting work on this."
```

Update `TODO.md`: change the relevant `- [ ]` items to `- [~]` (in progress).

### 2. While working

After each meaningful change (a fix, a new function, a passing test):

```bash
gh issue comment <N> --body "<one-line description of what was just done>"
```

Update the relevant `- [ ]` sub-items in `TODO.md` to `- [x]` as they are completed.

### 3. Completing an issue

When all sub-tasks are done and tests pass:

```bash
# Close the issue with a summary comment
gh issue close <N> --comment "Done. <brief summary of what was implemented/fixed>"

# Remove in-progress label
gh issue edit <N> --remove-label "in-progress"
```

Update `TODO.md`:
- Mark all sub-items `- [x]`
- Mark the parent item `- [x]`

Commit `TODO.md` together with the code changes:

```bash
git add TODO.md <changed source files>
git commit -m "<summary of work>

Closes #<N>

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"
```

### 4. Blocked or scope change

If the issue turns out to be blocked or larger than expected:

```bash
gh issue comment <N> --body "Blocked: <reason>. <next step or dependency>"
gh issue edit <N> --add-label "blocked"
```

Update `TODO.md` to note the blocker inline.

---

## TODO.md Conventions

| Marker | Meaning |
|--------|---------|
| `- [ ]` | Not started |
| `- [~]` | In progress |
| `- [x]` | Done |

Always keep `TODO.md` in sync with GitHub issues. The file is the local source of truth;
GitHub issues are the shared/public view.

---

## Labels

| Label | Meaning |
|-------|---------|
| `phase-1` … `phase-7` | Which implementation phase |
| `in-progress` | Currently being worked on |
| `blocked` | Waiting on a dependency |
