## git-workflow

This benchmark simulates a developer's Git workflow on the Chromium repository, focusing on metadata-heavy operations. It reproduces the process of applying a sequence of commits by extracting and applying patches in reverse chronological order, and re-committing them. This closely imitates the editing, staging, and committing cycle used during software development.

### Inputs

- `inputs/chromium/`: A working Git clone of the Chromium repository.
- `inputs/commits/`: A directory used to store commits and patches.

### Running

The benchmark performs the following steps:

1. **Reset Repository State**:
   - Checks out the `main` branch.
   - Deletes any previous benchmark branch (`bench_branch`) and creates a new one.
   - Performs a hard reset and cleans untracked files to ensure a clean starting state.

2. **Generate Patch Sequence**:
   - Extracts the latest N commits (default: 20) from the Chromium repository using `git rev-list`.
   - Stores the base commit (`HEAD~N`) and then:
     - Creates diff patches between each successive commit.
     - Captures their commit messages.
   - Patches and messages are stored in `inputs/commits/` as numbered `.diff` and `.commit` files.

3. **Replay Commits**:
   - Resets the working directory to the base commit.
   - Applies each patch in reverse order (oldest to newest).
   - After each patch:
     - Runs `git status`
     - Stages changes using `git add -A`
     - Commits with the original commit message and author identity.


### Validation

- This benchmark doesn't produce traditional outputs, but its correctness is implied by the successful application and recreation of all commits in order.

### References

- https://www.usenix.org/system/files/atc20-raghavan.pdf
- https://github.com/chromium/chromium