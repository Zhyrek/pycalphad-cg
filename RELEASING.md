# Releasing to PyPI

Releases are built and uploaded by GitHub Actions (`.github/workflows/release.yml`): wheels for Linux, macOS and
Windows (Python 3.11 to 3.14) and a source distribution are built on every push. They are uploaded to PyPI only
when a GitHub Release is published. The upload uses PyPI's *trusted publishing*, so no API token is stored anywhere.

## One-time setup

1. **PyPI account.** Create one at <https://pypi.org/account/register/> and turn on two-factor authentication
   (PyPI requires it for uploads).
2. **Pending publisher.** Register the GitHub workflow as allowed to create the project, before its first release:
   <https://pypi.org/manage/account/publishing/>, "Add a new pending publisher", GitHub tab:
   - PyPI project name: `pycalphad-cg`
   - Owner: `Zhyrek`
   - Repository name: `pycalphad-cg`
   - Workflow name: `release.yml`
   - Environment name: `pypi`
3. **GitHub environment.** In the repository: Settings, Environments, "New environment", name `pypi`. Optionally
   add yourself under "Required reviewers": every upload then waits for your approval in the Actions tab.

## Each release

1. Make sure the Tests workflow is green on `main` (both pycalphad 0.11.2 and develop).
2. Set the version in `pycalphad_cg/_version.py` (e.g. `0.1.0`), commit and push.
3. Create the release on GitHub: Releases, "Draft a new release", tag `v0.1.0` (created on publish), title
   `0.1.0`, notes. Publishing it starts the workflow.
4. Watch the "Build and release" run in the Actions tab (approve the `pypi` environment if you required a
   reviewer). When it is done, check <https://pypi.org/project/pycalphad-cg/> and, in a fresh environment,
   `pip install pycalphad-cg`.

A version number can be uploaded only once. A broken release is fixed with a new version, after
"yanking" the bad one on PyPI if needed.

## Optional: a dry run on TestPyPI

TestPyPI (<https://test.pypi.org>) is a separate site with its own account and pending publisher (same fields as
above, environment e.g. `testpypi`). To try the upload there first, add a job like `upload_pypi` with
`repository-url: https://test.pypi.org/legacy/` under `pypa/gh-action-pypi-publish`, and
`environment: testpypi`. Then install with
`pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ pycalphad-cg`.
