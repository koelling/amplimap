# Packaging
https://packaging.python.org/tutorials/packaging-projects/

## Prepare

    python3 -m pip install --upgrade build

# Release a new version

## Update `version.py`

## Tag

    git tag v0.4.15 && git push origin v0.4.15

## Build

    python3 -m build

## Upload (can use user/pass to login)

    python3 -m twine upload dist/amplimap-0.4.15.tar.gz