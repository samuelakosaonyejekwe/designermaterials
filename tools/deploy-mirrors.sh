#!/bin/sh
# Push the current branch to every configured git remote so each static host serves the same build.
# Add mirrors with:  git remote add <name> <url>   (for example GitLab, Codeberg or Bitbucket Pages)
set -e
node tools/build-assets.mjs
branch=$(git rev-parse --abbrev-ref HEAD)
for remote in $(git remote); do
  echo "→ $remote"
  git push "$remote" "$branch"
done
