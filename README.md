# sinanyavuz.github.io

Personal site: projects, browser tools, and (soon) a blog. Built with
[Jekyll](https://jekyllrb.com/), which GitHub Pages builds automatically on every push.
There's no build step and no dependencies to manage.

## Common edits

| I want to…                 | Edit                                                    |
|----------------------------|---------------------------------------------------------|
| Change name, bio links, CV | `_config.yml` (`author:` block)                         |
| Change the intro / About   | `index.html`                                            |
| Add a project              | Add an entry to `_data/projects.yml`                    |
| Add a browser tool         | Create `tools/<name>/index.html` + entry with `type: tool` |
| Write a blog post          | Add `_posts/YYYY-MM-DD-title.md` (see `_drafts/` for an example) |
| Add a CV                   | Put the PDF in `assets/`, set `author.resume` in `_config.yml` |

The **Blog** nav link appears automatically once the first post exists.

### Blog post template

```markdown
---
title: My post title
description: One-line summary shown in the list and in link previews.
tags: [data-science]
---

Write in Markdown here.
```

## Preview locally (optional)

```bash
bundle install
bundle exec jekyll serve --livereload   # http://localhost:4000
```

## Publishing

1. In the repo, open **Settings → Pages**.
2. Under *Build and deployment*, choose **Deploy from a branch**, branch `main`, folder `/ (root)`.
3. The site is live a minute or two after each push to `main`.
