
# How to Update Prof. Defeng Sun's Academic Homepage

This guide explains how to safely update and deploy the academic homepage of Prof. Defeng Sun.

---

## 📁Main Directory Structure Overview

```
defengwebsite.github.io/
├── _includes/
├── _layouts/
├── _pages/              ← Editable markdown content
├── assets/              ← CSS, fonts, images
├── images/              ← Profile and favicon images
├── index.md             ← The main homepage content (symlink or copy from _pages)
├── _config.yml          ← Site configuration
└── _data/navigation.yml ← Navigation bar items (order, names, URLs)
```

---

## ✏️ Steps to Update the Homepage

### 1. Edit Page Content

To update text, contact info, teaching, or publications, edit the markdown files in:

```
_pages/
```

Examples:
- `_pages/home.md` → main homepage section
- `_pages/publications.md` → publications list
- `_pages/recognitions.md` → honors and awards

You may use basic Markdown and HTML tags to format content.

---

### 2. Update Navigation Bar

To reorder or rename links at the top:

Edit:
```yaml
_data/navigation.yml
```

You can change the order like this:

```yaml
main:
  - title: "Publications"
    url: /publications/
  - title: "Professional Activities"
    url: /professional-activities/
```

---

### 3. Update Icons and Fonts (Optional)

Icons like ORCID and Google Scholar are supported via:

- **Academicons**: loaded in `custom.html`
- **Font Awesome**: for general icons

Edit:
```
_includes/head/custom.html
```

Make sure to preload fonts and include CSS like this:

```html
<link rel="preload" href="/assets/fonts/academicons.woff" as="font" type="font/woff" crossorigin="anonymous">
<link rel="stylesheet" href="/assets/css/academicons.css">
<link rel="stylesheet" href="https://cdnjs.cloudflare.com/ajax/libs/font-awesome/6.5.0/css/all.min.css">
```

---

### 4. Re-generate Homepage (index.md or index.html)

To set the main page:

- Option A: Copy the updated file as `index.md`
- Option B: Symlink `_pages/home.md` → `index.md`

Make sure `index.md` starts with YAML frontmatter:

```markdown
---
layout: home
title: "Welcome to Defeng Sun’s Academic Homepage"
---
```

---

### 5. Deploy to GitHub Pages

After editing, save and push:

```bash
git add .
git commit -m "update: homepage content and navigation"
git push
```

GitHub will auto-build and deploy the homepage at:

```
https://defengwebsite.github.io/
```

---

## 📤 Upload to Remote Server (if needed)

If you want to upload the site to a university server (via WinSCP or SSH):

1. Run Jekyll build (if applicable):
   ```bash
   bundle exec jekyll build
   ```
   This will output static files in `_site/`

2. Zip the `_site` directory:
   ```bash
   zip -r site.zip _site/
   ```

3. Upload `site.zip` to the server and unzip:

```bash
unzip site.zip
sudo rm -rf /var/www/html/*
sudo cp -r _site/* /var/www/html/
sudo chown -R www-data:www-data /var/www/html/
sudo chmod -R 755 /var/www/html/
```

---

## ✅ Notes

- Always preview locally with:
  ```bash
  bundle exec jekyll serve
  ```
- Clear browser cache (`Ctrl+Shift+R`) after each update.
- Commit each major edit with a meaningful message.

