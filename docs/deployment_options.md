# Quick GitHub Pages Deployment Guide

## Three Deployment Options

### Option 1: Simple JSON Hosting (Easiest - No installations needed)

**What it does:** Hosts your JSON files and uses nextstrain.org to display them

**Steps:**

1. Run the simple deployment script:

   ```powershell
   .\scripts\deploy_simple.ps1
   ```

2. Copy contents of `deploy\` folder to GitHub Pages:

   ```bash
   # Copy deploy folder contents to docs folder for GitHub Pages
   cp -r deploy/* docs/
   ```

3. Enable GitHub Pages in repository settings (source: docs folder)

4. Update the URLs in `docs/index.html` to match your GitHub Pages URL

**Pros:**

- ✅ No Node.js or Auspice installation required
- ✅ Lightweight deployment
- ✅ Uses official Nextstrain viewer

**Cons:**

- ❌ Requires internet connection to view
- ❌ Depends on nextstrain.org service

### Option 2: Full Auspice Build (Professional)

**What it does:** Creates a complete self-contained dashboard

**Requirements:** Node.js installed

**Steps:**

1. Install Node.js from https://nodejs.org/
2. Run the full build script:
   ```powershell
   .\scripts\deploy_full.ps1
   ```
3. Deploy the `build/` folder to any web hosting

**Pros:**

- ✅ Self-contained (works offline)
- ✅ Professional appearance
- ✅ Full Auspice functionality

**Cons:**

- ❌ Requires Node.js installation
- ❌ Larger file size

### Option 3: Manual Upload to Auspice.us (Simplest)

**What it does:** Users manually upload your JSON files to auspice.us

**Steps:**

1. Host your JSON files anywhere (GitHub Pages, Dropbox, etc.)
2. Users download JSON files and drag them to https://auspice.us

**Pros:**

- ✅ No installation required
- ✅ Always uses latest Auspice version
- ✅ No hosting complexity

**Cons:**

- ❌ Manual process for users
- ❌ Requires internet connection

## Recommended for You: Option 1 (Simple JSON Hosting)

Since you want to deploy without installing packages, run:

```powershell
.\scripts\deploy_simple.ps1
```

This will create a `deploy/` folder ready for GitHub Pages hosting!
