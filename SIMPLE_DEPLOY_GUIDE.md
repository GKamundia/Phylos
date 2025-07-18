# Simple Auspice Dashboard Deployment Guide

## Option 1: GitHub Pages (Static Files) - SIMPLEST

### Step 1: Create a simple index.html file

Since your Auspice JSON files are ready, we'll create a basic static site.

### Step 2: Enable GitHub Pages

1. Go to: https://github.com/GKamundia/Phylos/settings/pages
2. Under "Source", select "Deploy from a branch"
3. Choose "main" branch and "/ (root)" folder
4. Click "Save"

### Step 3: Your dashboard will be live at:

**https://gkamundia.github.io/Phylos/**

## Option 2: Netlify Drop (Instant Deploy)

### Step 1: Build static site locally

1. Install Auspice: `npm install -g auspice`
2. Build: `auspice build --extend ./auspice --outputDir ./netlify-deploy`

### Step 2: Deploy to Netlify

1. Go to: https://app.netlify.com/drop
2. Drag and drop the `netlify-deploy` folder
3. Get instant URL like: `https://wonderful-name-123456.netlify.app`

## Option 3: Use Nextstrain View (Immediate)

### Test your dashboard right now:

1. Go to: https://nextstrain.org/community/
2. Upload your JSON files directly
3. View immediately without any setup

This is the fastest way to see your dashboard working!

## Files needed for GitHub Pages:

- Your existing auspice/\*.json files ✅
- A simple index.html file (I'll create this)
- That's it!

## Next Steps:

Choose Option 1 (GitHub Pages) for a permanent free solution, or Option 3 (Nextstrain Community) to see it working immediately.
