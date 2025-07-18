# 📖 Complete GitHub Pages Setup Guide

## Step 1: Enable GitHub Pages in Repository Settings

### 🌐 Access Repository Settings

1. Go to your repository: `https://github.com/GKamundia/Phylos`
2. Click the **"Settings"** tab (top right of repository page)
3. Scroll down to **"Pages"** in the left sidebar

### ⚙️ Configure GitHub Pages

1. Under **"Source"**, select **"GitHub Actions"**
2. Click **"Save"**

![GitHub Pages Settings](https://docs.github.com/assets/cb-21851/images/help/pages/publishing-source-drop-down.png)

## Step 2: Commit and Push Your Changes

Since you already have the workflow file created, commit your changes:

```bash
# Add all files
git add .

# Commit with descriptive message
git commit -m "Add GitHub Pages deployment for Nextstrain dashboard"

# Push to main branch
git push origin main
```

## Step 3: Monitor the Deployment

### 📊 Watch GitHub Actions

1. Go to your repository
2. Click the **"Actions"** tab
3. You should see the workflow **"Deploy Simple Nextstrain Dashboard to GitHub Pages"** running

### 🔍 Check Deployment Status

The workflow will:

- ✅ Validate your JSON files
- ✅ Create simple deployment
- ✅ Deploy to GitHub Pages
- ✅ Provide you with a live URL

## Step 4: Access Your Live Dashboard

Once deployment completes (usually 2-3 minutes):

**Your dashboard will be live at:**

```
https://gkamundia.github.io/Phylos/
```

## Step 5: Verify Everything Works

### 🧪 Test Your Dashboard

1. Visit your GitHub Pages URL
2. Try downloading one of the JSON files
3. Go to [auspice.us](https://auspice.us)
4. Drag and drop the JSON file to verify it displays correctly

## Troubleshooting

### ❌ If Deployment Fails

**Check the Actions log:**

1. Go to Actions tab
2. Click on the failed workflow
3. Expand the failing step to see error details

**Common issues:**

- JSON files not found in `auspice/` folder
- Invalid JSON format
- Repository permissions

### 🔧 Manual Trigger

If you need to redeploy manually:

1. Go to Actions tab
2. Click "Deploy Simple Nextstrain Dashboard to GitHub Pages"
3. Click "Run workflow" button
4. Select "main" branch and click "Run workflow"

## What Happens Next

### 🔄 Automatic Updates

The workflow will automatically redeploy when you:

- Update files in the `auspice/` folder
- Push changes to the `main` branch

### 📁 File Structure on GitHub Pages

Your deployed site will contain:

```
https://gkamundia.github.io/Phylos/
├── index.html          # Landing page
├── rvf_L.json         # Large segment data
├── rvf_M.json         # Medium segment data
└── rvf_S.json         # Small segment data
```

### 👥 How Users Will Use Your Dashboard

1. **Visit your GitHub Pages URL**
2. **Download JSON files** from the dashboard
3. **Go to auspice.us**
4. **Drag and drop** JSON files for visualization

## Custom Domain (Optional)

### 🌐 Add Custom Domain

If you have a domain name:

1. **In repository settings:**

   - Go to Pages section
   - Add your domain in "Custom domain"
   - Save

2. **In your DNS settings:**

   - Add CNAME record: `your-domain.com` → `gkamundia.github.io`

3. **Enable HTTPS:**
   - Check "Enforce HTTPS" in Pages settings

## Next Steps

✅ **Immediate:** Commit and push to trigger first deployment
✅ **Test:** Verify dashboard works at your GitHub Pages URL  
✅ **Share:** Your dashboard will be publicly accessible
✅ **Update:** Any changes to `auspice/*.json` will auto-deploy

**Expected Timeline:**

- Setup: 5 minutes
- First deployment: 2-3 minutes
- Future updates: 1-2 minutes (automatic)

Your Nextstrain dashboard will be **live and free forever** on GitHub Pages! 🎉
