# GitHub Pages Deployment Setup Guide

## Quick Start (5 minutes)

### Step 1: Enable GitHub Pages

1. Go to your repository settings: `https://github.com/GKamundia/Phylos/settings`
2. Scroll down to "Pages" section
3. Under "Source", select **"GitHub Actions"**
4. Click "Save"

### Step 2: Test Local Build (Optional but Recommended)

```powershell
# On Windows PowerShell
.\scripts\deploy_local.ps1

# Or on WSL/Linux
bash scripts/deploy_local.sh
```

### Step 3: Trigger Deployment

The GitHub Action will automatically run when you:

- Push changes to the `main` branch that affect files in `auspice/` folder
- Or manually trigger it from the Actions tab

### Step 4: Access Your Dashboard

Once deployed, your dashboard will be available at:
**https://gkamundia.github.io/Phylos/**

## Manual Deployment (Alternative)

If you prefer manual control:

### 1. Build Static Site Locally

```bash
# Install Auspice if not already installed
npm install -g auspice

# Build static site
auspice build --extend ./auspice --outputDir ./docs

# Create CNAME file for custom domain (optional)
echo "your-domain.com" > ./docs/CNAME
```

### 2. Commit and Push

```bash
git add docs/
git commit -m "Deploy Auspice dashboard"
git push origin main
```

### 3. Configure GitHub Pages

- Repository Settings → Pages
- Source: "Deploy from a branch"
- Branch: `main`
- Folder: `/docs`

## Custom Domain Setup (Optional)

### 1. Add CNAME Record

In your DNS settings, add a CNAME record:

```
CNAME: dashboard.yourdomain.com → gkamundia.github.io
```

### 2. Update Repository Settings

- Go to Pages settings
- Add your custom domain: `dashboard.yourdomain.com`
- Enable "Enforce HTTPS"

## Automated Updates

The included GitHub Action (`.github/workflows/deploy.yml`) will:

1. **Trigger on changes** to Auspice JSON files
2. **Validate** JSON file integrity
3. **Build** static Auspice site
4. **Deploy** to GitHub Pages automatically

## Monitoring Deployments

- **Actions tab**: Monitor build status and logs
- **Environments**: View deployment history and status
- **Pages settings**: Check deployment URL and custom domain status

## Troubleshooting

### Build Fails

1. Check Actions tab for error logs
2. Ensure JSON files are valid
3. Verify Auspice can build locally

### Site Not Loading

1. Check if deployment completed successfully
2. Verify GitHub Pages is enabled
3. Check browser console for JavaScript errors

### Custom Domain Issues

1. Verify DNS propagation (can take up to 24 hours)
2. Check CNAME record is correct
3. Ensure HTTPS is enforced in settings

## Cost: **FREE**

- ✅ Unlimited bandwidth for public repos
- ✅ Custom domain support
- ✅ Automatic HTTPS
- ✅ Global CDN
- ✅ No server maintenance

Your Nextstrain dashboard will be live and accessible worldwide at no cost!
