# Free Deployment Options for RVF Nextstrain Dashboard

## Option 1: GitHub Pages with Static Auspice (Recommended)

### Setup Steps:

1. **Create a GitHub Pages deployment branch:**

   ```bash
   # Create a new orphan branch for GitHub Pages
   git checkout --orphan gh-pages
   git rm -rf .
   ```

2. **Set up Auspice static site:**

   ```bash
   # Install Auspice locally
   npm install -g auspice

   # Build static site from your JSON files
   auspice build --extend ./auspice --outputDir ./docs
   ```

3. **Configure GitHub Pages:**
   - Go to your repository settings
   - Enable GitHub Pages from the `gh-pages` branch
   - Your dashboard will be available at: `https://gkamundia.github.io/Phylos/`

### Automated Updates:

Create `.github/workflows/deploy.yml`:

```yaml
name: Deploy Auspice Dashboard

on:
  push:
    branches: [main]
    paths: ["auspice/**"]

jobs:
  deploy:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v3

      - name: Setup Node.js
        uses: actions/setup-node@v3
        with:
          node-version: "18"

      - name: Install Auspice
        run: npm install -g auspice

      - name: Build static site
        run: auspice build --extend ./auspice --outputDir ./build

      - name: Deploy to GitHub Pages
        uses: peaceiris/actions-gh-pages@v3
        with:
          github_token: ${{ secrets.GITHUB_TOKEN }}
          publish_dir: ./build
```

## Option 2: Netlify (Easy Drag & Drop)

1. **Build static site locally:**

   ```bash
   npm install -g auspice
   auspice build --extend ./auspice --outputDir ./netlify-build
   ```

2. **Deploy to Netlify:**

   - Go to https://netlify.com
   - Drag and drop the `netlify-build` folder
   - Get instant URL like: `https://wonderful-name-123456.netlify.app`

3. **Automated updates:** Connect your GitHub repo to Netlify for automatic rebuilds.

## Option 3: Vercel

Similar to Netlify but with different interface:

1. Install Vercel CLI: `npm install -g vercel`
2. Build static site: `auspice build --extend ./auspice --outputDir ./vercel-build`
3. Deploy: `vercel --prod ./vercel-build`

## Option 4: Railway.app (Container Deployment)

For your Docker-based setup:

1. **Create `railway.dockerfile`:**

   ```dockerfile
   FROM nginx:alpine
   COPY auspice/ /usr/share/nginx/html/
   COPY deployment/nginx.conf /etc/nginx/conf.d/default.conf
   EXPOSE 80
   ```

2. **Deploy to Railway:**
   - Connect GitHub repo to Railway
   - Deploy from Dockerfile
   - Get URL like: `https://your-app.railway.app`

## Option 5: Render.com (Free Static Sites)

1. **Build static site**
2. **Connect GitHub repo to Render**
3. **Set build command:** `npm install -g auspice && auspice build --extend ./auspice --outputDir ./build`
4. **Set publish directory:** `build`

## Recommended Approach: GitHub Pages

For your use case, **GitHub Pages is the best option** because:

- ✅ **Free forever** for public repositories
- ✅ **Custom domain support** (optional)
- ✅ **Automatic HTTPS**
- ✅ **Easy CI/CD integration**
- ✅ **Perfect for static Auspice dashboards**
- ✅ **Good performance** for genomic data visualization

## Cost Comparison:

| Platform         | Cost              | Custom Domain | Auto-Deploy | Docker Support |
| ---------------- | ----------------- | ------------- | ----------- | -------------- |
| **GitHub Pages** | Free              | Yes           | Yes         | No             |
| **Netlify**      | Free (100GB/mo)   | Yes           | Yes         | No             |
| **Vercel**       | Free (100GB/mo)   | Yes           | Yes         | Limited        |
| **Railway**      | $5/mo after trial | Yes           | Yes         | Yes            |
| **Render**       | Free (100GB/mo)   | Yes           | Yes         | Yes            |

## Next Steps:

1. Choose your preferred platform
2. Build static Auspice site from your JSON files
3. Set up automated deployment
4. Configure custom domain (optional)
5. Test dashboard functionality

Would you like me to help you set up any of these deployment options?
