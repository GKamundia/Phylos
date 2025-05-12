# filepath: c:\Users\Anarchy\Documents\Data_Science\NextStrain\rvf-nextstrain\docs\embedding_capabilities.md

# RVF-Nextstrain Dashboard Embedding Documentation

This document outlines how to embed the RVF-Nextstrain dashboard or its components into external websites.

## Embedding Options

The RVF-Nextstrain dashboard supports three embedding options:

1. **Full Dashboard Embedding** - Embed the complete Auspice interface
2. **Component Embedding** - Embed specific visualization components
3. **Data Integration** - Use the JSON API to integrate with custom visualizations

## 1. Full Dashboard Embedding

### Using an iframe

The simplest way to embed the complete dashboard is using an iframe:

```html
<iframe
  src="https://rvf-nextstrain.example.org/?c=country&r=region"
  width="100%"
  height="700px"
  style="border: none; overflow: hidden;"
  title="RVF-Nextstrain Dashboard"
  loading="lazy"
  allowfullscreen
></iframe>
```

URL Parameters for Customization
The embedded dashboard can be customized using URL parameters:

d - Dataset to load, e.g., ?d=rvf_L for L segment data
c - Color by attribute, e.g., ?c=country
r - Filter by region, e.g., ?r=africa
m - Map view, e.g., ?m=country
p - Panel layout, e.g., ?p=tree-map
s - Sidebar, e.g., ?s=closed to hide the sidebar
Example with combined parameters:
https://rvf-nextstrain.example.org/?d=rvf_L&c=host&r=africa&p=tree-map

## 2. Embedding Script

Add this JavaScript to your website for responsive iframe sizing:

```html
<script>
  function resizeAuspiceIframe(id) {
    const iframe = document.getElementById(id);
    const resize = () => {
      const width = iframe.parentElement.clientWidth;
      const height = window.innerHeight * 0.8;
      iframe.style.width = `${width}px`;
      iframe.style.height = `${height}px`;
    };
    window.addEventListener("resize", resize);
    resize();
  }

  document.addEventListener("DOMContentLoaded", () => {
    resizeAuspiceIframe("nextstrain-iframe");
  });
</script>
```

## 2. Component Embedding

For more integration flexibility, we provide embeddable components that can be placed anywhere in your web page.

### Tree Component

```html
<div id="rvf-tree" data-dataset="rvf" data-component="tree"></div>

<script src="https://rvf-nextstrain.example.org/js/embed.js"></script>
<script>
  RVFNexstrain.embedComponent("rvf-tree", {
    height: 500,
    colorBy: "host",
    filters: {
      region: "africa",
      segment: "L",
    },
  });
</script>
```

### Map Component

```html
<div id="rvf-map" data-dataset="rvf" data-component="map"></div>

<script src="https://rvf-nextstrain.example.org/js/embed.js"></script>
<script>
  RVFNexstrain.embedComponent("rvf-map", {
    height: 400,
    resolution: "country",
    showTransmissions: true,
  });
</script>
```

### Component Options

Common options for all components:

| Option  | Type          | Description                        |
| ------- | ------------- | ---------------------------------- |
| height  | Number        | Height in pixels                   |
| width   | Number/String | Width in pixels or percentage      |
| colorBy | String        | Attribute to color by              |
| filters | Object        | Key-value pairs for filtering data |

Component-specific options:

### Tree Component

- `layout`: Tree layout (`rectangular, radial, unrooted, clock`)
- `showBranches`: Show branch lengths (boolean)

### Map Component

- `resolution`: Geographic resolution (region, country, division)
- `showTransmissions`: Show transmission lines (boolean)

## 3. Data Integration via API

For advanced integration, use our API to fetch data for custom visualizations.

### API Endpoints

- `GET /api/v1/dataset/{dataset_name}` - Get full Auspice JSON data
- `GET /api/v1/dataset/{dataset_name}?format=tree_only` - Get only tree data
- `GET /api/v1/dataset/{dataset_name}?format=metadata_only` - Get only metadata

Example: Custom D3.js Visualization

```bash
<div id="custom-rvf-viz"></div>

<script src="https://d3js.org/d3.v7.min.js"></script>
<script>
async function createCustomViz() {
  // Fetch tree data
  const response = await fetch(
    'https://rvf-nextstrain.example.org/api/v1/dataset/rvf_L?format=tree_only'
  );
  const data = await response.json();

  // Create custom visualization with D3
  const svg = d3.select('#custom-rvf-viz')
    .append('svg')
    .attr('width', 800)
    .attr('height', 600);

  // Custom visualization code...
}

document.addEventListener('DOMContentLoaded', createCustomViz);
</script>
```

### Implementation Example

Create an embeddable JavaScript library to enable these embedding capabilities:

```bash
/**
 * RVF-Nextstrain Embedding Library
 * Allows embedding Nextstrain visualizations in external websites
 */

(function(window) {
  'use strict';

  // Base URL for the Nextstrain instance
  const BASE_URL = 'https://rvf-nextstrain.example.org';
  const API_BASE = `${BASE_URL}/api/v1`;

  // Default configuration
  const DEFAULT_CONFIG = {
    height: 400,
    width: '100%',
    colorBy: 'country',
    layout: 'rectangular',
    showBranches: true,
    resolution: 'country',
    showTransmissions: true,
    filters: {}
  };

  /**
   * Converts options object to URL parameters
   */
  function optionsToParams(options) {
    const params = new URLSearchParams();

    if (options.colorBy) params.set('c', options.colorBy);

    if (options.filters) {
      for (const [key, value] of Object.entries(options.filters)) {
        if (key === 'region') params.set('r', value);
        else if (key === 'segment') params.set('s', value);
        // Add more filter mappings as needed
      }
    }

    return params.toString();
  }

  /**
   * Embeds a full Auspice dashboard in an iframe
   */
  function embedDashboard(elementId, options = {}) {
    const element = document.getElementById(elementId);
    if (!element) {
      console.error(`Element with ID "${elementId}" not found`);
      return;
    }

    const config = { ...DEFAULT_CONFIG, ...options };
    const dataset = options.dataset || element.dataset.dataset || 'rvf';
    const params = optionsToParams(config);

    const iframe = document.createElement('iframe');
    iframe.src = `${BASE_URL}/?d=${dataset}&${params}`;
    iframe.width = config.width;
    iframe.height = config.height;
    iframe.style.border = 'none';
    iframe.style.overflow = 'hidden';
    iframe.title = `${dataset.toUpperCase()} Nextstrain Dashboard`;
    iframe.setAttribute('loading', 'lazy');
    iframe.setAttribute('allowfullscreen', '');

    element.appendChild(iframe);

    // Add responsive behavior if requested
    if (config.responsive) {
      const resizeObserver = new ResizeObserver(entries => {
        for (const entry of entries) {
          iframe.style.width = `${entry.contentRect.width}px`;
        }
      });

      resizeObserver.observe(element);
    }

    return iframe;
  }

  /**
   * Loads the Auspice component script if not already loaded
   */
  function loadComponentScript() {
    return new Promise((resolve, reject) => {
      if (window.auspiceComponents) {
        resolve(window.auspiceComponents);
        return;
      }

      const script = document.createElement('script');
      script.src = `${BASE_URL}/js/auspice-components.js`;
      script.async = true;
      script.onload = () => resolve(window.auspiceComponents);
      script.onerror = () => reject(new Error("Could not load Auspice components"));

      document.head.appendChild(script);
    });
  }

  /**
   * Embeds a specific visualization component
   */
  async function embedComponent(elementId, options = {}) {
    const element = document.getElementById(elementId);
    if (!element) {
      console.error(`Element with ID "${elementId}" not found`);
      return;
    }

    const config = { ...DEFAULT_CONFIG, ...options };
    const dataset = options.dataset || element.dataset.dataset || 'rvf';
    const component = options.component || element.dataset.component || 'tree';

    try {
      // Load auspice-components.js if not already loaded
      const auspiceComponents = await loadComponentScript();

      // Set element dimensions
      element.style.width = config.width;
      element.style.height = `${config.height}px`;

      // Fetch data for the component
      const format = component === 'tree' ? 'tree_only' :
                    component === 'map' ? 'metadata_only' : 'full';

      const response = await fetch(`${API_BASE}/dataset/${dataset}?format=${format}`);
      if (!response.ok) {
        throw new Error(`Failed to load data: ${response.status} ${response.statusText}`);
      }

      const data = await response.json();

      // Initialize the component
      switch (component) {
        case 'tree':
          new auspiceComponents.Tree({
            target: element,
            props: {
              data: data.tree,
              layout: config.layout,
              colorBy: config.colorBy,
              showBranches: config.showBranches,
              filters: config.filters
            }
          });
          break;

        case 'map':
          new auspiceComponents.Map({
            target: element,
            props: {
              data: data,
              geoResolution: config.resolution,
              colorBy: config.colorBy,
              showTransmissions: config.showTransmissions,
              filters: config.filters
            }
          });
          break;

        case 'entropy':
          new auspiceComponents.Entropy({
            target: element,
            props: {
              data: data,
              geneMap: data.genome_annotations || {},
              fill: config.fill || "#1f77b4"
            }
          });
          break;

        default:
          console.error(`Unknown component type: ${component}`);
      }

      return element;
    } catch (error) {
      console.error(`Error embedding component: ${error}`);
      element.innerHTML = `<div class="error">Failed to load ${component}: ${error.message}</div>`;
    }
  }

  /**
   * Public API
   */
  window.RVFNextstrain = {
    embedDashboard,
    embedComponent,

    // Utility to load dataset info (versions, etc.)
    loadDatasetInfo: async (datasetName) => {
      const response = await fetch(`${API_BASE}/versions/${datasetName}`);
      return response.json();
    }
  };

})(window);
```
