# Metadata-Driven Auspice Configuration

This document explains how the RVF-Nextstrain project uses metadata to dynamically generate and select Auspice configurations based on different pathogens and dataset characteristics.

## Architecture Overview

The system uses a three-tiered approach to Auspice configuration:

1. **Base Configuration**: Static JSON files stored in `pathogens/<pathogen>/config/auspice_config.json`
2. **Dynamic Enhancement**: Runtime generation of enhanced configurations based on actual metadata
3. **Multi-Pathogen Support**: Selection of appropriate configurations via `master_config.yaml`

## How It Works

### Configuration Selection

The workflow automatically selects the appropriate base configuration using this process:

1. The `active_pathogen` value in `config/master_config.yaml` determines which pathogen is being analyzed
2. The system looks up pathogen-specific settings in the `pathogens` section:
   ```yaml
   pathogens:
     rvf:
       name: "Rift Valley Fever virus"
       config_path: "pathogens/rvf/config/config.yaml"
       auspice_config_path: "pathogens/rvf/config/auspice_config.json"
       lat_longs_path: "pathogens/rvf/config/lat_longs.tsv"
       has_segments: true
   ```
3. The `auspice_config_path` is used as the base configuration file

### Dynamic Enhancement

For each build, the system enhances the base configuration with metadata-specific optimizations:

1. The `generate_dynamic_auspice_config` rule processes:

- Base configuration from the pathogen-specific JSON
- Actual metadata from the current dataset
- Settings from master_config.yaml

2. The enhancements include:

- Adding colorings for metadata fields present in the actual data
- Ensuring filters match available metadata fields
- Optimizing defaults based on dataset characteristics
- Adding genome segment displays for segmented pathogens
- Generating appropriate color scales for categorical data

3. The dynamically enhanced configuration is used in the export rule

# Adding a New Pathogen

1. To add support for a new pathogen:
   pathogens/
   └── your_pathogen/
   └── config/
   ├── config.yaml
   ├── auspice_config.json
   └── lat_longs.tsv

2. Add the pathogen to config/master_config.yaml:

```yaml
pathogens:
  your_pathogen:
    name: "Your Pathogen Scientific Name"
    config_path: "pathogens/your_pathogen/config/config.yaml"
    auspice_config_path: "pathogens/your_pathogen/config/auspice_config.json"
    lat_longs_path: "pathogens/your_pathogen/config/lat_longs.tsv"
    has_segments: false # Set to true if your pathogen has genome segments
```

3. Create a base `auspice_config.json` with pathogen-specific settings:

- Define core colorings relevant to your pathogen
- Add pathogen-specific annotation features
- Set appropriate defaults

4. Set `active_pathogen: "your_pathogen"` in `master_config.yaml` to analyze this pathogen

The metadata-driven system will handle the rest, adjusting visualizations based on available data.

# Customizing Dynamic Generation

To adjust how metadata drives configuration:

1. Edit `scripts/generate_auspice_config.py` to add custom logic for your pathogen
2. Add pathogen-specific settings in the `enhance_auspice_config` function
3. For complex customizations, add conditional logic based on the pathogen name

# Debugging

If you encounter issues with Auspice configuration:

1. Check the dynamically generated configuration in `results/configs/`
2. Review logs at `logs/generate_dynamic_auspice_config_*.log`
3. Test with a simple dataset to isolate if issues are data or configuration related

## How This Solution Works

The implemented solution provides a robust, metadata-driven approach to Auspice configuration by:

1. **Maintaining Pathogen-Specific Base Configurations**: Each pathogen has its own base configuration providing essential settings.

2. **Dynamic Enhancement During Build**: The script analyzes actual metadata and enhances the base configuration to reflect the specific dataset characteristics.

3. **Smart Defaults**: The system automatically adjusts visualization defaults based on dataset characteristics (e.g., setting segment coloring for segmented viruses).

4. **Seamless Integration**: The Snakemake workflow transparently incorporates this process, requiring no changes to how users run the pipeline.

This approach provides the flexibility to support multiple pathogens while ensuring each visualization is optimized for its specific data characteristics.
