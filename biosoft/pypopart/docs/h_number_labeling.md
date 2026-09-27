# Custom H Number Labeling

## Overview

The PyPopART Dash app allows users to customize the labels assigned to haplotype groups in network visualizations. By default, haplotypes are labeled sequentially as H1, H2, H3, etc. This feature enables you to replace these default labels with meaningful names for publication or analysis.

## Use Cases

- **Publication-ready figures**: Replace H1, H2, H3 with descriptive names like "Central", "Branch_A", "Rare_variant"
- **Geographic populations**: Label haplotypes by their predominant location (e.g., "Pacific", "Atlantic", "Mediterranean")
- **Functional groups**: Assign labels based on phenotype or functional characteristics
- **Standardization**: Use consistent labeling across multiple analyses

## Workflow

### Step 1: Generate the Network

1. Upload your sequence alignment (FASTA, NEXUS, or PHYLIP format)
2. Optionally upload metadata (population, location data)
3. Select a network algorithm (MST, MSN, TCS, or MJN)
4. Click "Compute Network"

The network will be displayed with default labels (H1, H2, H3, etc.)

### Step 2: Download Label Template

1. Navigate to the **Haplotype Summary** tab
2. Click the **"⬇️ Download Label Template"** button
3. Save the CSV file (default name: `h_number_mapping_template.csv`)

The template CSV contains two columns:

- `current_h_number`: The current H number (H1, H2, H3, etc.)
- `new_label`: Initially set to the same as current_h_number (for you to edit)

### Step 3: Edit Labels

Open the downloaded CSV in any spreadsheet application or text editor:

```csv
current_h_number,new_label
H1,H1
H2,H2
H3,H3
H4,H4
```

Edit the `new_label` column with your desired labels:

```csv
current_h_number,new_label
H1,Central
H2,Branch_A
H3,Branch_B
H4,Rare_variant
```

**Important Guidelines:**

- Keep the `current_h_number` column unchanged
- Do not use duplicate values in the `new_label` column
- Avoid special characters that might cause issues in visualization
- Keep labels concise for better readability in plots
- Do not leave `new_label` cells empty

### Step 4: Upload Custom Mapping

1. Return to the **Haplotype Summary** tab in the app
2. Click the **"⬆️ Upload Label Mapping"** button
3. Select your edited CSV file
4. The system will validate your mapping

**If successful:**

- A green success message appears
- The network graph updates with your custom labels
- The haplotype summary table shows your custom labels
- Search dropdown reflects your custom labels
- Tooltips display your custom labels

**If errors are detected:**

- A red error message appears with specific issues
- Fix the errors in your CSV and try uploading again

## Validation Rules

The system performs the following validations:

### 1. Column Format

- CSV must have exactly two columns: `current_h_number` and `new_label`
- Column headers must match exactly (case-sensitive)

### 2. Unknown H Numbers

- All `current_h_number` values must exist in the current network
- Error example: _"Row 5: Unknown H number 'H99'"_

### 3. Duplicate Labels

- Each `new_label` must be unique
- Error example: _"Duplicate label 'Central' for node1 and node2"_

### 4. Missing Values

- Both columns must have values for each row
- Error example: _"Row 3: Missing new_label for H2"_

### 5. Whitespace

- Leading and trailing whitespace is automatically trimmed
- Labels like `" Central "` become `"Central"`

## Example Validation Errors

```
❌ Validation Errors:
• Row 3: Missing new_label for H2
• Row 5: Unknown H number "H99"
• Duplicate label "Central" for H1 and H3
(3 total errors)
```

## Tips

1. **Start with the template**: Always download the current template before editing to ensure you have the correct H numbers
2. **Backup your CSV**: Keep a copy of your custom mapping for reuse
3. **Test with small changes**: Try uploading with just a few label changes first
4. **Use consistent naming**: Establish a naming convention for your labels
5. **Avoid unicode**: Stick to ASCII characters for maximum compatibility

## Persistence

- Custom labels persist across layout changes (e.g., switching from hierarchical to spring)
- Custom labels are **not** saved with the network export (they're display-only)
- To reset to default labels, compute a new network or reload the page
- To update labels, download a fresh template and upload a new mapping

## Integration with Other Features

- **Search**: Custom labels appear in the search dropdown
- **Tooltips**: Hovering over nodes shows custom labels
- **Summary Export**: Downloaded CSV includes custom labels
- **Network Export**: Exported GraphML/JSON uses original node IDs (not custom labels)

## Troubleshooting

**Problem**: Upload button doesn't respond

- **Solution**: Ensure you've computed a network first

**Problem**: Template download is empty

- **Solution**: Compute a network before downloading the template

**Problem**: Custom labels don't appear after upload

- **Solution**: Check the feedback message for validation errors

**Problem**: Labels are too long and overlap in the plot

- **Solution**: Use shorter labels or increase network spacing

**Problem**: Need to start over

- **Solution**: Refresh the page or compute a new network to reset labels

## API Reference

For programmatic access, see the callback documentation in `src/pypopart/gui/app.py`:

- `download_h_number_template()`: Generates template CSV
- `upload_h_number_mapping()`: Processes uploaded mapping

## Example Mappings

### Geographic Labeling

```csv
current_h_number,new_label
H1,Pacific_Core
H2,Atlantic_North
H3,Atlantic_South
H4,Mediterranean
```

### Functional Labeling

```csv
current_h_number,new_label
H1,Wildtype
H2,Resistant
H3,Tolerant
H4,Susceptible
```

### Alphanumeric Labeling

```csv
current_h_number,new_label
H1,A
H2,B
H3,C
H4,D
```
