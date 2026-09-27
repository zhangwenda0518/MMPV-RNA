# Custom H Number Labeling Example

This directory contains an example CSV file demonstrating the custom H number labeling feature in the PyPopART Dash GUI.

## File: h_number_mapping_example.csv

This example shows how to customize haplotype labels by replacing default H numbers (H1, H2, H3...) with meaningful names.

### Example Content

```csv
current_h_number,new_label
H1,Central
H2,Branch_A
H3,Branch_B
H4,Outlier
```

### Usage in Dash GUI

1. **Compute a Network**
   - Upload your sequence alignment (FASTA, NEXUS, or PHYLIP)
   - Select a network algorithm (MST, MSN, TCS, or MJN)
   - Click "Compute Network"
   - The network will display with default labels: H1, H2, H3, H4...

2. **Download Template**
   - Navigate to the "Haplotype Summary" tab
   - Click "⬇️ Download Label Template"
   - This generates a CSV with all current H numbers

3. **Edit Labels**
   - Open the downloaded CSV in Excel, Google Sheets, or a text editor
   - Modify the `new_label` column with your desired names
   - Save the file
   - Use the example above as a reference

4. **Upload Custom Mapping**
   - In the "Haplotype Summary" tab, click "⬆️ Upload Label Mapping"
   - Select your edited CSV file
   - The app will validate your mapping and show any errors
   - On success, the graph updates instantly with your custom labels

### Example Use Cases

#### Geographic Labels

```csv
current_h_number,new_label
H1,Pacific_Central
H2,Atlantic_North
H3,Atlantic_South
H4,Mediterranean
```

#### Functional/Phenotypic Labels

```csv
current_h_number,new_label
H1,Wildtype
H2,Resistant_A
H3,Resistant_B
H4,Susceptible
```

#### Simple Alphanumeric Labels

```csv
current_h_number,new_label
H1,A
H2,B
H3,C
H4,D
```

### Important Notes

- **Keep Column Names**: The CSV must have exactly two columns: `current_h_number` and `new_label`
- **No Duplicates**: Each `new_label` must be unique
- **Match H Numbers**: All `current_h_number` values must exist in your current network
- **No Empty Values**: Both columns must have values for each row
- **Whitespace**: Leading and trailing spaces are automatically trimmed

### Validation

The app performs automatic validation and will show specific error messages if:

- CSV has wrong column names
- H numbers don't exist in the network (e.g., H99 when you only have 10 haplotypes)
- New labels are duplicated
- Any values are missing

Example error message:

```
❌ Validation Errors:
• Row 3: Missing new_label for H2
• Row 5: Unknown H number "H99"
• Duplicate label "Central" for H1 and H3
(3 total errors)
```

### Persistence

- Custom labels persist across layout changes
- Custom labels appear in all UI elements (graph, tooltips, search, summary)
- Custom labels are included in CSV exports
- Custom labels reset when you compute a new network

### Getting Help

For detailed documentation, see:

- User Guide: `docs/h_number_labeling.md`
- Implementation Details: `H_NUMBER_LABELING_IMPLEMENTATION.md`

### Example Workflow

This example shows the complete workflow:

1. Start with a network of 4 haplotypes (H1-H4)
2. Download template:

   ```csv
   current_h_number,new_label
   H1,H1
   H2,H2
   H3,H3
   H4,H4
   ```

3. Edit to meaningful names:

   ```csv
   current_h_number,new_label
   H1,Central
   H2,Branch_A
   H3,Branch_B
   H4,Outlier
   ```

4. Upload the edited file
5. View the updated network with custom labels

The graph nodes will now display "Central", "Branch_A", "Branch_B", and "Outlier" instead of H1-H4.
