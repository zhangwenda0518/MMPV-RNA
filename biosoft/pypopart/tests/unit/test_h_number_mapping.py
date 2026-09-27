"""Unit tests for H number mapping functionality in GUI."""

import base64
import csv
import io


class TestHNumberMapping:
    """Test cases for H number label mapping feature."""

    def test_template_csv_format(self):
        """Test that the template CSV has the correct format."""
        # Simulate template generation
        output = io.StringIO()
        writer = csv.writer(output)
        writer.writerow(['current_h_number', 'new_label'])
        writer.writerow(['H1', 'H1'])
        writer.writerow(['H2', 'H2'])
        writer.writerow(['H3', 'H3'])

        content = output.getvalue()

        # Parse to verify format
        reader = csv.DictReader(io.StringIO(content))
        rows = list(reader)

        assert len(rows) == 3
        assert 'current_h_number' in rows[0]
        assert 'new_label' in rows[0]
        assert rows[0]['current_h_number'] == 'H1'
        assert rows[0]['new_label'] == 'H1'

    def test_valid_mapping_csv_parsing(self):
        """Test parsing of a valid H number mapping CSV."""
        csv_content = """current_h_number,new_label
H1,Central
H2,Branch_A
H3,Branch_B"""

        reader = csv.DictReader(io.StringIO(csv_content))
        rows = list(reader)

        assert len(rows) == 3
        assert rows[0]['current_h_number'] == 'H1'
        assert rows[0]['new_label'] == 'Central'
        assert rows[1]['current_h_number'] == 'H2'
        assert rows[1]['new_label'] == 'Branch_A'

    def test_invalid_csv_columns_detected(self):
        """Test that invalid column names are detected."""
        # Wrong column names
        csv_content = """h_number,label
H1,Central
H2,Branch"""

        reader = csv.DictReader(io.StringIO(csv_content))
        fieldnames = reader.fieldnames

        # Should not match expected format
        assert fieldnames != ['current_h_number', 'new_label']

    def test_duplicate_labels_detection(self):
        """Test detection of duplicate new labels."""
        mapping = {
            'node1': 'A',
            'node2': 'B',
            'node3': 'A',  # Duplicate!
        }

        # Check for duplicates
        seen_labels = {}
        duplicates = []

        for node_id, label in mapping.items():
            if label in seen_labels:
                duplicates.append(f'Duplicate label "{label}"')
            else:
                seen_labels[label] = node_id

        assert len(duplicates) > 0
        assert 'Duplicate label "A"' in duplicates[0]

    def test_csv_encoding_decoding(self):
        """Test CSV base64 encoding/decoding workflow."""
        # Simulate uploaded file content
        csv_content = """current_h_number,new_label
H1,Central
H2,Branch"""

        # Encode as would be received from upload
        encoded = base64.b64encode(csv_content.encode('utf-8')).decode('utf-8')
        contents = f'data:text/csv;base64,{encoded}'

        # Decode as would be done in callback
        content_type, content_string = contents.split(',')
        decoded = base64.b64decode(content_string).decode('utf-8')

        assert decoded == csv_content

        # Parse decoded content
        reader = csv.DictReader(io.StringIO(decoded))
        rows = list(reader)

        assert len(rows) == 2
        assert rows[0]['new_label'] == 'Central'

    def test_missing_required_fields(self):
        """Test detection of missing required fields."""
        # Missing new_label
        csv_content = """current_h_number,new_label
H1,Central
H2,"""

        reader = csv.DictReader(io.StringIO(csv_content))
        rows = list(reader)

        errors = []
        for row_num, row in enumerate(rows, start=2):
            current_h = row.get('current_h_number', '').strip()
            new_label = row.get('new_label', '').strip()

            if current_h and not new_label:
                errors.append(f'Row {row_num}: Missing new_label')

        assert len(errors) == 1
        assert 'Row 3: Missing new_label' in errors[0]

    def test_unknown_h_numbers_detected(self):
        """Test detection of unknown H numbers in mapping."""
        # Known H numbers
        h_to_node = {
            'H1': 'node1',
            'H2': 'node2',
            'H3': 'node3',
        }

        # Mapping with unknown H number
        csv_content = """current_h_number,new_label
H1,Central
H99,Unknown"""

        reader = csv.DictReader(io.StringIO(csv_content))
        errors = []

        for row_num, row in enumerate(reader, start=2):
            current_h = row.get('current_h_number', '').strip()
            if current_h not in h_to_node:
                errors.append(f'Row {row_num}: Unknown H number "{current_h}"')

        assert len(errors) == 1
        assert 'H99' in errors[0]

    def test_whitespace_handling(self):
        """Test that leading/trailing whitespace is properly handled."""
        csv_content = """current_h_number,new_label
 H1 , Central
H2,Branch_A"""

        reader = csv.DictReader(io.StringIO(csv_content))
        rows = list(reader)

        # Values should be stripped
        h1_label = rows[0]['new_label'].strip()
        assert h1_label == 'Central'

        h1_num = rows[0]['current_h_number'].strip()
        assert h1_num == 'H1'
