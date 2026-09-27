"""Unit tests for geographic metadata handling."""

import pytest

from pypopart.io.metadata import (
    extract_coordinates,
    parse_coordinate,
    validate_latitude,
    validate_longitude,
)


class TestCoordinateParsing:
    """Test coordinate parsing and validation."""

    def test_parse_decimal_coordinate(self):
        """Test parsing decimal coordinate."""
        assert parse_coordinate('45.5') == 45.5
        assert parse_coordinate('-123.4') == -123.4

    def test_parse_coordinate_with_degree_symbol(self):
        """Test parsing coordinate with degree symbol."""
        assert parse_coordinate('45.5°') == 45.5
        assert parse_coordinate('-123.4°') == -123.4

    def test_parse_coordinate_with_whitespace(self):
        """Test parsing coordinate with whitespace."""
        assert parse_coordinate(' 45.5 ') == 45.5
        assert parse_coordinate(' -123.4° ') == -123.4

    def test_parse_invalid_coordinate(self):
        """Test error handling for invalid coordinates."""
        with pytest.raises(ValueError, match='Invalid coordinate'):
            parse_coordinate('invalid')
        with pytest.raises(ValueError, match='Invalid coordinate'):
            parse_coordinate('N45.5')

    def test_validate_latitude_valid(self):
        """Test validation of valid latitude values."""
        validate_latitude(0)
        validate_latitude(45.5)
        validate_latitude(-45.5)
        validate_latitude(90)
        validate_latitude(-90)

    def test_validate_latitude_invalid(self):
        """Test validation of invalid latitude values."""
        with pytest.raises(ValueError, match='Latitude must be between'):
            validate_latitude(91)
        with pytest.raises(ValueError, match='Latitude must be between'):
            validate_latitude(-91)

    def test_validate_longitude_valid(self):
        """Test validation of valid longitude values."""
        validate_longitude(0)
        validate_longitude(123.4)
        validate_longitude(-123.4)
        validate_longitude(180)
        validate_longitude(-180)

    def test_validate_longitude_invalid(self):
        """Test validation of invalid longitude values."""
        with pytest.raises(ValueError, match='Longitude must be between'):
            validate_longitude(181)
        with pytest.raises(ValueError, match='Longitude must be between'):
            validate_longitude(-181)


class TestExtractCoordinates:
    """Test coordinate extraction from metadata."""

    def test_extract_coordinates_success(self):
        """Test successful coordinate extraction."""
        metadata = {'latitude': '45.5', 'longitude': '-123.4', 'other': 'data'}
        coords = extract_coordinates(metadata)
        assert coords == (45.5, -123.4)

    def test_extract_coordinates_missing(self):
        """Test extraction when coordinates are missing."""
        metadata = {'other': 'data'}
        coords = extract_coordinates(metadata)
        assert coords is None

    def test_extract_coordinates_partial(self):
        """Test extraction when only one coordinate is present."""
        metadata = {'latitude': '45.5', 'other': 'data'}
        coords = extract_coordinates(metadata)
        assert coords is None

    def test_extract_coordinates_custom_columns(self):
        """Test extraction with custom column names."""
        metadata = {'lat': '45.5', 'lon': '-123.4'}
        coords = extract_coordinates(metadata, lat_column='lat', lon_column='lon')
        assert coords == (45.5, -123.4)

    def test_extract_coordinates_invalid_values(self):
        """Test extraction with invalid coordinate values."""
        metadata = {'latitude': '95', 'longitude': '-123.4'}
        with pytest.raises(ValueError):
            extract_coordinates(metadata, validate=True)

    def test_extract_coordinates_no_validation(self):
        """Test extraction without validation."""
        metadata = {'latitude': '95', 'longitude': '-123.4'}
        coords = extract_coordinates(metadata, validate=False)
        assert coords == (95, -123.4)
