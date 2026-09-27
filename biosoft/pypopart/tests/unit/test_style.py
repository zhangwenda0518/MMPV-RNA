"""Tests for the shared palette and population colour assignment."""

from pathlib import Path

import pytest

from pypopart.visualization.style import (
    MIN_NODE_LUMA,
    PALETTE,
    POP_ART_PALETTE,
    generate_population_colors,
    luma,
)

THEME_CSS = (
    Path(__file__).resolve().parents[2]
    / 'src'
    / 'pypopart'
    / 'gui'
    / 'assets'
    / 'theme.css'
)


class TestPaletteDoesNotDrift:
    """CSS cannot import the Python constants, so pin them together."""

    def test_theme_css_ships(self):
        """The stylesheet has to exist for Dash to auto-serve it."""
        assert THEME_CSS.is_file()

    @pytest.mark.parametrize('color', PALETTE)
    def test_palette_colour_is_in_the_stylesheet(self, color):
        """Every palette colour is declared as a CSS custom property."""
        assert color in THEME_CSS.read_text().lower()

    @pytest.mark.parametrize('color', POP_ART_PALETTE)
    def test_population_colours_are_distinct(self, color):
        """Two populations sharing a colour would be unreadable."""
        assert POP_ART_PALETTE.count(color) == 1


class TestGeneratePopulationColors:
    """Assigning colours to populations."""

    def test_uses_the_palette_while_it_lasts(self):
        """Small studies get the hand-picked colours."""
        colors = generate_population_colors(['PopB', 'PopA'])

        # Sorted, so PopA takes the first palette entry.
        assert colors['PopA'] == POP_ART_PALETTE[0]
        assert colors['PopB'] == POP_ART_PALETTE[1]

    def test_exactly_palette_length_still_uses_palette(self):
        """The boundary case must not tip into the HSV fallback."""
        pops = [f'Pop{i:02d}' for i in range(len(POP_ART_PALETTE))]
        colors = generate_population_colors(pops)

        assert set(colors.values()) == set(POP_ART_PALETTE)

    def test_falls_back_beyond_the_palette(self):
        """Cycling would repeat a colour, so switch to spaced hues."""
        pops = [f'Pop{i:02d}' for i in range(len(POP_ART_PALETTE) + 1)]
        colors = generate_population_colors(pops)

        assert len(set(colors.values())) == len(pops)

    @pytest.mark.parametrize('n', [1, 2, 5, 8, 9, 20, 50])
    def test_colours_are_always_unique(self, n):
        """The guarantee that matters, at every size."""
        colors = generate_population_colors([f'Pop{i:02d}' for i in range(n)])

        assert len(set(colors.values())) == n
        assert all(c.startswith('#') for c in colors.values())

    def test_empty_input(self):
        """No populations, no colours."""
        assert generate_population_colors([]) == {}


class TestNodeColoursAreLightEnough:
    """Node labels are ink on the fill, so dark fills are unreadable."""

    @pytest.mark.parametrize('color', POP_ART_PALETTE)
    def test_palette_clears_the_floor(self, color):
        """Every hand-picked population colour is light enough."""
        assert luma(color) >= MIN_NODE_LUMA

    @pytest.mark.parametrize('n', [1, 8, 16, 17, 25, 50, 120])
    def test_generated_colours_clear_the_floor(self, n):
        """Including the HSV fallback, whose blues used to come out dark."""
        colors = generate_population_colors([f'Pop{i:03d}' for i in range(n)])

        assert colors
        for color in colors.values():
            assert luma(color) >= MIN_NODE_LUMA, f'{color} is too dark'

    @pytest.mark.parametrize(
        'color,expected',
        [('#ffffff', 255.0), ('#000000', 0.0)],
    )
    def test_luma_endpoints(self, color, expected):
        """Sanity-check the measure itself."""
        assert luma(color) == pytest.approx(expected)

    @pytest.mark.parametrize('bad', ['', 'nonsense', '#fff', '#gggggg'])
    def test_unparseable_reads_as_dark(self, bad):
        """Failing closed means a bad value can never pass a light check."""
        assert luma(bad) == 0.0
