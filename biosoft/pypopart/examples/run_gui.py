#!/usr/bin/env python
"""
Example script to launch the PyPopART GUI.

This script demonstrates how to start the PyPopART web interface.
"""

from pypopart.gui import main

if __name__ == '__main__':
    print('=' * 60)
    print('PyPopART GUI - Haplotype Network Analysis')
    print('=' * 60)
    print('\nStarting web server...')
    print('Once started, open your browser to: http://localhost:8050')
    print('\nTo stop the server, press Ctrl+C in this terminal')
    print('=' * 60)
    print()

    # Launch the GUI
    # Set debug=True for development, debug=False for production
    main(debug=True, port=8050)
