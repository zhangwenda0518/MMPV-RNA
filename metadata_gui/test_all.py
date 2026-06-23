"""Quick integration test for metadata GUI."""
import sys, os, traceback

try:
    from models.data_store import MetadataStore
    store = MetadataStore()

    info = r'D:\桌面\延伸基因组\MMPV-RNA\public_metadata_pipeline\public_data_pipeline_output\info'
    store.load(os.path.join(info, 'Global_Unified_Metadata_Core13.tsv'))
    print(f'1. Loaded {store.row_count} rows')

    vals = store.get_unique_values('Tissue', 20)
    print(f'2. Unique Tissue: {vals[:5]}')

    r = store.search('leaf')
    print(f'3. Search leaf: {len(r)} results')

    s = store.get_stats()
    print(f'4. Stats: {s["total_rows"]} rows')

    from views.search_view import SearchPanel
    sp = SearchPanel(store)
    print('5. SearchPanel OK')

    from views.detail_panel import DetailEditPanel
    dp = DetailEditPanel(store)
    dp.load_record(0)
    print('6. DetailEditPanel OK')

    from views.viz_panel import VisualizationPanel
    vp = VisualizationPanel(store)
    print('7. VisualizationPanel OK')

    from views.metadata_table import MetadataTableView
    tv = MetadataTableView(store)
    print('8. MetadataTableView OK')

    from views.main_window import MainWindow
    from PySide6.QtWidgets import QApplication
    qapp = QApplication.instance()
    if qapp is None:
        qapp = QApplication(sys.argv)
    w = MainWindow()
    print(f'9. MainWindow OK, rows={w.store.row_count}')

    print('\n=== ALL TESTS PASSED ===')
except Exception:
    traceback.print_exc()
