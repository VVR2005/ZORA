"""
Example ZORA plugin: GC Content Calculator
Demonstrates the plugin system.
Drop .py files into the plugins/ directory with a register_plugin() function.
"""

def register_plugin():
    return {
        'name': 'GC Stats Calculator',
        'version': '1.0',
        'description': 'Calculates detailed GC statistics for loaded sequences',
        'callback': run_gc_stats,
    }

def run_gc_stats():
    """Called when user activates this plugin from menu."""
    from PySide6.QtWidgets import QMessageBox, QApplication
    from zora_main import SeqUtils

    # Find main window
    app = QApplication.instance()
    if not app:
        return

    # Try to find the main window
    main_windows = [w for w in app.topLevelWidgets() if hasattr(w, 'sequences')]
    if not main_windows or not main_windows[0].sequences:
        QMessageBox.information(None, "GC Stats", "No sequences loaded. Load sequences first.")
        return

    win = main_windows[0]
    lines = ["GC Statistics Report", "=" * 40]
    for rec in win.sequences:
        seq = rec.sequence
        gc = SeqUtils.gc_content(seq)
        skew = SeqUtils.gc_skew(seq)
        lines.append(f"\n{rec.name}:")
        lines.append(f"  Length: {len(seq)} bp")
        lines.append(f"  GC%: {gc:.2f}%")
        lines.append(f"  GC Skew: {skew:.4f}")
        lines.append(f"  G count: {seq.upper().count('G')}")
        lines.append(f"  C count: {seq.upper().count('C')}")

    QMessageBox.information(None, "GC Stats Results", '\n'.join(lines))
