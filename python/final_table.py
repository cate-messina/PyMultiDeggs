import pandas as pd

import contextlib
import os
import tempfile
import webbrowser

_BOOL_COLS = ["is_stimulation", "consensus_stimulation", "is_inhibition",
              "consensus_inhibition", "is_directed", "consensus_direction"]

_PAGE_TEMPLATE = """<html><head><meta charset="utf-8">
<style>
table { border-collapse: collapse; font-family: sans-serif; font-size: 14px; }
th, td { border: 1px solid #ccc; padding: 6px 10px; text-align: left; }
th { background-color: #f2f2f2; position: sticky; top: 0; }
tr:nth-child(even) { background-color: #fafafa; }
</style></head><body>__TABLE__</body></html>"""


def _wrap_html_page(table_html):
    """Wrap a pre-rendered HTML table inside a complete styled page.

    Args:
        table_html (str): HTML string for the table body.

    Returns:
        str: Full HTML document containing the table and default styling.
    """
    return _PAGE_TEMPLATE.replace("__TABLE__", table_html)


def _truncate_pmdi_html(value, max_ids=5):
    """Shorten long semicolon-separated PMDI values for HTML display.

    Args:
        value: Raw PMDI interaction value.
        max_ids (int, optional): Maximum number of IDs to display before the
            rest is collapsed. Defaults to 5.

    Returns:
        str: The original value if it is short enough; otherwise, a collapsible
        HTML summary block.
    """
    if pd.isna(value) or str(value).strip() == "":
        return ""
    ids = str(value).split(";")
    if len(ids) <= max_ids:
        return value
    preview = "; ".join(ids[:max_ids])
    return (f"<details><summary>{preview} ... (+{len(ids)-max_ids} more)</summary>"
            f"{value}</details>")


class MultiDEGGsTable(pd.DataFrame):
    """DataFrame subclass with enhanced HTML rendering for results.

    This class behaves like a standard pandas DataFrame, while converting
    boolean columns to ``True``/``False`` strings and collapsing lengthy
    ``PMDI_interaction`` values into clickable details blocks in Jupyter.
    """
    _metadata = []

    @property
    def _constructor(self):
        return MultiDEGGsTable

    def _to_truefalse(self):
        df = pd.DataFrame(self).copy()
        for col in _BOOL_COLS:
            if col in df.columns:
                df[col] = df[col].astype(int).map({1: "True", 0: "False"})
        return df

    def _repr_html_(self):
        df = self._to_truefalse()
        if "PMDI_interaction" in df.columns:
            df["PMDI_interaction"] = df["PMDI_interaction"].apply(_truncate_pmdi_html)
        return df.to_html(escape=False, index=False)

    def save(self, path="multideggs_results.html"):
        """Save the full table as a standalone HTML document.

        Args:
            path (str, optional): Destination path for the generated file.
                Defaults to "multideggs_results.html".

        Returns:
            str: The path where the HTML file was saved.
        """
        df = self._to_truefalse()  # No truncation: full table
        table_html = df.to_html(escape=False, index=False)
        full_page = _wrap_html_page(table_html)
        with open(path, "w", encoding="utf-8") as f:
            f.write(full_page)
        print(f"Complete table saved to: {path}")
        return path


class MultiDEGGsTableDict(dict):
    """Dictionary wrapper for multi-omics result tables.

    This behaves like a regular dictionary while also providing HTML rendering
    and export support for concatenated category-level tables.
    """

    def _combined(self):
        frames = []
        for cat, df in self.items():
            if isinstance(df, pd.DataFrame) and not df.empty:
                tmp = pd.DataFrame(df).copy()
                tmp.insert(0, "category", cat)
                for col in _BOOL_COLS:
                    if col in tmp.columns:
                        tmp[col] = tmp[col].astype(int).map({1: "True", 0: "False"})
                frames.append(tmp)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()

    def _repr_html_(self):
        df = self._combined()
        if "PMDI_interaction" in df.columns:
            df["PMDI_interaction"] = df["PMDI_interaction"].apply(_truncate_pmdi_html)
        return df.to_html(escape=False, index=False)

    def save(self, path="multideggs_results.html"):
        """Save the full concatenated multi-omics table as HTML.

        Args:
            path (str, optional): Destination path for the generated file.
                Defaults to "multideggs_results.html".

        Returns:
            str: The path where the HTML file was saved.
        """
        df = self._combined()  # No truncation: full table
        table_html = df.to_html(escape=False, index=False)
        full_page = _wrap_html_page(table_html)
        with open(path, "w", encoding="utf-8") as f:
            f.write(full_page)
        print(f"Complete table saved to: {path}")
        return path


def show_html_results(final_outputs):
    """Open the full results table in the browser using a temporary HTML file.

    The table is rendered as a standalone HTML page and opened in the default
    browser without writing the output to a user-selected location. This is
    intended for local Jupyter sessions where a browser is available on the
    same machine.

    Args:
        final_outputs: A DataFrame, a tuple of ``(results, final_outputs)``, or a
            dictionary keyed by category.

    Returns:
        None: Opens the HTML output in the browser.
    """

    # Accept either the full tuple (results, final_outputs) or final_outputs alone.
    if isinstance(final_outputs, tuple) and len(final_outputs) == 2:
        final_outputs = final_outputs[1]

    if isinstance(final_outputs, pd.DataFrame):
        df = pd.DataFrame(final_outputs).copy()
        for col in _BOOL_COLS:
            if col in df.columns:
                df[col] = df[col].astype(int).map({1: "True", 0: "False"})
        table_html = df.to_html(escape=False, index=False)
    else:
        frames = []
        for cat, df in final_outputs.items():
            if isinstance(df, pd.DataFrame) and not df.empty:
                tmp = df.copy()
                tmp.insert(0, "category", cat)
                for col in _BOOL_COLS:
                    if col in tmp.columns:
                        tmp[col] = tmp[col].astype(int).map({1: "True", 0: "False"})
                frames.append(tmp)
        combined = pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
        table_html = combined.to_html(escape=False, index=False)

    full_page = _wrap_html_page(table_html)

    # Temporary file. This is needed for webbrowser.open() can open a real file more reliably than a raw data: URL.
    with tempfile.NamedTemporaryFile(mode="w", suffix=".html",
                                      delete=False, encoding="utf-8") as f:
        f.write(full_page)
        temp_path = f.name

    # Suppress harmless browser warnings 
    with open(os.devnull, "w") as devnull:
        with contextlib.redirect_stderr(devnull):
            webbrowser.open(f"file://{temp_path}")