"""Package exactly what the Google Drive code folder needs -> code_for_drive.zip.

The notebook runs from this folder as-is (locally, or in Google Colab after
uploading the whole folder): it imports markov.py and viz.py and reads
data/mumbai_daily_rain.csv. Figures are written to figures/.
"""

import os
import zipfile

FILES = [
    "analysis.ipynb", "markov.py", "viz.py", "data/mumbai_daily_rain.csv",
    "download_data.py", "prepare_data.py", "requirements.txt",
]

README = """Mumbai monsoon Markov chain - code

Run analysis.ipynb from top to bottom (Jupyter or Google Colab).
It needs markov.py, viz.py and data/mumbai_daily_rain.csv in the same folder.
Install packages with:  pip install -r requirements.txt

The random seed is fixed, so every number in the report is reproduced exactly.
download_data.py + prepare_data.py rebuild the CSV from the original IMD and NOAA
sources (about 1 GB download); they are not needed to run the notebook.
"""

with zipfile.ZipFile("code_for_drive.zip", "w", zipfile.ZIP_DEFLATED) as z:
    for f in FILES:
        z.write(f, f"mumbai_monsoon_markov/{f}")
    z.writestr("mumbai_monsoon_markov/figures/.keep", "")
    z.writestr("mumbai_monsoon_markov/README.txt", README)
print("code_for_drive.zip:", os.path.getsize("code_for_drive.zip") // 1024, "KB")
