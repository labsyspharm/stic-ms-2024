1. Download
   [ome-omero-roitool](https://github.com/glencoesoftware/ome-omero-roitool/releases/tag/v0.2.5)
   (v0.2.5 was used for testing `ome-omero-roitool-0.2.5.zip`) and de-compress
   the zip file. The executable is
   `ome-omero-roitool-0.2.5/bin/ome-omero-roitool`

1. Install [palom](https://github.com/labsyspharm/palom) 

1. Construct a configuration CSV file for running the mapping script
   (`map-geomx-roi.py` in the next step). Each line in the CSV will be run
   sequentially. For example, make a `pair.csv` file -

    ```csv
    geomx_path,cycif_path
    hits-server/path/to/moving-image/LSP16076_P54_A31_C100_HMS_Orion7@20230401_030302_784959.ome.tiff,hits-server/path/to/reference-image/LSP16076_P54_A31_C100_HMS_Orion7@20230401_030302_784959.ome.tiff
    ```

    The three columns shown above are:

    - `geomx_path`: filepath to the moving (GeoMx OME-TIFF) image
    - `cycif_path`: filepath to the reference image (Pyramidal OME-TIFF)

1. Run [`map-geomx-roi.py`
   script](https://github.com/labsyspharm/stic-ms-2024/blob/main/CYCIF_analysis/map-roi/map-geomx-roi.py)
   in the conda env where `palom` is installed. For example -

    ```bash
    python map-geomx-roi.py run-batch pairs.csv --out_dir hits-server/path/to/output/rois --refine --num_processes 2
    ```

    Note that the `--refine` flag isn't recommended for large ROI annotations
    nor sections that were not adjacent.

1. After running the script, one new XML file containing transformed ROIs will
   be generated for each line in the CSV file. QC images and logs will also be
   generated under the output directory specified. For example, the following
   XML is generated

    ```text
    LSP16076_P54_A31_C100_HMS_Orion7@20230401_030302_784959-rois-from-LSP16076@20240109_195658_063915-coarse.ome.xml
    ```

1. On OMERO, find the image ID of the reference image and import the mapped ROI.

    ```bash
    ome-omero-roitool import 6797 hits-server/path/to/output/rois/LSP16076_P54_A31_C100_HMS_Orion7@20230401_030302_784959-rois-from-LSP16076@20240109_195658_063915-coarse.ome.xml --server idp.tissue-atlas.org --key <SESSION-TOKEN>
    ```