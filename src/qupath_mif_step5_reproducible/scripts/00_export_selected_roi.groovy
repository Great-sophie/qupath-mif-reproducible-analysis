/**
 * QuPath 0.7.0 — STEP 00: export exact manually drawn analysis annotation.
 * READ ONLY with respect to image objects. Run on ORIGINAL project.
 * Select exactly ONE annotation first; do not select detected cells.
 */

String ROOT = System.getenv('QUPATH_MIF_ROOT')
if (!ROOT) throw new IllegalStateException('Set QUPATH_MIF_ROOT to the local project directory')
def imageData = getCurrentImageData()
if (imageData == null) {
    println "STOP: Open LuCa-7color image first."
    return
}

def selected = getSelectedObjects().toList()
if (selected.size() != 1 || !selected[0].isAnnotation()) {
    println "STOP: Select exactly ONE analysis annotation (yellow ROI), not cells."
    println "Selected objects: ${selected.size()}"
    return
}

def roi = selected[0].getROI()
if (roi == null || !roi.isArea()) {
    println "STOP: Selected annotation must be an AREA ROI."
    return
}

def outDir = new File(ROOT, "results/step5")
if (!outDir.exists() && !outDir.mkdirs()) {
    println "STOP: Unable to create ${outDir.absolutePath}"
    return
}

def output = new File(outDir, "original_roi_reference.geojson")
if (output.exists()) {
    println "STOP: ROI reference file already exists; refusing to overwrite it."
    println "Existing: ${output.absolutePath}"
    return
}

exportObjectsToGeoJson(selected, output.absolutePath, "FEATURE_COLLECTION")

def cal = imageData.getServer().getPixelCalibration()
def area_um2 = roi.getArea() * cal.getPixelWidthMicrons() * cal.getPixelHeightMicrons()

println "===== ORIGINAL ROI BACKUP ====="
println "Image: ${imageData.getServer().getMetadata().getName()}"
println "Selected ROI: ${roi.getClass().getSimpleName()}"
println "ROI bounding box px: x=${roi.getBoundsX()}, y=${roi.getBoundsY()}, width=${roi.getBoundsWidth()}, height=${roi.getBoundsHeight()}"
println String.format(java.util.Locale.US, "ROI area: %.6f mm2", area_um2 / 1_000_000.0)
println "Existing cells (not modified): ${getCellObjects().size()}"
println "GeoJSON: ${output.absolutePath}"
println "READ-ONLY EXPORT COMPLETE"
