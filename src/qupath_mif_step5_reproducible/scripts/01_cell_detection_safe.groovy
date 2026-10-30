/**
 * QuPath 0.7.0 — STEP 01: SAFE single-run DAPI watershed detection.
 * Run ONLY in a FRESH SEPARATE QuPath PROJECT with reference ROI imported.
 * Protects any image with existing detections; never clears or overwrites them.
 * Source: 2026-10-30 command history, FINAL DETECTION (threshold 2.0).
 */

String ROOT = System.getenv('QUPATH_MIF_ROOT')
if (!ROOT) throw new IllegalStateException('Set QUPATH_MIF_ROOT to the local project directory')
final int EXPECTED_REFERENCE_CELLS = 1310

def imageData = getCurrentImageData()
if (imageData == null) {
    println "STOP: No image is open."
    return
}

def server = imageData.getServer()
def metadata = server.getMetadata()
def channels = metadata.getChannels().collect { it.getName() }
def cal = server.getPixelCalibration()

// Deliberately restricted to this LuCa image, not a general batch processor.
if (server.getWidth() != 1868 || server.getHeight() != 1400) {
    println "STOP: Wrong image dimensions. Expected 1868 x 1400."
    return
}
if (!channels.contains("DAPI") || channels.size() != 8) {
    println "STOP: Expected 8-channel LuCa image with DAPI. Found ${channels}"
    return
}
if (!cal.hasPixelSizeMicrons() ||
        Math.abs(cal.getPixelWidthMicrons() - 0.4979944789) > 0.005 ||
        Math.abs(cal.getPixelHeightMicrons() - 0.4979944789) > 0.005) {
    println "STOP: Pixel calibration differs from reference."
    return
}

// Critical safety guard: no existing detections of ANY kind.
int nExisting = getDetectionObjects().size()
if (nExisting != 0) {
    println "STOP: ${nExisting} detections already exist. Nothing was changed."
    println "Do not run detection on the original 1310-cell image."
    return
}

def selected = getSelectedObjects().toList()
if (selected.size() != 1 || !selected[0].isAnnotation()) {
    println "STOP: Select EXACTLY ONE imported reference annotation."
    println "Selected objects: ${selected.size()}"
    return
}
if (selected[0].getROI() == null || !selected[0].getROI().isArea()) {
    println "STOP: Selected annotation must have an area ROI."
    return
}

println "===== SAFE CELL DETECTION START ====="
println "Image: ${metadata.getName()}"
println "Annotation area px2: ${selected[0].getROI().getArea()}"
println "Before detection: ${nExisting} detections"
setImageType('FLUORESCENCE')

// Exactly one detection call, copied from the LAST runPlugin in the history.
runPlugin('qupath.imagej.detect.cells.WatershedCellDetection', '{"detectionImage":"DAPI","requestedPixelSizeMicrons":0.5,"backgroundRadiusMicrons":8.0,"backgroundByReconstruction":true,"medianRadiusMicrons":0.0,"sigmaMicrons":1.5,"minAreaMicrons":10.0,"maxAreaMicrons":200.0,"threshold":2.0,"watershedPostProcess":true,"cellExpansionMicrons":3.0,"includeNuclei":true,"smoothBoundaries":true,"makeMeasurements":true}')

int nCells = getCellObjects().size()
println "===== SAFE CELL DETECTION RESULT ====="
println "Cell count: ${nCells}"
println "Reference count: ${EXPECTED_REFERENCE_CELLS}"
println "Count match: ${nCells == EXPECTED_REFERENCE_CELLS}"
println "NOTE: Exact reproducibility also requires centroid-level matching to reference CSV."
println "NOTE: This script does NOT save the project automatically; inspect then Ctrl+S."
println "===== COMPLETE ====="
