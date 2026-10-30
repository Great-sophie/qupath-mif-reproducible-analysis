// QuPath 0.7: Export all detected cells, compartment-selected measurements,
// and autofluorescence to a separate CSV. Does not change QuPath objects.
import java.text.SimpleDateFormat
import java.util.Locale

def ROOT = System.getenv('QUPATH_MIF_ROOT')
if (!ROOT) throw new IllegalStateException('Set QUPATH_MIF_ROOT to the local project directory')
def imageData = getCurrentImageData()
if (imageData == null) {
    println 'ERROR: Open the LuCa-7color image before running this script.'
    return
}

def cells = getCellObjects().toList()
if (cells.isEmpty()) {
    println 'ERROR: No cell detections. Open the image with your 1310 saved cells.'
    return
}

def expected = [
    PDL1: 'Cell: PDL1 (Opal 520) mean',
    CD8: 'Cell: CD8 (Opal 540) mean',
    FoxP3: 'Nucleus: FoxP3 (Opal 570) mean',
    CD68: 'Cytoplasm: CD68 (Opal 620) mean',
    PD1: 'Cell: PD1 (Opal 650) mean',
    CK: 'Cell: CK (Opal 690) mean',
    DAPI: 'Nucleus: DAPI mean',
    Autofluorescence: 'Cell: Autofluorescence mean'
]

def names = cells[0].getMeasurementList().getNames()
def missing = expected.findAll { label, measurement -> !names.contains(measurement) }
if (!missing.isEmpty()) {
    println 'ERROR: Missing measurement columns, refusing to export misleading values:'
    missing.each { k, v -> println "  ${k} -> ${v}" }
    return
}

def cal = imageData.getServer().getPixelCalibration()
double sx = cal.getPixelWidthMicrons()
double sy = cal.getPixelHeightMicrons()
if (!Double.isFinite(sx) || !Double.isFinite(sy) || sx <= 0 || sy <= 0) {
    println 'ERROR: No valid physical pixel calibration.'
    return
}

// This file can be used alongside the full QuPath measurements export.
// A unique timestamp avoids silently overwriting prior results.
def outDir = new File(ROOT, 'results/day4')
if (!outDir.isDirectory() && !outDir.mkdirs()) {
    println 'ERROR: Could not create ' + outDir.absolutePath
    return
}
def stamp = new SimpleDateFormat('yyyyMMdd_HHmmss').format(new Date())
def output = new File(outDir, "LuCa7_full_AF_${cells.size()}_${stamp}.csv")
def fmt = { double value -> Double.isFinite(value) ? String.format(Locale.US, '%.8f', value) : '' }
def keys = expected.keySet().toList()

output.withWriter('UTF-8') { w ->
    w.writeLine((['cell_id','x_um','y_um','nucleus_area_um2'] + keys).join(','))
    cells.each { cell ->
        def roi = cell.getROI()
        def nroi = cell.getNucleusROI()
        double narea = nroi == null ? Double.NaN : nroi.getArea() * sx * sy
        def ml = cell.getMeasurementList()
        def line = [cell.getID().toString(),
                    fmt(roi.getCentroidX() * sx),
                    fmt(roi.getCentroidY() * sy),
                    fmt(narea)]
        keys.each { k -> line.add(fmt(ml.get(expected[k]))) }
        w.writeLine(line.join(','))
    }
}

println '===== FULL mIF + AF EXPORT ====='
println "Image type: ${imageData.getImageType()}"
println "Cells exported: ${cells.size()}"
println "AF measurement: ${expected.Autofluorescence}"
println "Pixel calibration: ${sx} x ${sy} um"
println "CSV: ${output.absolutePath}"
println 'No QuPath cells or classifications were changed.'
println '===== COMPLETE ====='
