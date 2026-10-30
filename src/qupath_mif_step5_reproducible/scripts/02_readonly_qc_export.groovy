/**
 * QuPath 0.7.0 — STEP 02: read-only report and CSV export.
 * Safe on both original image and separate reproduction project.
 * Does NOT modify detections or classifications.
 */
import java.text.SimpleDateFormat
import java.util.Locale

String ROOT = System.getenv('QUPATH_MIF_ROOT')
if (!ROOT) throw new IllegalStateException('Set QUPATH_MIF_ROOT to the local project directory')
def imageData = getCurrentImageData()
if (imageData == null) {
    println "STOP: No image is open."
    return
}
def cells = getCellObjects().toList()
if (cells.isEmpty()) {
    println "STOP: No detected cells to report."
    return
}

def cal = imageData.getServer().getPixelCalibration()
if (!cal.hasPixelSizeMicrons()) {
    println "STOP: No calibrated pixel size."
    return
}
double pxW = cal.getPixelWidthMicrons()
double pxH = cal.getPixelHeightMicrons()

def areas = cells.collect { c ->
    def roi = c.getNucleusROI()
    roi == null ? null : roi.getArea() * pxW * pxH
}.findAll { it != null }.sort()

def outDir = new File(ROOT, "results/step5/reproduction")
if (!outDir.exists() && !outDir.mkdirs()) {
    println "STOP: Cannot create output directory: ${outDir.absolutePath}"
    return
}

def stamp = new SimpleDateFormat("yyyyMMdd_HHmmss_SSS").format(new Date())
def out = new File(outDir, "cells_${cells.size()}_${stamp}.csv")

def fmt = { double x -> String.format(Locale.US, "%.8f", x) }

out.withWriter("UTF-8") { writer ->
    writer.writeLine("cell_id,x_um,y_um,nucleus_area_um2,CK,CD8,DAPI")
    cells.each { cell ->
        def roi = cell.getROI()
        def nucleus = cell.getNucleusROI()
        def m = cell.getMeasurementList()
        def line = [
            cell.getID().toString(),
            fmt(roi.getCentroidX() * pxW),
            fmt(roi.getCentroidY() * pxH),
            nucleus == null ? "" : fmt(nucleus.getArea() * pxW * pxH),
            fmt(m.get("Cell: CK (Opal 690) mean")),
            fmt(m.get("Cell: CD8 (Opal 540) mean")),
            fmt(m.get("Nucleus: DAPI mean"))
        ]
        writer.writeLine(line.join(","))
    }
}

println "===== REPRODUCTION QC ====="
println "Image type: ${imageData.getImageType()}"
println "Cells: ${cells.size()}"
println "Cells with nucleus ROI: ${areas.size()}"
if (!areas.isEmpty()) {
    println String.format(Locale.US, "Median nucleus area: %.4f um2", areas[areas.size().intdiv(2)])
}
println "CSV: ${out.absolutePath}"
println "No QuPath objects changed."
println "===== COMPLETE ====="
