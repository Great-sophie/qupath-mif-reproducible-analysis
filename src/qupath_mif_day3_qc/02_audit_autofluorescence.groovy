// Read-only QuPath 0.7 inspection of autofluorescence measurements.
// Run with LuCa7 image open and detected cells present.
def cells = getCellObjects()
if (cells.isEmpty()) {
    println 'ERROR: No cells open or detected.'
    return
}
def names = cells[0].getMeasurementList().getNames()
def af = names.findAll { it.toLowerCase().contains('autofluorescence') }
println '===== AUTOFLUORESCENCE AUDIT ====='
println "Cells: ${cells.size()}"
println "AF measurement columns: ${af.size()}"
af.each { println "  ${it}" }
if (af.isEmpty()) {
    println 'WARNING: AF measurement not found; check all 125 measurement names and channel naming.'
    return
}
af.each { name ->
    def xs = cells.collect { it.getMeasurementList().get(name) }
        .findAll { Double.isFinite(it) }.sort()
    if (xs.size() > 0) {
        def q = { double p -> xs[(int)Math.round((xs.size() - 1) * p)] }
        println String.format('AF: %s N=%d P05=%.4f P50=%.4f P95=%.4f',
            name, xs.size(), q(0.05), q(0.5), q(0.95))
    }
}
println '===== DONE ====='
