/**
 * Export several full-resolution regions of one series as OME-TIFF, in one QuPath launch.
 *
 *   QuPath script -a <image.vsi> -a <series> -a x,y,w,h,out.ome.tif [-a x,y,w,h,out2.ome.tif ...] \
 *       scripts/qupath/export_regions.groovy
 *
 * Same output as `QuPath convert-ome --series=N -r x,y,w,h -c ZLIB`, but the image is
 * opened once for all regions instead of once per region (fus.cells uses it for the
 * six target ROIs of a section).
 */
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.images.writers.ome.OMEPyramidWriter
import qupath.lib.regions.ImageRegion
import java.awt.image.BufferedImage

def path = args[0]
def series = args[1] as int
def support = ImageServerProvider.getPreferredUriImageSupport(BufferedImage, path)
def server = support.getBuilders()[series].build()
for (spec in args[2..-1]) {
    def parts = spec.split(",", 5)
    def (x, y, w, h) = parts[0..3].collect { it as int }
    def region = ImageRegion.createInstance(x, y, w, h, 0, 0)
    new OMEPyramidWriter.Builder(server)
            .region(region)
            .downsamples(1.0d)
            .compression(OMEPyramidWriter.CompressionType.ZLIB)
            .parallelize()
            .build()
            .writePyramid(parts[4])
    println "WROTE " + parts[4]
}
server.close()
