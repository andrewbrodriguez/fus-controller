/**
 * List every image series in a whole-slide file, as JSON on stdout.
 *
 *   QuPath script -a <image.vsi> scripts/qupath/list_series.groovy
 *
 * Olympus .vsi files hold several series: overview/label images plus one
 * series per scanned section. This prints enough to tell them apart.
 */
import qupath.lib.images.servers.ImageServerProvider
import qupath.lib.io.GsonTools
import java.awt.image.BufferedImage

def path = args[0]
def support = ImageServerProvider.getPreferredUriImageSupport(BufferedImage, path)
def rows = []
support.getBuilders().eachWithIndex { builder, i ->
    def server = builder.build()
    def cal = server.getPixelCalibration()
    rows << [
        index      : i,
        name       : server.getMetadata().getName(),
        width      : server.getWidth(),
        height     : server.getHeight(),
        pixel_um   : cal.hasPixelSizeMicrons() ? cal.getAveragedPixelSizeMicrons() : null,
        channels   : server.getMetadata().getChannels().collect { it.getName() },
        pixel_type : server.getPixelType().toString(),
        z_slices   : server.nZSlices(),
        downsamples: server.getPreferredDownsamples() as List,
        args       : builder.getURIs().collect { it.toString() } + [builder.getArgs() as List],
    ]
    server.close()
}
println "SERIES_JSON " + GsonTools.getInstance().toJson(rows)
