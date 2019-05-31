import bioformats
import javabridge
import numpy as np

# img has shape (channels, height, width)
def write_ome_tiff(path, img, name, channel_names):
    omexml = bioformats.omexml.OMEXML()
    omexml.image(0).Name = name
    omexml.image(0).Pixels.PixelType = bioformats.omexml.PT_UINT16
    omexml.image(0).Pixels.DimensionOrder = bioformats.omexml.DO_XYCZT
    omexml.image(0).Pixels.SizeX = img.shape[2]
    omexml.image(0).Pixels.SizeY = img.shape[1]
    omexml.image(0).Pixels.SizeC = img.shape[0]
    omexml.image(0).Pixels.SizeZ = 1
    omexml.image(0).Pixels.SizeT = 1
    assert len(channel_names) == img.shape[0]
    omexml.image(0).Pixels.channel_count = len(channel_names)
    for i, channel_name in enumerate(channel_names):
        omexml.image(0).Pixels.Channel(i).Name = channel_name
        omexml.image(0).Pixels.Channel(i).SamplesPerPixel = 1
    javabridge.start_vm(class_path=bioformats.JARS)
    try:
        env = javabridge.get_env()
        buffer = env.make_object_array(img.shape[0], env.find_class('[B'))
        for i, channel_img in enumerate(img):
            channel_data = np.frombuffer(np.ascontiguousarray(channel_img, "<u2").data, np.uint8)
            env.set_object_array_element(buffer, i, env.make_byte_array(channel_data))
        javabridge.run_script(
            """
            importClass(Packages.loci.common.DebugTools);
            importClass(Packages.loci.common.services.ServiceFactory);
            importClass(Packages.loci.formats.services.OMEXMLService);
            importClass(Packages.loci.formats.ImageWriter);
            DebugTools.enableLogging("INFO");
            var service = new ServiceFactory().getInstance(OMEXMLService);
            var metadata = service.createOMEXMLMetadata(xml);
            var writer = new ImageWriter();
            writer.setMetadataRetrieve(metadata);
            writer.setWriteSequentially(true);
            writer.setId(path);
            for (var i = 0; i < buffer.length; ++i) {
                writer.saveBytes(i, buffer[i]);
            }
            writer.close();
            """,
            dict(path=path, xml=omexml.to_xml(), buffer=buffer)
        )
    finally:
        javabridge.kill_vm()
