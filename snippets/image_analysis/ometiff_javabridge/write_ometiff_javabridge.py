import bioformats
import javabridge
import numpy as np

# img has shape (channels, height, width)
def write_ome_tiff(path, img, name, fluor_names, channel_names=None, pixel_type=bioformats.omexml.PT_UINT16):
    FIELD_FLUOR = 'Fluor'
    omexml = bioformats.omexml.OMEXML()
    omexml.image(0).Name = name
    omexml.image(0).Pixels.PixelType = pixel_type
    omexml.image(0).Pixels.DimensionOrder = bioformats.omexml.DO_XYCZT
    omexml.image(0).Pixels.SizeX = img.shape[2]
    omexml.image(0).Pixels.SizeY = img.shape[1]
    omexml.image(0).Pixels.SizeC = img.shape[0]
    omexml.image(0).Pixels.SizeZ = 1
    omexml.image(0).Pixels.SizeT = 1
    assert len(fluor_names) == img.shape[0]
    if channel_names is None:
        channel_names = fluor_names
    assert (len(fluor_names) == len(channel_names))
    omexml.image(0).Pixels.channel_count = len(fluor_names)
    for i, (fluor_name, channel_name) in enumerate(zip(fluor_names, channel_names)):
        omexml.image(0).Pixels.Channel(i).Name = channel_name
        omexml.image(0).Pixels.Channel(i).SamplesPerPixel = 1
        omexml.image(0).Pixels.Channel(i).node.set(FIELD_FLUOR, fluor_name)
    javabridge.start_vm(class_path=bioformats.JARS)
    try:
        env = javabridge.get_env()
        buffer = env.make_object_array(img.shape[0], env.find_class('[B'))
        for i, channel_img in enumerate(img):
            channel_data = np.frombuffer(np.ascontiguousarray(channel_img).data, np.uint8)
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
