clearDetections();
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.57568 0.58068 0.57568 ", "Stain 2" : "DAB", "Values 2" : "0.39716 0.58723 0.70528 ", "Background" : " 255 255 255 "}');
selectAnnotations();
runPlugin('qupath.imagej.superpixels.DoGSuperpixelsPlugin', '{"downsampleFactor": 1.0,  "sigmaPixels": 10.0,  "minThreshold": 10.0,  "maxThreshold": 230.0,  "noiseThreshold": 1.0}');
selectDetections();
runPlugin('qupath.lib.algorithms.IntensityFeaturesPlugin', '{"downsample": 1.0,  "region": "ROI",  "tileSizePixels": 200.0,  "colorOD": false,  "colorStain1": true,  "colorStain2": true,  "colorStain3": false,  "colorRed": false,  "colorGreen": false,  "colorBlue": false,  "colorHue": false,  "colorSaturation": false,  "colorBrightness": false,  "doMean": true,  "doStdDev": true,  "doMinMax": true,  "doMedian": true,  "doHaralick": false,  "haralickDistance": 1,  "haralickBins": 32}');
