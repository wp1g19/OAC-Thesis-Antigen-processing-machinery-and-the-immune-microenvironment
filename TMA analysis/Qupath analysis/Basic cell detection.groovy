setImageType('BRIGHTFIELD_H_DAB');
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.65111 0.70119 0.29049", "Stain 2" : "DAB", "Values 2" : "0.26917 0.56824 0.77759", "Background" : " 255 255 255"}');
selectAnnotations();
runPlugin('qupath.opencv.CellCountsCV', '{"stainChannel": "Hematoxylin + DAB",  "gaussianSigmaMicrons": 1.0,  "backgroundRadiusMicrons": 10.0,  "doDoG": false,  "threshold": 0.05,  "thresholdDAB": 0.1,  "detectionDiameter": 25.0}');
