// script to detect P-cadherin positive cells and score them as high (3+), medium (2+) 
// and low (1+) based on the average DAB signal assocaited with the cell membrane 

setImageType('BRIGHTFIELD_H_DAB');
clearDetections();
setColorDeconvolutionStains('{"Name" : "H-DAB default", "Stain 1" : "Hematoxylin", "Values 1" : "0.68055 0.60778 0.40921 ", "Stain 2" : "DAB", "Values 2" : "0.24903 0.44956 0.85784 ", "Background" : " 149 112 85 "}');
selectAnnotations();
runPlugin('qupath.imagej.detect.cells.WatershedCellMembraneDetection', '{"detectionImageBrightfield": "Hematoxylin",  "requestedPixelSizeMicrons": 0.1,  "backgroundRadiusMicrons": 6.0,  "medianRadiusMicrons": 0.0,  "sigmaMicrons": 1.5,  "minAreaMicrons": 5.0,  "maxAreaMicrons": 500.0,  "threshold": 0.1,  "maxBackground": 4.0,  "watershedPostProcess": true,  "excludeDAB": true,  "cellExpansionMicrons": 5.0,  "limitExpansionByNucleusSize": false,  "includeNuclei": true,  "smoothBoundaries": true,  "makeMeasurements": true}');

// Delete objects that do not have a nucleus
def noNuclei = getCellObjects().findAll {it.getNucleusROI() == null}
removeObjects(noNuclei, true)

high = getPathClass('high')
medium = getPathClass('medium')
low = getPathClass('low')
negative = getPathClass('negative')

for (cell in getCellObjects()) {
    ch1 = measurement(cell, 'Membrane: DAB OD mean')

    if (ch1 > 1.0)
        cell.setPathClass(high)
        
    if (ch1 > 0.6 && ch1 <=1.0)
        cell.setPathClass(medium)
    
    if (ch1 > 0.3 && ch1 <=0.6)
        cell.setPathClass(low)
    
    if (ch1 <=0.3)
        cell.setPathClass(negative)
}
fireHierarchyUpdate()