Anatomo-functional mapping of the barrel cortex: Step by step user guide

Toolbox to reproduce the numerical results of the paper:

L. Perronnet, M.E. Vilarchao, G. Hucher, D.E. Shulz, G. Peyré, I. Ferezou. [An automated workflow for the anatomo-functional mapping of the barrel cortex](https://hal.archives-ouvertes.fr/hal-01196436). Journal of Neuroscience Methods, 263(1), pp. 145–154, 2016.

![Example of registration of barrels](img/registration.png)

Copyright (c) 2016 L. Perronnet, M.E. Vilarchao, G. Hucher, D.E. Shulz, G. Peyré, I. Ferezou


=================
Requirements

Matlab with its Statistics Toolbox. (I don’t know if the code needs the image processing toolbox or not, it was the case at one point but maybe it was for the segmentation that has finally been removed... I don’t know how to check this…).

=================
Get started

Open the file Run_registration_GUI.m with Matlab. By running the code, you will open the first page of the graphic user interface (GUI).

=================
Load images

Clicking on the “load images” button will allow you to browse your computer to choose the set of images you want to work on and to load the corresponding files.
Images will then appear as thumbnails at the bottom part of the window. A simple click on an image triggers its display in the main window. Checkboxes allow you to select the images you want to use for the subsequent steps. Click on the “Next” button to continue.

=================
Segmentation

First check the number of “regions” included in your images (2 if the image contains only the histological slice on a homogeneous background, 3 if the background contains an additional black region in the periphery due to the optical configuration). A click on the “segment” button will trigger the segmentation which can take several seconds. Click on the “Next” button to continue.

=================
Registration

The left control panel allows you to select the couple of images you want to register and to adjust the parameters (a precise description of the parameters can be found in the full text article: An automated workflow for the anatomo-functional mapping of the barrel cortex). The initialization “Angles” relate to the relative rotations of the images initially tested by the algorithm. To run the registration, click the “Register!” button.
Single images are displayed on the right (top), together with the result of the registration (bottom).
If you want to select a region of interest (ROI) to restrict the analysis on the barrel field, first click on the button indicated by a square, draw a rectangle on its corresponding image and double click on it to validate the ROI.
Two additional tools appear next to the second image on the right. The upper one allows you to reinitialize the image after an unsuccessful registration attempt (in this case it is recommended to repeat the operation with another set of parameters). The other one gives you the opportunity to flip the image in case the histological slice would have been mounted in the wrong orientation.
Next to the display of the registration result, a slider allows you to change the transparency of the upper image.
Following the registration of all the pairs of images, you can explore the results by selecting a chosen pair from the upper left list. Finally, you can select (bottom right) a target folder in which the transformed images will be stored, and click on the “next” button to continue.

=================
Fusion

Use the checkboxes to select the images which contain some apparent barrels and click on the “Merge” button. The fusion process might take some time, but after about a minute the result should appear as an image on the bottom left. You can then set the post-processing parameters and click on the “Post-process” button.
Finally, you can set the target folder where the two images will be stored, and click on the “Save data” button.
