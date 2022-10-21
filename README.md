# Snake Model based Tongue Segmentation

<b>Update (September 23, 2017)This repository is the MATLAB implementation of our paper [Adaptive active contour model based automatic tongue image segmentation](http://ieeexplore.ieee.org/document/7852933/) which contributes to the improvement on active contour model by designing a self-adaptive coefficient and a searching region.

<br/>

## Procedure:
- Extract the rough binary image of the tounge based on Otsu's threshold method, wavelet filtering, thresholding, Morphological operator and kmeans.

- Extract and order the initial boundary from the binary image extracted before.

- Boundary refinement based on modified snake model by iteration