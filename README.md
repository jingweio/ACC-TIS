## Adaptive Active Contour Model based Automatic Tongue Image Segmentation

Code for CISP-BMEI 2016 paper [Adaptive active contour model based automatic tongue image segmentation](http://ieeexplore.ieee.org/document/7852933/), which contributes to the improvement on active contour model by designing a self-adaptive coefficient and a searching region.


## Overall Procedures
- Extract the rough binary image of the tounge based on Otsu's threshold method, wavelet filtering, thresholding, Morphological operator and kmeans.
- Extract and order the initial boundary from the binary image extracted before.
- Boundary refinement based on modified snake model by iteration


## Citation
```
@INPROCEEDINGS{7852933,
  author={Guo, Jingwei and Yang, Yikang and Wu, Qingwei and Su, Jionglong and Ma, Fei},
  booktitle={2016 9th International Congress on Image and Signal Processing, BioMedical Engineering and Informatics (CISP-BMEI)}, 
  title={Adaptive active contour model based automatic tongue image segmentation}, 
  year={2016},
  volume={},
  number={},
  pages={1386-1390},
  doi={10.1109/CISP-BMEI.2016.7852933}}
```
