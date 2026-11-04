**Guide to Navigating Files**


_Igor Pro_

Code written in the Igor Pro programming language is stored in .ipf format (Igor Procedure Files) and must be opened in Igor Pro to be read. Most analysis functions used for the dissertation are contained within MaxwellFreemanFunctions.ipf, namely RPCropAnalysis2(), RPDCF(), and RPStats(). Unrelated and helper functions are also contained in that file. Other helper functions, mostly DeleteNaNs() and Regrid_Param2() are contained in AlexTools_v16.5, which belongs to Dr. J. Alex Huffman.


_R_

Code written in R was primarily used for application of the random forest classification and is stored in .R files. Note that the files are not designed to simply be run in their entirety--sections are designated for separate purposes and lines from those sections are run individually as necessary. REBS_RF.R corresponds to the classification of Raman spectral data; RPFM_RF_2026-11-02.R corresponds to classification of fluorescence data.


_Python_

Python 3 code is in .py files. As mentioned in Appendix A of the dissertation, pre-processing scripts necessary to convert raw images from the REBS into spectra belongs to Dr. David Doughty from the Army Research Laboratory and is not published here. Contact Dr. David Doughty or Dr. J. Alex Huffman (University of Denver) with inquiries into those scripts. Supplementary scripts written by Maxwell Freeman and the graphical user-interface script for the novel fluorescence microspectrometer are included here in log_metrics.py and PiCameraGUI_02-08-22.py, respectively.
