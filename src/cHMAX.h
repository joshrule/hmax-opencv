#ifndef __CHMAX__
#define __CHMAX__

#include "opencv2/opencv.hpp"  
#include <list>

using namespace std;

//=======================================================================================

static const int TYPE = CV_64F;
static const int debug = 0;
static const int jniDebug = 0;

//=======================================================================================

class CHMAX {
	//==============================================================================
	//                 		VARIABLES
	//==============================================================================
private:
	cv::Mat mC1SpaceSS; // spatial Space Steps: Vector defining the spatial pooling range for each scale band
	cv::Mat mC1ScaleSS; //scaling Scale Steps: Vector defining scale bands
	cv::Mat mPatchSizes;
	
    int mNumPatches; // 4075 in our tests, read from config file or computed from xml file data
	int mNumPatchSizes; // dfd: was double. Number of patch SIZES, = 4 in our tests (4, 8, 12, 16)
	int mNumScales;	// dfd: number of unique scales, = 17 (but 16 are used). Number of Gabor filters = mNumScales x mNumOrientations
	int mNumScaleBands;
	double mS1C1Suppress; 
	int mC1Overlap; // Defines overlap between C1 units. Not used
	
	// Filtering
	int mNumOrientations; // dfd: was called numSimpleFilters
	vector<cv::Mat> mGaborFilters; // vector of 2D square filter arrays
	int mNumFilters; // Number of Gabor filters, = mNumScales x mNumOrientations 
	vector<int> mFilterSizes; // non-unique filter sizes. dfd: Unique filter sizes should be a data member
	
	vector<cv::Mat> mC1Res; // dfd: Remove. Pass arguments instead
	vector<cv::Mat> mS2Target; // dfd: Remove. Pass arguments instead
	vector<vector<cv::Mat> > mPatchesVec; // vector of patches. dfd: OK to keep as data member because it is constant for large nb of images
	
	int mDoSkipSquareRoot; // dfd: flag to turn on use of exp(-minDistance^2) instead of minDistance to a patch as feature component for that patch
	
public:
	
	//==============================================================================
	//                 		 METHODS
	//==============================================================================
private:
	
	
public:
	//Constructor
	CHMAX() {
		cout << "HMAX default constructor- hmax parameters not fully initialized" << endl;
		mDoSkipSquareRoot = 0;
	};
	
	// Overloaded Constructor initializes necessary fields based 
	CHMAX(const cv::Mat& filters, const cv::Mat& fSiz) { 
		mDoSkipSquareRoot = 0;
		//Extract filters into 3D matrices
		CHMAX::extractFilters(filters,fSiz);
		
		cout << "Number of Filters: " << mNumFilters << endl;
		
		cout << "Filters initialized. Exiting hmax constructor" << endl;
  	};
	
	CHMAX(const string xmlFile) {
		mDoSkipSquareRoot = 0;
		CHMAX::initialize(xmlFile);
		cout << "doSkipSqrt: " << mDoSkipSquareRoot << endl;
		cout << "Exiting hmax constructor" << endl;
	};
	
	//Default Destructor
	~CHMAX() {};
	
	//==============================================================================================
	void initialize(const string xmlFile);
	void initializeWithoutPatches(const string xmlFile,  int numPatches); // dfd: ignore precomputed patches when reading xml file made by Matlab
	
	//===================== Wrapper ========================
	void runHMAX(const vector<cv::Mat>& imgs, cv::Mat& features); //Vector of images
	void runHMAX(const cv::Mat& img, cv::Mat& features);  //Single image version
	
	//===================== S1 Layer =======================
	// Translates Matlab: function [c1res,s1res] = c1(stim, filters, fSiz, c1SpaceSS, c1ScaleSS, c1OL,s1c1suppress,INCLUDEBORDERS)
	void procS1C1(const cv::Mat& inIm);
	void extractFilters(const cv::Mat& filters, const cv::Mat& filterSizes); //Convert filters from column vectors to 3D matrices
	
	void makeS1Layer(const cv::Mat& inIm, vector<cv::Mat>& s1Layer);  // dfd: Version without use of internal object states
	void makeC1Layer(const vector<cv::Mat>& s1Layer, vector<cv::Mat>& c1Layer);  // dfd: Version without use of variable internal object states
	
	//getSimpleFilterReponse: S1 Layer filter response
	void getSimpleFilterResponse(const cv::Mat& inIm, const cv::Mat& s1Norm,int iUFilterIndex, cv::Mat& response, bool REMOVEBORDERS = false); 
    
    // Streamlined version of getSimpleFilterResponse
    void findFilterResponse(const cv::Mat& inIm, const cv::Mat& s1Norm, int filterIndex, cv::Mat& normalizedFilteredIm);
    
	void addC1ScaleSS(int width); //Adds another band with specified width
	
	void makePatchesVec(const vector<cv::Mat>& imgs, vector<vector<cv::Mat> >& patchesVec);  // dfd: make random patches instead of reading them
	
	//===================== S2 Layer =======================
	void procS2C2(cv::Mat& c2res_out);
	void extractPatch(const int iS2, const int iCenter, cv::Mat& Patch);
	void extractAllPatches();
	
	void makeS2C2Features(vector<cv::Mat>& c1Layer, cv::Mat& c2Features); // dfd: Version without use of variable internal object states
	
	// GETTERS & SETTERS //
	vector<int> getFilterSizes() {return mFilterSizes;};
	cv::Mat getC1SpaceSS(){ return mC1SpaceSS; };
	cv::Mat getC1ScaleSS() {return mC1ScaleSS; };
	int getNumScaleBands(){	return mNumScaleBands; };
	int getNumScales(){ return mNumScales; }; 
	int getGlobalFilterIndex(int iBand, int iScale, int iFilt); //Ordering is important  
	vector<cv::Mat> getC1Res() {return mC1Res;};   
	
	int getNumOrientations(){return mNumOrientations;}; 
	int getNumPatches() {return mNumPatches;};
	void getScalesInThisBand(int band, vector<int>& bandScales);
	
	void setC1SpaceSS(const cv::Mat c1Space) {mC1SpaceSS = c1Space;};
	void setC1ScaleSS(const cv::Mat c1Scale) {mC1ScaleSS = c1Scale;};
	void setSuppressThresh(const cv::Mat suppressThreshMat) {mS1C1Suppress = suppressThreshMat.at<double>(0,0);};
	void setPatchSizes(const cv::Mat patchSizes) {mPatchSizes = patchSizes;};
	void setPatchesVec(const vector<vector<cv::Mat> > patchesVec) {mPatchesVec = patchesVec;};
	void setSkipSqrtFlag(int doSkipSquareRoot) {mDoSkipSquareRoot = doSkipSquareRoot;};
	//	void setS2Target(const cv::Mat s2Tgt) {s2Target.push_back(s2Tgt.clone());} // dfd: not used: s2Target is now filled in initialize function
};

//==============================================================================================
//==============================================================================================

#endif

