/*
 *  cHMAX.cpp
 *  HMAX classifier
 *
 *  Created by Rodriguez, Pedro A., March 2012 -- July 2012
 *  Copyright 2012 Pedro A. Rodriguez and JHU/APL. All rights reserved.
 *
 */
#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h>

#if 0
#include <list>
#include <map>
#include <sys/time.h>
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#endif

#include "cHMAX.h"
#include "Filters.h"
#include "mouseAndTimer.hpp" // dfd Sept. 2012

#include "opencv2/opencv.hpp"

using namespace std;
using namespace cv;

static const float maxFloat = std::numeric_limits<float>::max();
static const double maxDouble = std::numeric_limits<double>::max();
static const int debugPts = 20;

void removeBorders(const cv::Mat& src, cv::Mat& dst);

//==============================================================================
//             Main processing path for S1,C1 Layer
//==============================================================================
void CHMAX::runHMAX(const vector<cv::Mat>& imgs, cv::Mat& features){
	// This method runs through the C1,C2 layers of the HMAX classifier. It requires
	// vector of grayscale images (could be a single image, but must be a vector).
	// It returns a matrix of feature vectors where each row is a feature vector for
	// the corresponding image.     
	//	cv::Mat img;
	cv::Mat output; // dfd: feature vector for single image. Poor name choice...
	cv::Mat tempFeat;
    	MyTimer hmaxTimer;
	for (int i=0; i<imgs.size(); i++){
        	hmaxTimer.start();

		vector<cv::Mat> s1Layer, c1Layer;
		cv::Mat c2Features;
		makeS1Layer(imgs[i], s1Layer);
		makeC1Layer(s1Layer, c1Layer);
		makeS2C2Features(c1Layer, c2Features);
		features.push_back(c2Features);

        hmaxTimer.stop();        
        cout << "Finding features for image " << i+1 << " took " << hmaxTimer.time() << " seconds" << endl;
	}
}

void CHMAX::runHMAX(const cv::Mat& img, cv::Mat& features){
	// This method runs through the C1,C2 layers of the HMAX classifier. It requires
	// vector of grayscale images (could be a single image, but must be a vector).
	// It returns a matrix of feature vectors where each row is a feature vector for
	// the corresponding image.     
	//	cv::Mat img;
	cv::Mat output; // dfd: feature vector for single image. Poor name choice...
	cv::Mat tempFeat;
    	MyTimer hmaxTimer;
	
    	hmaxTimer.start();

	vector<cv::Mat> s1Layer, c1Layer;
	cv::Mat c2Features;
	makeS1Layer(img, s1Layer);
	makeC1Layer(s1Layer, c1Layer);
	makeS2C2Features(c1Layer, c2Features);
	features.push_back(c2Features);

	hmaxTimer.stop();        
	cout << "Finding features for image took " << hmaxTimer.time() << " seconds" << endl;
	
}

//==============================================================================
//             Main processing path for S1,C1 Layer
//==============================================================================
void CHMAX::procS1C1(const cv::Mat& inIm)
{
	cv::FileStorage fs;
	if (debug) {
		cout << "Processing in C1 Layer" << endl;
		fs.open("c1Debug.xml", cv::FileStorage::WRITE); 
	}
	cv::Mat img;
    	inIm.convertTo(img,TYPE);

	if (debug) cout << "Num Bands, Filter Orientations: " << mNumScaleBands << "," << mNumOrientations << endl;
	
	// Get image squared
	//    cv::Mat imsq = inIm.clone();
	//	imsq.convertTo(imsq,TYPE);  
    	cv::Mat imsq;
	imsq = img.mul(img);
	
	if (debug) cout << "IMSQ (r,c,ch,tp): " << imsq.rows << "," << imsq.cols << "," << imsq.channels() << "," << imsq.type() << endl;
	
	
	//******************** S1 Layer ***********************
	if (debug) cout << "Starting S1 Layer" << endl;
	// Initialize and compute s1Norm //
	cv::Mat s1Norm;
	cv::Mat zeroElmt;	
	vector<cv::Mat> s1NormVector;
	s1NormVector.resize(mNumFilters);
	if (debug) cout << "Computing S1 Norm" << endl;
	//-------- TESTING --------
	vector<int> uniqueFilterSizes(mFilterSizes); // dfd: Copy mFilterSizes vector to uniqueFilterSizes
	if (debug) cout << "Original fSize: " << uniqueFilterSizes.size() << endl;
	sort(uniqueFilterSizes.begin(), uniqueFilterSizes.end());
	vector<int>::iterator new_end = unique(uniqueFilterSizes.begin(), uniqueFilterSizes.end());
	uniqueFilterSizes.erase(new_end, uniqueFilterSizes.end());
	if (debug) cout << "Unique fSize: " << uniqueFilterSizes.size() << endl;
	
	for (int iFilter=0; iFilter < uniqueFilterSizes.size(); iFilter++){ //TODO: This loop may be excessive, ok for now
		// Precalculate normalizations
		cv::Scalar radius(0,0,0,(uniqueFilterSizes[iFilter]-1)/2);
		sumfilter(&imsq,radius,&s1Norm);
		// -- Remove small zero values
		cv::Mat tempmask;
		double thresh = 1e-14;
		cv::compare(s1Norm,thresh,tempmask,cv::CMP_GE);
		cv::divide(tempmask,255.0,tempmask);
		tempmask.convertTo(tempmask,TYPE);
		cv::multiply(s1Norm,tempmask,s1Norm);
		// ----
		if (debug && (iFilter == 1)) fs << "sumfiltS1" << s1Norm;
		cv::sqrt(s1Norm,s1Norm);
		
		// Avoid divide by zero
		cv::compare(s1Norm,0,zeroElmt,cv::CMP_EQ);
		cv::divide(zeroElmt,255.0,zeroElmt);
		zeroElmt.convertTo(zeroElmt,TYPE);
		s1Norm = s1Norm + zeroElmt;
		
		// Add current filter s1Norm to vector
		s1Norm.copyTo(s1NormVector[uniqueFilterSizes[iFilter]]);
		if (debug && (iFilter == 1)) fs << "s1Norm" << s1Norm;
	}
	zeroElmt.release();   
	
	if (debug) cout << "Unique s1Norm: " << s1NormVector.size() << endl;
	if (debug){
		cv::Mat s1n;
		s1NormVector[0].copyTo(s1n);    
		
	}
	
	// Loop over all bands and scales
	vector<int> bandScales; //Will be ordered list of scales in particular band
	cv::Mat mins,maxs,mins2,maxs2;
	cv::Mat response,mask; //Initialize for thresholding later
	cv::Mat suppressthresh = cv::Mat::zeros(img.rows,img.cols,TYPE);
	vector<cv::Mat> s1res;
	
	//Loop over bands
	for (int iBand=0; iBand<mNumScaleBands; iBand++){ 
		if (bandScales.size() > 0) bandScales.clear();
		getScalesInThisBand(iBand,bandScales); // dfd: Only used here to get the bandScales size, which is: 2.
		
		//Loop over scales within band
		for(int iScale=0; iScale< bandScales.size(); iScale++) {
			
			//Loop over orientations
			for (int iFilter=0; iFilter<mNumOrientations; iFilter++){
				int iUFilterIndex = getGlobalFilterIndex(iBand,iScale,iFilter);
				// Retrieve s1Norm for this filter
				s1NormVector[mFilterSizes[iUFilterIndex]].copyTo(s1Norm);
				getSimpleFilterResponse(img,s1Norm,iUFilterIndex,response);
				
				// Update filter response mins/maxs
				if (iFilter == 0) {
					maxs = response.clone();
					mins = response.clone();
				} else {                    
					cv::min(mins,response,mins2);
					mins2.copyTo(mins); 
					cv::max(maxs,response,maxs2);
					maxs2.copyTo(maxs);
				}
				
				if (debug && iBand == 0 && iScale == 1 && iFilter < 6) {
					cout << iUFilterIndex << ",";
				}
				s1res.push_back(response); // dfd: move this next to getSimpleFilterResponse. Remove clone and release
										   //				response.release();
			} //end for simple filter
			
			if (mS1C1Suppress > 0){ //For this band/scale only
				cv::Mat s1res_work;
				suppressthresh = mins + mS1C1Suppress*(maxs-mins);
				int start = getGlobalFilterIndex(iBand,iScale,0);
				
				for (int i=start; i< start+mNumOrientations; i++) {
					s1res[i].copyTo(s1res_work);
			
					cv::compare(s1res_work,suppressthresh,mask,cv::CMP_GE); //Mask has values of 0 or 255
					cv::divide(mask,255.0,mask);
					mask.convertTo(mask,TYPE);
					cv::multiply(s1res_work,mask,s1res[i]); // Elements with values less than threshold get 0
															//s1res[i] = s1res_work.clone(); //TODO: Is this right?                         
				}                    
			}
			maxs.release();
			mins.release();
		} //end for scale
	} //end for band
	
	if (debug) cout << "s1res: " << s1res.size() << endl;
	
	//******************** C1 Layer ***********************
	if (debug) cout << "Starting C1 Layer" << endl; //cv::waitKey(-1);
	vector<cv::Mat> tempC1;
	mC1Res.clear(); // Make sure that mC1Res is an empty vector
	cv::Mat s1Mat, c1Mat, mmfMat, c1Mat2;
	cv::Mat c1MatMC;
	cv::Mat c1mins, c1mins2, c1maxs, c1maxs2;     
	
	//Loop over scale bands
	for (int iBand=0; iBand<mNumScaleBands; iBand++){
		double temp = mC1SpaceSS.at<double>(0,iBand);
		cv::Scalar poolRange;
		poolRange[3] = (int) temp;
		
		//Loop over orientations
		for (int iFilter=0; iFilter<mNumOrientations; iFilter++){ // dfd: numOrientations (12)
			if (bandScales.size() > 0) bandScales.clear();
	     	getScalesInThisBand(iBand,bandScales);
			
			//Get first filter in scale band for c1Mat initialization               
			int iUFilter = getGlobalFilterIndex(iBand,0,iFilter);
			c1Mat = cv::Mat::zeros(s1res[iUFilter].rows,s1res[iUFilter].cols,TYPE);
			
			//Loop over band scales
			for(int iScale=0; iScale< (int)bandScales.size(); iScale++) {
				iUFilter = getGlobalFilterIndex(iBand,iScale,iFilter);
				s1res[iUFilter].copyTo(s1Mat); 
				cv::max(c1Mat,s1Mat,c1Mat2); // dfd: We are doing this operation twice. We did the same thing after calculating response. 
				c1Mat = c1Mat2.clone();
			} //end for scale
			
			if (debug && iFilter == 0 & iBand == 0) cout << "numscales: " << bandScales.size() << endl;
			mymaxfilter(c1Mat,poolRange,mmfMat); // dfd: for each orientation, find scale and space max           
			
			if (iFilter == 0) {
				c1mins = mmfMat.clone(); //cv::Mat::ones(mmfMat.rows,mmfMat.cols,TYPE)*(maxDouble);
				c1maxs = mmfMat.clone(); //cv::Mat::ones(mmfMat.rows,mmfMat.cols,TYPE)*(-maxDouble);
			} else {
				cv::min(c1mins,mmfMat,c1mins2); //Will be minimum across all filters
				c1mins = c1mins2.clone();
				cv::max(c1maxs,mmfMat,c1maxs2); //Will be maximum across all filters
				c1maxs = c1maxs2.clone();
			}
			tempC1.push_back(mmfMat);    //To be used for suppression
			
		} //end for simple filters
		
		//////// Suppression (within band) /////////
		if (mS1C1Suppress > 0){
			cv::Mat tempC1_work;
			if (debug) cout << "Suppression..." << endl;
			suppressthresh = cv::Mat::zeros(mmfMat.rows,mmfMat.cols,CV_64FC(mNumOrientations));
			suppressthresh = c1mins + mS1C1Suppress*(c1maxs-c1mins);
					
			for (int iFilter=0; iFilter<mNumOrientations; iFilter++){   
				tempC1[iFilter].copyTo(tempC1_work);                 
				cv::compare(tempC1_work,suppressthresh,mask,cv::CMP_GE); //Mask has values of 0 or 255
				cv::divide(mask,255.0,mask); //Convert to 0 or 1
				mask.convertTo(mask,TYPE);
				cv::multiply(tempC1_work,mask,tempC1_work); // Elements with values less than threshold get 0
				tempC1[iFilter] = tempC1_work.clone();
			}		
		} 
		
		cv::merge(tempC1,c1MatMC);
		mC1Res.push_back(c1MatMC);
		tempC1.clear(); //Should be empty for processing the next band
	} //end for band
	
	if (debug) cout << "mC1Res: " << mC1Res.size() << endl;
	//**************** Memory management ******************
	//Vectors
	tempC1.clear();
	s1res.clear();
	//Mats
#if 0
	imsq.release();
	suppressthresh.release();
	mask.release();
	mins.release();
	maxs.release();
	c1mins.release();
	c1maxs.release();
	response.release();
#endif
	fs.release();
	
} //end procS1C1

//==============================================================================
//             Main processing path for S2,C2 Layer
//==============================================================================
void CHMAX::procS2C2(cv::Mat& c2Res_out) {//using vector version
										  // Inputs:
										  //mC1Res = is a std vector of 3D cv::Mats
										  //mS2Target is a std vector of 3D matrix patches
										  //patchSize=Indicates the patch sizes
	
	//Estimate relevant dimensions of c1res and s2Target:
	int nbBands = mC1Res.size();
	int nbOrientations = mC1Res[0].channels(); //nfiltes need to be identical to s2Target[0].size.p[2]; 
	
	if (mPatchesVec[0][0].depth() != mC1Res[0].depth()) {
        cout << "The types of mC1Res and mPatchVec must be the same." << endl;//make this an error later. dfd: exit
        exit(1);
    }
	
	
	if (debug) cout << "nbBands: " << nbBands << endl;
	if (debug) cout << "nbOrientations: " << nbOrientations << endl;
	
	int totalNbPatches = 0;
	//    int nbPatchSizes = mPatchVec.size();
	for (int i=0; i<mNumPatchSizes; i++) {
        int nbPatchesPerSize = mPatchesVec[i].size();
        totalNbPatches += nbPatchesPerSize;
	}
    
	cv::Mat tempC2res(1,totalNbPatches,TYPE);
	cv::Mat tempcols;
	
	int curind, lastind = 0; //For forming one feature vector
	
    // dfd: we use mPatchVec dimensions instead of s2Target dimensions because s2Target exists only if we read patches from xml file, 
    // mS2Target is not defined when we create patches from calib images
	for (int iPatchSize = 0; iPatchSize < mNumPatchSizes; iPatchSize++) { 
		int nbPatchesPerSize = mPatchesVec[iPatchSize].size();
		cv::Mat c2Res(1,nbPatchesPerSize,TYPE,0.0);
		
		for (int iCenter = 0; iCenter < nbPatchesPerSize; iCenter++) {
			//extract serialized patch from s2Target and reshape:
			vector<cv::Mat> patch;
			
			cv::split(mPatchesVec[iPatchSize][iCenter], patch);
			
			for (int iBand= 0; iBand < nbBands; iBand++) {
				cv::Mat dst_vec_out;
				double minVal, maxVal;	
				cv::Point minLoc,maxLoc;
				
				vector<cv::Mat> c1Vec;
				cv::split(mC1Res[iBand],c1Vec);
				WindowedPatchDistance(c1Vec, patch, dst_vec_out, mDoSkipSquareRoot);
				cv::minMaxLoc(dst_vec_out, &minVal, &maxVal, &minLoc, &maxLoc);
	
                if (iBand == 0){
					c2Res.at<double>(0,iCenter) = minVal;
				}
				
				if (minVal < c2Res.at<double>(iCenter)) {
					c2Res.at<double>(0,iCenter) = minVal;
				}
			} //end for iBand
		} // for iCenter
		curind = lastind;
		lastind += nbPatchesPerSize;
		tempcols = tempC2res.colRange(curind,lastind); //The assignment is just a Mat header, not the column data.
		c2Res.copyTo(tempcols);		
    } // for iPatchSize
	
	tempC2res.copyTo(c2Res_out);
} //end procS2C2

//========================================================================================================


//==============================================================================
//  Functionality of procS1C1 using two functions, makeS1Layer and makeC1Layer 
//  where the S1 Layer and C1 Layer are explicit arguments
//==============================================================================

void CHMAX::makeS1Layer(const cv::Mat& inIm, vector<cv::Mat>& s1Layer)
{
	cv::FileStorage fs;

	if (debug) fs.open("s1Debug.xml", cv::FileStorage::WRITE); 
	
	cv::Mat img;
    inIm.convertTo(img,TYPE);
	
	// Get image squared
	cv::Mat imsq;
	imsq = img.mul(img);
	
	if (debug) cout << "IMSQ (r,c,ch,tp): " << imsq.rows << "," << imsq.cols << "," << imsq.channels() << "," << imsq.type() << endl;
	
	
	//******************** S1 Layer ***********************
	if (debug) cout << "Starting S1 Layer" << endl;
	// Initialize and compute s1Norm //
	cv::Mat s1Norm;
	cv::Mat zeroElmt;	
	vector<cv::Mat> s1NormVector;
	s1NormVector.resize(mNumFilters); 
	if (debug) cout << "Computing S1 Norm" << endl;
	//-------- TESTING --------
	vector<int> uFiltSizes(mFilterSizes); // dfd: Copy mFilterSizes vector to uFiltSizes
	if (debug) cout << "Original fSize: " << uFiltSizes.size() << endl;
	sort(uFiltSizes.begin(), uFiltSizes.end());
	vector<int>::iterator new_end = unique(uFiltSizes.begin(), uFiltSizes.end());
	uFiltSizes.erase(new_end, uFiltSizes.end());
	if (debug) cout << "Unique fSize: " << uFiltSizes.size() << endl;
	
	for (int iFilter=0; iFilter < uFiltSizes.size(); iFilter++){
		// Precalculate normalizations
		cv::Scalar radius(0,0,0,(uFiltSizes[iFilter]-1)/2);
		sumfilter(&imsq,radius,&s1Norm);
		// -- Remove small zero values
		cv::Mat tempmask;
		double thresh = 1e-14;
		cv::compare(s1Norm,thresh,tempmask,cv::CMP_GE);
		cv::divide(tempmask,255.0,tempmask);
		tempmask.convertTo(tempmask,TYPE);
		cv::multiply(s1Norm,tempmask,s1Norm);
		// ----
		if (debug && (iFilter == 1)) fs << "sumfiltS1" << s1Norm.clone();
		cv::sqrt(s1Norm,s1Norm);
		
		// Avoid divide by zero
		cv::compare(s1Norm,0,zeroElmt,cv::CMP_EQ);
		cv::divide(zeroElmt,255.0,zeroElmt);
		zeroElmt.convertTo(zeroElmt,TYPE);
		s1Norm = s1Norm + zeroElmt;
		
		// Add current filter s1Norm to vector
		s1Norm.copyTo(s1NormVector[uFiltSizes[iFilter]]);
		if (debug && (iFilter == 1)) fs << "s1Norm" << s1Norm;
	}
    //	zeroElmt.release();   
	
	if (debug) cout << "Unique s1Norm: " << s1NormVector.size() << endl;
	if (debug){
		cv::Mat s1n;
		s1NormVector[0].copyTo(s1n);    
    }
	
	// Loop over all bands and scales
	cv::Mat mins,maxs,mins2,maxs2;
	cv::Mat response0, response, mask; //Initialize for thresholding later
	cv::Mat suppressthresh = cv::Mat::zeros(img.rows,img.cols,TYPE);
    
	//Loop over bands
	for (int iBand=0; iBand<mNumScaleBands; iBand++){ 
        //		if (bandScales.size() > 0) bandScales.clear();
		vector<int> bandScales; // dfd: define locally so you don't have to clear it
		getScalesInThisBand(iBand,bandScales); // dfd: Only used here to get the bandScales size, which is: 2.
		
		//Loop over scales within band
		for(int iScale=0; iScale< bandScales.size(); iScale++) {
			//Loop over simple filters
			for (int iFilter=0; iFilter<mNumOrientations; iFilter++){  // dfd: numOrientations
				int iUFilterIndex = getGlobalFilterIndex(iBand,iScale,iFilter);
				// Retrieve s1Norm for this filter
				s1NormVector[mFilterSizes[iUFilterIndex]].copyTo(s1Norm);
#if 1    
				getSimpleFilterResponse(img,s1Norm,iUFilterIndex, response);
#else 
                // dfd: comment out #else and switch response0 and response to compare results of getSimpleFilterResponse and condensed function findFilterResponse. 
                // dfd: Mystery: same responses for getSimpleFilterResponse and findFilterResponse, but different features are found at the end, unless we push_back response.clone()
                findFilterResponse(img, s1Norm, iUFilterIndex, response);
                cv::Scalar debugSumAbsDiff = cv::sum(abs(response-response0));
                cout << "Sum of abs differences between response and response0: " << debugSumAbsDiff[0] << endl;
#endif
				s1Layer.push_back(response.clone()); // dfd: Part of mystery solved. response and response.clone give different results!				
													 //				s1Layer.push_back(response); // dfd: moved this next to getSimpleFilterResponse. Removed clone and release				
													 // Update filter response mins/maxs
#if 0
				if (iFilter == 0) {
					maxs = response.clone();
					mins = response.clone();
				} else {                    
					cv::min(mins,response,mins2);
					mins2.copyTo(mins); 
					cv::max(maxs,response,maxs2);
					maxs2.copyTo(maxs);
				}
#else
				if (iFilter == 0) {
					maxs = response.clone();
					mins = response.clone();
				} else {                    
					cv::min(mins,response,mins);
					cv::max(maxs,response,maxs);
				}
#endif
				if (debug && iBand == 0 && iScale == 1 && iFilter < 6) {
					//fs << "s1Norm" << s1Norm;
					cout << iUFilterIndex << ",";
					//fs << "response" << response;
				}
                	//response.release();
			} //end for simple filter
			
			if (mS1C1Suppress > 0){ //For this band/scale only
				cv::Mat s1res_work;
				suppressthresh = mins + mS1C1Suppress*(maxs-mins);
				int start = getGlobalFilterIndex(iBand,iScale,0);
				
				for (int i=start; i< start+mNumOrientations; i++) {
					s1Layer[i].copyTo(s1res_work);
					cv::compare(s1res_work,suppressthresh,mask,cv::CMP_GE); //Mask has values of 0 or 255
					cv::divide(mask,255.0,mask);
					mask.convertTo(mask,TYPE);
					cv::multiply(s1res_work,mask,s1Layer[i]); // Elements with values less than threshold get 0
					//s1res[i] = s1res_work.clone(); //TODO: Is this right?                         
				}                    
			}
    	} //end for scale
	} //end for band
	
	if (debug) cout << "s1Layer Size: " << s1Layer.size() << endl;
}

//=======================================================================================

//******************** C1 Layer ***********************

void CHMAX::makeC1Layer(const vector<cv::Mat>& s1Layer, vector<cv::Mat>& c1Layer)
{
	cv::FileStorage fs;
    	
	if (debug) {
		cout << "Starting C1 Layer" << endl; //cv::waitKey(-1);
		fs.open("c1Debug.xml", cv::FileStorage::WRITE); 
	}
	
	cv::Mat s1Mat, c1Mat, mmfMat, c1Mat2;
	cv::Mat c1MatMC;
	cv::Mat c1mins, c1mins2, c1maxs, c1maxs2;     
	
	if (jniDebug) {
		cout << "Enter makeC1Layer() ..." << endl;
		cout << "Num Scale Bands: " << mNumScaleBands << endl;
	}

	//Loop over scale bands
	for (int iBand=0; iBand<mNumScaleBands; iBand++){
		vector<cv::Mat> tempC1; // dfd: moved inside loop so that it is redefined at every band, so clear is not needed
		double temp = mC1SpaceSS.at<double>(0,iBand);
		cv::Scalar poolRange;
		poolRange[3] = (int) temp;
		
		//Loop over simple filters
		for (int iFilter=0; iFilter<mNumOrientations; iFilter++){ // dfd: numOrientations (12)
//			if (bandScales.size() > 0) bandScales.clear();
			vector<int> bandScales;  // dfd: define locally so you don't have to clear it
	     	getScalesInThisBand(iBand,bandScales);
						
			//Get first filter in scale band for c1Mat initialization               
			int iUFilter = getGlobalFilterIndex(iBand,0,iFilter);
			c1Mat = cv::Mat::zeros(s1Layer[iUFilter].rows,s1Layer[iUFilter].cols,TYPE);
			
			//Loop over band scales
			for(int iScale=0; iScale< (int)bandScales.size(); iScale++) {
				iUFilter = getGlobalFilterIndex(iBand,iScale,iFilter);
				s1Layer[iUFilter].copyTo(s1Mat); 
				cv::max(c1Mat,s1Mat,c1Mat2); // dfd: We are doing this operation twice. We did the same thing after calculating response. 
				c1Mat = c1Mat2.clone();
			} //end for scale
			
			if (jniDebug && iFilter == 0 & iBand == 0) {
				cout << "c1Mat: " << endl;
				for (int i=0; i<debugPts; i++) {
					cout << c1Mat.at<double>(0,i) << ",";
				}
				cout << endl;
			}

			if (jniDebug && iFilter == 0 && iBand == 0) cout << "Num Scales: " << bandScales.size() << endl;
			mymaxfilter(c1Mat,poolRange,mmfMat); // dfd: for each orientation, find scale and space max           
			
			if (iFilter == 0) {
				c1mins = mmfMat.clone(); 
				c1maxs = mmfMat.clone(); 
			} else {
				cv::min(c1mins,mmfMat,c1mins2); //Will be minimum across all filters
				c1mins = c1mins2.clone();
				cv::max(c1maxs,mmfMat,c1maxs2); //Will be maximum across all filters
				c1maxs = c1maxs2.clone();
			}

			if (jniDebug && iFilter == 0 && iBand == 0) {
				cout << "c1mins: " << endl;
				for (int i=0; i<debugPts; i++) {
					cout << c1mins.at<double>(0,i) << ",";
				}
				cout << endl;
				
				cout << "c1maxs: " << endl;
				for (int i=0; i<debugPts; i++) {
					cout << c1maxs.at<double>(0,i) << ",";
				}
				cout << endl;

				cout << "mmfMat: " << endl;
				for (int i=0; i<debugPts; i++) {
					cout << mmfMat.at<double>(0,i) << ",";
				}
				cout << endl;
			}
			
			tempC1.push_back(mmfMat.clone());    //To be used for suppression
			
		} //end for simple filters
		
		//////// Suppression (within band) /////////
		if (mS1C1Suppress > 0){
			cv::Mat tempC1_work;
			if (jniDebug) cout << "Suppression..." << endl;
			cv::Mat suppressthresh = cv::Mat::zeros(mmfMat.rows,mmfMat.cols,CV_64FC(mNumOrientations));
			suppressthresh = c1mins + mS1C1Suppress*(c1maxs-c1mins);
			
			for (int iFilter=0; iFilter<mNumOrientations; iFilter++){   
				cv::Mat mask;
				tempC1[iFilter].copyTo(tempC1_work);                 
				cv::compare(tempC1_work,suppressthresh,mask,cv::CMP_GE); //Mask has values of 0 or 255
				cv::divide(mask,255.0,mask); //Convert to 0 or 1
				mask.convertTo(mask,TYPE);
				cv::multiply(tempC1_work,mask,tempC1_work); // Elements with values less than threshold get 0
				tempC1[iFilter] = tempC1_work.clone();
			}
		} 
		
		cv::merge(tempC1,c1MatMC);
		c1Layer.push_back(c1MatMC);
        //		tempC1.clear(); //Should be empty for processing the next band. dfd: clear not needed
	} //end for band
	
	if (jniDebug) cout << "c1Layer Size: " << c1Layer.size() << endl;
    
	fs.release();
} //end makeC1Layer

//=======================================================================================

void CHMAX::makeS2C2Features(vector<cv::Mat>& c1Layer, cv::Mat& c2Features) 
{
	// Inputs:
	//c1Layer = is a std vector of 3D cv::Mats
	//patchSize=Indicates the patch sizes
	cout << "Entering makeS2C2Features() ..." << endl;

	int numBands = c1Layer.size();
	int numOrientations = c1Layer[0].channels(); //nfiltes need to be identical to s2Target[0].size.p[2]; 
	
	if (mPatchesVec[0][0].depth() != c1Layer[0].depth()) {
        cout << "The types of c1Layer and mPatchesVec must be the same." << endl;//make this an error later. dfd: exit
        exit(1);
    }
	if (debug) cout << "numBands: " << numBands << endl;
	if (debug) cout << "numOrientations: " << numOrientations << endl;
	
	int totalNumPatches = 0;
	for (int i=0; i<mNumPatchSizes; i++) {
        int numPatchesPerSize = mPatchesVec[i].size();
        totalNumPatches += numPatchesPerSize; // dfd: We should just read mNumPatches instead, or assert the result should be the same
	}
    
	if (jniDebug) {
		cout << "Num Bands: " << numBands << endl;
		cout << "Num Orientations: " << numOrientations << endl;
		cout << "Num Patch Sizes: " << mNumPatchSizes << endl;
		cout << "Total Num Patches: " << totalNumPatches << endl;
	}

	cv::Mat tempC2res(1,totalNumPatches,TYPE);
	cv::Mat tempcols;
	
	int curind, lastind = 0; //For forming one feature vector
	
    // dfd: we use mPatchesVec dimensions instead of s2Target dimensions because s2Target exists only if we read patches from xml file, 
    // mS2Target is not defined when we create patches from calib images
	for (int iPatchSize = 0; iPatchSize < mNumPatchSizes; iPatchSize++) { 
		int numPatchesPerSize = mPatchesVec[iPatchSize].size();
		cv::Mat c2Res(1,numPatchesPerSize,TYPE,0.0);
		if (jniDebug) cout << "Num Patches Per Size: " << numPatchesPerSize << endl;
		for (int iCenter = 0; iCenter < numPatchesPerSize; iCenter++) {
			//extract serialized patch from s2Target and reshape:
			vector<cv::Mat> patch;
			cv::split(mPatchesVec[iPatchSize][iCenter], patch);

			for (int iBand= 0; iBand < numBands; iBand++) {
				cv::Mat dst_vec_out;
				double minVal, maxVal;	
				cv::Point minLoc,maxLoc;
				vector<cv::Mat> c1Vec;
				cv::split(c1Layer[iBand],c1Vec);
				if (jniDebug && iBand == 0 && iPatchSize == 0 && iCenter == 0) {
				WindowedPatchDistance(c1Vec, patch, dst_vec_out, mDoSkipSquareRoot);
				} else {
				WindowedPatchDistance(c1Vec, patch, dst_vec_out, mDoSkipSquareRoot);
				}

				if (jniDebug && iBand == 0 && iPatchSize == 0 && iCenter == 0) {
					cout << "c1Vec[0]: " << endl;
					for (int i=0; i<debugPts; i++) {
						cout << c1Vec[0].at<double>(0,i) << ",";
					}				
					cout << endl;

					cout << "Patch[0]: " << endl;
					for (int i=0; i<debugPts; i++) {
						cout << patch[0].at<double>(0,i) << ",";
					}				
					cout << endl;

					cout << "Windowed Patch Distance Vals: " << endl;
					for (int i=0; i<debugPts; i++) {
						cout << dst_vec_out.at<double>(0,i) << ",";
					}				
					cout << endl;
				}

				cv::minMaxLoc(dst_vec_out, &minVal, &maxVal, &minLoc, &maxLoc);
			
                if (iBand == 0){
					c2Res.at<double>(0,iCenter) = minVal;
				}
				
				if (minVal < c2Res.at<double>(0,iCenter)) {
					c2Res.at<double>(0,iCenter) = minVal;
				}
			} //end for iBand
		} // end for iCenter
		curind = lastind;
		lastind += numPatchesPerSize;
		tempcols = tempC2res.colRange(curind,lastind); //The assignment is just a Mat header, not the column data.
		if (!mDoSkipSquareRoot) {
			c2Res.copyTo(tempcols);	
		}
		else {
			exp(-c2Res, tempcols); // dfd: use exp(-min squared distance) as feature component
		}
    } // for iPatchSize
	
	tempC2res.copyTo(c2Features);

	if (jniDebug) {
		cout << "c2 Features: " << endl;
		for (int i=0; i<debugPts; i++) {
			cout << c2Features.at<double>(0,i) << ",";
		}
		cout << endl;
	}
	cout << "Exiting makeS2C2Features()" << endl;
} //end makeS2C2Features

//=======================================================================================

void CHMAX::makePatchesVec(const vector<cv::Mat>& imgs, vector<vector<cv::Mat> >& patchesVec)
// This method runs through every calibration image and makes a C1 layer for each image
// then it collects a number of patches from each image C1 and adds them to the vector of patches patchesVec
{
    MyTimer patchMakingTimer;
    patchMakingTimer.start();
    
    patchesVec.resize(mNumPatchSizes); //dfd: initialize bins for patches of given sizes
    // runs through the C1 layer of the    
    srand((unsigned)31415927); // initialize seeds. Use time() if need for time varying randomness
	vector<int> randBandIndices(mNumPatches);
    vector<int> randPatchSizeIndices(mNumPatches);
    
	int numImages = imgs.size();
	int numPatchesPerImage;
	int totalPatchCount = 0;
	int imageLimit = numImages<mNumPatches? numImages:mNumPatches; // dfd: if mNumPatches less than numImages, don't process all the images
	int patchLimit = floor(mNumPatches / numImages); // dfd: zero if requested mNumPatches is less than numImages
	int extraPatches = mNumPatches % numImages; // dfd: remainder of mNumPatches/numImages. mNumPatches if mNumPatches<numImages
	for (int imageIndex=0; imageIndex<imageLimit; imageIndex++){
#if 0
		procS1C1(imgs[imageIndex]);
		vector<cv::Mat> c1Layer = getC1Res(); //dfd: I don't like this.
#else
		vector<cv::Mat> s1Layer, c1Layer;
		makeS1Layer(imgs[imageIndex], s1Layer);
		makeC1Layer(s1Layer, c1Layer);
#endif
		if (imageIndex < extraPatches) {
			numPatchesPerImage = patchLimit + 1;
		}
		else {
			numPatchesPerImage = patchLimit;
		}
		for (int patchCount=0; patchCount<numPatchesPerImage; patchCount++) {
            int randBandIndex, randPatchSize, randPatchSizeIndex, numRows, numCols;
            while(1){
                randBandIndex = rand() % mNumScaleBands; // dfd: ( value % 100 ) is in the range 0 to 99
                randPatchSizeIndex = rand() % mNumPatchSizes;
                randPatchSize = mPatchSizes.at<double>(0,randPatchSizeIndex);
                
                numRows = c1Layer[randBandIndex].rows;
                numCols = c1Layer[randBandIndex].cols;
                int reducedNumRows = numRows-randPatchSize; // dfd: range of rows and cols that top left corner of patch can belong to. Can be 0, but not negative
                int reducedNumCols = numCols-randPatchSize;
                if ((reducedNumRows>=0) && (reducedNumCols>=0)) { // dfd: then patch can fit in image
                    break;
                }
                else {
                    cout << "Patch is larger than image. Throw dice again..." << endl;
                }
                
            }
            vector<cv::Mat> c1Vec;
            cv::split(c1Layer[randBandIndex],c1Vec); // find the c1 images for different orientations
			int randRow = rand() % (numRows-randPatchSize+1); // between 0 and numRows-randPatchSize included
			int randCol = rand() % (numCols-randPatchSize+1); // between 0 and numRows-randPatchSize included
            cv::Range range1(randRow, randRow+randPatchSize); // right bound is EXCLUSIVE
			cv::Range range2(randCol, randCol+randPatchSize);
			cv::Range ranges[] = {range1, range2};
			
			vector<cv::Mat> patch;
			for (int iOrient=0; iOrient<mNumOrientations; iOrient++){
                //				cv::Mat c1ForOneOrient = c1Vec[iOrient];
				patch.push_back(c1Vec[iOrient](ranges)); // dfd: grab a patch in the right range in each orientation
			}
			cv::Mat mergedPatch;
			cv::merge(patch, mergedPatch); // dfd: merge all orientations of patch vector into channels of single mergedPatch
			patchesVec[randPatchSizeIndex].push_back(mergedPatch);
            totalPatchCount++;
            cout << "Added patch number " << totalPatchCount << endl;
		}
	}	
    patchMakingTimer.stop();
    cout << "Finding " << mNumPatches << " patches from " << numImages << " images took " << patchMakingTimer.time() << " seconds" << endl;
}

//=======================================================================================


//==============================================================================
//                 Class helper functions
//==============================================================================
//METHOD: Returns integer array of scales for specified band
void CHMAX::getScalesInThisBand(int iBand,vector<int>& bandScales)
{
	//Matlab equivalent: mC1ScaleSS(iBand):mC1ScaleSS(iBand+1)-1
	int bandStart = (int) mC1ScaleSS.at<double>(0,iBand); 
	int bandStop = (int) (mC1ScaleSS.at<double>(0,iBand+1)-1); 
	
	for (int i=bandStart; i<bandStop+1; i++){
		bandScales.push_back(i);
	}
}

//========================================================================================================

//METHOD: Returns simple filter response for given image and filter
void CHMAX::getSimpleFilterResponse(const cv::Mat& inIm, const cv::Mat& s1Norm, int iUFilterIndex, cv::Mat& sfilter, bool REMOVEBORDERS)
{
	cv::Mat sqfilter;
	cv::Mat img = inIm.clone();
	
	// Retrieve Gabor Filter at desired index (filter is square kernel)     
	mGaborFilters[iUFilterIndex].copyTo(sqfilter);
	
	cv::FileStorage fs2;
	if (debug && iUFilterIndex == 12) fs2.open("filterDebug.xml", cv::FileStorage::WRITE); 
	if (debug && iUFilterIndex == 12) fs2 << "sqfilter_init" << sqfilter.clone();
	
	// Apply filter
	conv2(&img,&sqfilter,&sfilter);
	sfilter = abs(sfilter); // dfd: Check why this is necessary
	
	if (debug && iUFilterIndex == 12) fs2 << "sfilter_conv2" << sfilter.clone();
	
	cv::divide(sfilter,s1Norm,sfilter); 
	if (debug && iUFilterIndex == 12) fs2 << "s1Norm" << s1Norm.clone();
	if (debug && iUFilterIndex == 12) fs2 << "output" << sfilter.clone();
	//TODO: if (REMOVEBORDERS) {removeBorders(sfilter,filtSizes[iUFilterIndex],sfilter);}
	
	sqfilter.release();
	img.release();
	fs2.release();
}

//========================================================================================================

void CHMAX::findFilterResponse(const cv::Mat& inIm, const cv::Mat& s1Norm, int filterIndex, cv::Mat& absFilteredIm)
// dfd: Simplified version of getSimpleFilterResponse. 
// Unfortunately, response requires cloning when it is pushed into s1Layer!
// This function is not currently used.
{
    cv::Mat filteredIm, normalizedFilteredIm;
    
    int borderMode = 0;  // dfd: this is BORDER_CONSTANT, which set the border to zero. This is not BORDER_DEFAULT
    int sourceDepth = -1; // dfd: use source depth
    cv::filter2D(inIm, filteredIm, sourceDepth, mGaborFilters[filterIndex], cv::Point(-1,-1), 0, borderMode);  
    
	cv::divide(filteredIm, s1Norm, normalizedFilteredIm);
    absFilteredIm = abs(normalizedFilteredIm); // dfd: Eq. 2 in Mutch & Lowe
}

//========================================================================================================

//METHOD: Returns the global filter index for given band, scale, and filter
int CHMAX::getGlobalFilterIndex(int iBand,int iScale,int iFilt){
	// This function assumes the original s1res Mats were stored via 3 loops:
	// for iBand...; for iScale...; for iFilt...  Mixing this order will return undesired results.
	vector<int> bandScales;
	getScalesInThisBand(iBand,bandScales);
	
	int S,F;
	S = bandScales.size();  
	F = mNumOrientations;
	
	return ((iBand)*S*F + (iScale)*F + (iFilt));
} 

//========================================================================================================

//METHOD: Extracts Gabor filters into vector. Takes matrix of filters as input.
void CHMAX::extractFilters(const cv::Mat& gaborFilterMat, const cv::Mat& filterSizeMat)
{
	// dfd: There is only one filter per column in gaborFilterMat, so there are many zeros in columns of small filters	
	mNumFilters = gaborFilterMat.cols;
	mNumOrientations = mNumFilters / mNumScales; // dfd: also called numSimpleFilters, which is less informative
	
	if (debug) cout << "numFilters before: " << mNumFilters << endl;
	
	// Extract filter sizes into integer array
	mFilterSizes.clear();
	mFilterSizes.resize(mNumFilters); // dfd: repeated sizes: 7, 7, 7, ...., 9, 9, ...., 39, 39, ..., 39
	for (int iFilter=0; iFilter<mNumFilters; iFilter++){ 
		mFilterSizes[iFilter] = (int) filterSizeMat.at<double>(0,iFilter);
	} 
	
	// Extract single column of filters and reshape into a NxN matrix
	cv::Mat currentFilter,tempFilter;
	for (int iFilter=0; iFilter<mNumFilters; iFilter++){
		cv::Mat temp = gaborFilterMat.col(iFilter); // Only a header is returned
		currentFilter = temp.clone();
		int filterSize = mFilterSizes[iFilter];  
		int numElements = filterSize * filterSize; // dfd: was (int) pow((double)numel,2.0);
		tempFilter = currentFilter.rowRange(0,numElements);
		currentFilter = tempFilter.clone(); // dfd: Gabor filters are not 2D arrays at this point
		currentFilter = currentFilter.reshape(1, filterSize).t();// Transpose because Matlab is col-major and OpenCV is row-major
		cv::flip(currentFilter,currentFilter,-1); // Flip up/down and left/right //TODO: Is this right?
		mGaborFilters.push_back(currentFilter); 
	}
} 

//========================================================================================================

// METHOD: Extracts single patch from s2Target matrix
void CHMAX::extractPatch(const int iS2, const int iCenter, cv::Mat& patch) {
	//Assume patchSize has [3] and [2] already filled with values
	
	vector<cv::Mat> patches;
	int rows = (int) mPatchSizes.at<double>(0,iS2); // dfd: mPatchSizes has one column per patch, with first row showing nb of rows of patch and second row showing nb of cols
	int columns = (int) mPatchSizes.at<double>(1,iS2); 
	
	cv::Mat s2Target0;
	mS2Target[iS2].copyTo(s2Target0);
	cv::Mat s2Tgt = s2Target0.col(iCenter).clone();
	s2Tgt.convertTo(s2Tgt,TYPE);
	
	for (int i=0; i<mNumOrientations; i++){
		cv::Mat temp = s2Tgt.rowRange(i*(rows*columns),(i+1)*(rows*columns)).clone();
		//cout << temp.rows << "," << temp.cols << endl;
		temp = temp.reshape(1,rows).t(); //NOTE: transpose because opencv is row-major. 1 channel, nb of rows = rows. nb of cols implicit
		patches.push_back(temp);
		//		temp.release();
	}
	
	cv::merge(patches, patch); // make a 3D matrix which is split again when it is used
	patches.clear();
}

//========================================================================================================

//METHOD: Extracts all patches into single vector for a given s2Target matrix
void CHMAX::extractAllPatches() 
{
	mPatchesVec.clear();
	
	//Extract patches for multiple s2Target matrices
	for (int iS2=0; iS2 < mS2Target.size(); iS2++) {
		int numRBFCenters = mS2Target[iS2].cols;
		vector<cv::Mat> tempMatVec;
		for (int iCenter=0; iCenter < numRBFCenters; iCenter++) {
			cv::Mat temp;
			extractPatch(iS2,iCenter,temp); 
			tempMatVec.push_back(temp);
			temp.release();
		}
		mPatchesVec.push_back(tempMatVec);
	}	
}

//========================================================================================================

//========= (RE)INITIALIZE HMAX OBJECT ==========
//METHOD: (Re)Initializes HMAX object parameters from XML file
void CHMAX::initialize(const string XMLfile) {
	cout << "Initializing HMAX" << endl;
	mDoSkipSquareRoot = 0;
	cout << "doSkipSquareRoot: " << mDoSkipSquareRoot << endl;

	mC1ScaleSS.release();
	mC1SpaceSS.release();
	mPatchSizes.release();
	mS2Target.clear();
	mGaborFilters.clear();
	mFilterSizes.clear();
	
	cv::FileStorage read_fs_mat(XMLfile, cv::FileStorage::READ);
	cv::Mat gaborFilterMat, filterSizeMat, fSizes,c1OL,cPatches,numPatchSizes, s1C1Suppress;
	read_fs_mat["filters"] >> gaborFilterMat;
	cout << "Num gabor filters: " << gaborFilterMat.cols << endl;
	read_fs_mat["fSiz"] >> filterSizeMat;
	read_fs_mat["c1ScaleSS"] >> mC1ScaleSS;
	cout << "c1ScaleSS: " << mC1ScaleSS << endl;
	read_fs_mat["c1SpaceSS"] >> mC1SpaceSS;
	cout << "c1SpaceSS: " << mC1SpaceSS << endl;
	read_fs_mat["c1OL"] >> c1OL;
	mC1Overlap = (int) c1OL.at<double>(0,0);
	read_fs_mat["suppressionparams"] >> s1C1Suppress;
	mS1C1Suppress = s1C1Suppress.at<double>(0,0);
	cout << "Suppression threshold: " << mS1C1Suppress << endl;
	read_fs_mat["patchSizes"] >> mPatchSizes;
	cout << "Patch sizes: " << mPatchSizes << endl;
	read_fs_mat["num_patches"] >> numPatchSizes;
	
	mNumPatchSizes = (int) numPatchSizes.at<double>(0,0);
	cout << "Num Patch Sizes: " << numPatchSizes << endl;
	mNumScaleBands = mC1ScaleSS.cols-1;
	mNumScales = mC1ScaleSS.at<double>(0,mNumScaleBands)-1;
    string patchbase = "cPatches_";
	char numstr[50];
	mNumPatches = 0;
	for (int i=0; i<mNumPatchSizes; i++) { 
		sprintf(numstr,"%d",i+1);
		cout << (patchbase + numstr) << endl;
		read_fs_mat[patchbase + numstr] >> cPatches;
		mNumPatches += cPatches.cols; // compute total number of patches as sum of numbers of each patch size
		mS2Target.push_back(cPatches);
	}
	read_fs_mat.release();
	
	cout << "Extracting gabor filters" << endl;
	extractFilters(gaborFilterMat, filterSizeMat); // mGaborFilters, mFilterSizes and mNumOrientations are defined there
	cout << "Finished extracting gabor filters" << endl;
	cout << "Extracting patches" << endl;
	extractAllPatches(); // dfd: mPatchVec is made here
	cout << "Finishing extracting patches" << endl;
	cout << "Exiting initialize()" << endl;
}

//========================================================================================================

void CHMAX::initializeWithoutPatches(const string xmlFileName,  int numPatches) 
// dfd: Initialize Gabor filters and sizes, but don't read patches
{
	cout << "Initializing HMAX without reading patches" << endl;
	mDoSkipSquareRoot = 0;
	cout << "doSkipSquareRoot: " << mDoSkipSquareRoot << endl;

	mC1ScaleSS.release();
	mC1SpaceSS.release();
	mPatchSizes.release();
	mGaborFilters.clear();
	mFilterSizes.clear();
	
	mNumPatches = numPatches;
	cv::FileStorage read_fs_mat(xmlFileName, cv::FileStorage::READ);
	cv::Mat gaborFilterMat, filterSizeMat, fSizes,c1OL,cPatches,numPatchSizes, s1C1Suppress;
	read_fs_mat["filters"] >> gaborFilterMat;
	read_fs_mat["fSiz"] >> filterSizeMat;
	read_fs_mat["c1ScaleSS"] >> mC1ScaleSS;
	read_fs_mat["c1SpaceSS"] >> mC1SpaceSS;
	read_fs_mat["c1OL"] >> c1OL;
	mC1Overlap = (int) c1OL.at<double>(0,0);
	read_fs_mat["suppressionparams"] >> s1C1Suppress;
	mS1C1Suppress = s1C1Suppress.at<double>(0,0);
	read_fs_mat["patchSizes"] >> mPatchSizes;
	read_fs_mat["num_patches"] >> numPatchSizes;
	
	mNumPatchSizes = (int) numPatchSizes.at<double>(0,0);
	mNumScaleBands = mC1ScaleSS.cols-1;
	mNumScales = mC1ScaleSS.at<double>(0,mNumScaleBands)-1;
    
	extractFilters(gaborFilterMat, filterSizeMat); // mGaborFilters, mFilterSizes and mNumOrientations are defined there
	
	read_fs_mat.release();
}

//=======================================================================================
//=======================================================================================
