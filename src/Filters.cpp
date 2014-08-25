/*
 *  Filters.cpp
 *  Applied various of the filters (e.g. convolution, sumfilter)
 *  needed by HMAX classifier
 *
 *  Created by Rodriguez, Pedro A., March 2012 -- July 2012
 *  Copyright 2012 Pedro A. Rodriguez and JHU/APL. All rights reserved.
 *  Modified by Daniel DeMenthon (dfd), Sept. 2012
 */


using namespace std;

#include <iostream>
#include <math.h>
#include <string.h>
#if 0
#include <sys/types.h>
#include <sys/param.h>
#include <sys/stat.h>
#include <sys/time.h>
#endif

#include "Filters.h"

static const int TYPE = CV_64F;

//========================================================================================================

/* subroutine to do 2D convolution like in Matlab, with the "same" option */
void conv2(cv::Mat* src, cv::Mat* kernel, cv::Mat* dst) {
    
    CvPoint anchor = cvPoint(kernel->cols - kernel->cols/2 - 1, kernel->rows - kernel->rows/2 - 1);
    
    cv::Mat rot_kernel2;
    cv::flip(*kernel, rot_kernel2, -1); // dfd: DO NOT FLIP. filter2D uses A CORRELATION, not a convolution
    
    int borderMode = 0;
    cv::Mat dst2;
    cv::filter2D(*src, dst2, src->depth(), rot_kernel2, anchor, 0, borderMode);
    *dst=dst2;
    
    //dst2.release();
#if 0
    rot_kernel2.release();
#endif
}

//========================================================================================================


/* subroutine to imitate sumfilter.m function in Matlab */
void sumfilter(cv::Mat* I, cv::Scalar radius, cv::Mat* dst) {
    
	// I is the input image
	// radius is the additional radius of the window, i.e., 5 means 11 x 11
	// if a four value vector is specified for radius, then any rectangular support may be used for max.
	// in the order left top right bottom.
    
	// if (size(I,3) > 1)
	//     error('Only single-channel images are allowed');
	// end
	//
    cv::Mat kernel;
    cv::Mat dst2;
	
    if (radius[0] == 0 && radius[1] == 0 && radius[2] == 0){
		//Assumes that only radius[3] has a value:
		
		// Matlab Equivalent:
		// mask = ones(2*radius+1,1);
		// I3_pedro = conv2(I,mask*mask','same');
		
		cv::Mat mask(2*radius[3]+1, 1, I->type(), 1.0);
		cv:: mulTransposed(mask, kernel, false);
		//kernel=mask*mask.t();
		
		mask.release();
		//cout << "verify radius values" << endl;
		
		
    } else {
		
		//  Matlab Equivalent:  I3_pedro=conv2(I,ones(radius(1)+radius(3)+1,1)*ones(1,radius(2)+radius(4)+1),'same');
        
        cv::Mat vect_1(radius[0]+radius[2]+1, 1, I->type(), 1.0);
        cv::Mat vect_2(1, radius[1]+radius[3]+1, I->type(), 1.0);
        kernel = vect_1*vect_2;
        
        cv::Size sz_kernel;
        sz_kernel= kernel.size();
		//cout << "Kernel First dimension " << sz_kernel.width << endl;
		//cout << "Kernel Second dimension " << sz_kernel.height << endl;
        
#if 0
        vect_1.release();
        vect_2.release();
#endif
        
    }
	//Do actual convolution
	conv2(I, &kernel, &dst2);
	*dst=dst2.clone();
	
	kernel.release();
	//dst2.release();
}

void WindowedPatchDistance(cv::Mat& Im, cv::Mat& Patch, cv::Mat& dst) {
	//void WindowedPatchDistance(cv::MatND* Im, cv::MatND* Patch) {
	
	//Estimate dimensions of Im and Patch:
	int numDim=Im.dims;
	int length = Im.size.p[0];
	int width = Im.size.p[1];
	int depth=Im.size.p[2];
	int matriz_sz = length * width;
	
	cout << "numDim: " << numDim << endl;
	cout << "depth: " << depth << endl;
	cout << "matriz_sz: " << matriz_sz << endl;
	
	//   int numDim_Patch=Patch.dims;
	int length_Patch = Patch.size.p[0];
	int width_Patch = Patch.size.p[1];
	int depth_Patch=Patch.size.p[2];
	// int matriz_sz_Patch = length_Patch * width_Patch;
	
    if (depth != depth_Patch) {
        cout << "The third dimension of IM and Patch must be the same." << endl;//make this an error later
    }
	
    //Initialize variables:
    cv::Scalar Psqr;
    Psqr[0]=0;
	
    //Maybe check to make sure they are both of the same type.
    int type_data=Im.type();
    cout << "IM type:  " << type_data << endl;
	
	//Take 3D IM and Patch and save them into Awork and Bwork and the rest is the same.
	
    cv::Mat Awork_cv(length, width, type_data, 0);
    cv::Mat Awork_square(length, width, type_data, 0);
	
    cv::Range range1(0,length);
    cv::Range range2(0,width);
	
    cv::Mat Imsq(length,width, type_data, 0);
    cv::Mat Imsq_temp(length,width, type_data, 0);
    cv::Mat conv_temp(length,width, type_data, 0);
	
    cv::Mat PI(length,width, type_data, 0);
    cv::Mat PI_temp(length,width, type_data, 0);
	
    cv::Mat Bwork_cv(length_Patch, width_Patch,type_data, 0);
    cv::Mat Bwork_square(length_Patch, width_Patch,type_data, 0);
    
    cv::Range range1_Patch(0,length_Patch);
    cv::Range range2_Patch(0,width_Patch);//Reverse direction
	
	
    for (int i_chip = 0; i_chip < depth; i_chip++) {
		
		//Read slices of IM:
		cv::Range range3(i_chip,i_chip+1);
		cv::Range ranges[] = {range1,range2,range3};
		Im(ranges).copyTo(Awork_cv);//OR Aword_cv=Im(ranges)
		
		
		//Awork_cv.create(length,width,type_data);
		cv::Mat miniCubeA;
		miniCubeA.create(length,width,type_data);
		cv::resize(Awork_cv,miniCubeA,miniCubeA.size(),0,0);
		Awork_cv=miniCubeA;
		miniCubeA.release();
		
        //cout << "dimension Awork_cv:  " << Awork_cv.dims << endl;
        //cout << "element 2,3: " << Awork_cv.at<double>(2,3) << endl;
		
		//Read slices of Patch:
		cv::Range range3_Patch(i_chip,i_chip+1);
		cv::Range ranges_Patch[] = {range1_Patch,range2_Patch,range3_Patch};
		Patch(ranges_Patch).copyTo(Bwork_cv);//OR Aword_cv=Im(ranges)
		
		//Bwork_cv.create(length_Patch,width_Patch,type_data);
		cv::Mat miniCubeB;
		miniCubeB.create(length_Patch,width_Patch,type_data);
		cv::resize(Bwork_cv,miniCubeB,miniCubeB.size(),0,0);
		Bwork_cv=miniCubeB;
		miniCubeB.release();
		
		
		cv::flip(Bwork_cv, Bwork_cv, -1);  // flip prior to convolution
		
		//cout << "dimension Bwork_cv:  " << Bwork_cv.dims << endl;
		
        //Calculate PSQR
        cv::multiply(Bwork_cv, Bwork_cv, Bwork_square);
        Psqr+= cv::sum(Bwork_square);
		
        //cv::multiply(Awork_cv, Awork_cv, Awork_square);
        Awork_square=Awork_cv.mul(Awork_cv);
        if (i_chip == 0){
            //Imsq=Awork_square.clone();
            Awork_square.copyTo(Imsq);
        } else {
            Imsq_temp=Imsq;
            cv::add(Imsq_temp, Awork_square, Imsq);
        }
		
		
        //Create PI (conv2 stuff) 
        conv2(&Awork_cv, &Bwork_cv, &conv_temp);
        
        if (i_chip == 0){
            PI=conv_temp;
        } else {
            PI_temp=PI;
            cv::add(PI_temp, conv_temp, PI);
        }        
		
    } // end for loop in 3rd dimension
	
    cv::Scalar sum_support;
    sum_support[0]=(int)((double)length_Patch/2.0+0.5)-1;
    sum_support[1]=(int)((double)width_Patch/2.0+0.5)-1;
    sum_support[2]=(int)((double)length_Patch/2.0);
    sum_support[3]=(int)((double)width_Patch/2.0);
	
    //if prune==0
    cv::Mat Imsq_filter;
    sumfilter(&Imsq, sum_support, &Imsq_filter);
    Imsq=Imsq_filter;
	
    //else
    //    II = zeros(size(Im, 1), size(Im, 2));
    //PatchOrientations = +(Patch~=0);
    //for i = 1:dIm
    //        II = II + conv2(Imsq(:, :, i), PatchOrientations(:, :, i), 'same');
    //        end
    //                Imsq = II;
    //end
	
    cv::Mat D_temp(length, width, type_data, 0);
    cv::Mat D(length, width, type_data, 0);
    
    //Calculating:  Imsq - 2 * PI + Psqr
    D=Imsq-2*PI;
    D+=Psqr;
	
	cv::Mat tempmask;
	double thresh = 1e-14;
	cv::compare(D,thresh,tempmask,cv::CMP_GE);
	cv::divide(tempmask,255.0,tempmask);
	tempmask.convertTo(tempmask,TYPE);
	cv::multiply(D,tempmask,D);
	
    sqrt(D, D_temp);
    // D(abs(D) < 10e-15) = 0;
    D=D_temp;
	
    //cout << "dimension D:  " << D.dims << endl;
	
    
	//     //double thresh = 10e-10;
	//     cv::Mat A = cv::abs(D) > thresh;    //(elements we want to keep are one, zeros otherwise)
	//     A.convertTo(A , CV_64FC1);          //(the > operator gives an 8-bit array, I believe
	//                                         //we have to convert it back to the data type of array D
	//                                         //for the next line to work)
	//     D = A.mul(D);
	
    dst=D; //OUTPUT
	
    //Release mat data
	//    PI.release();
	//    Imsq.release();
	//    
	//    Awork_cv.release();
	//    Bwork_cv.release();
	//    
	//    Awork_square.release();
	//    Bwork_square.release();
	//    
	//    Imsq_temp.release();
	//    PI_temp.release();
	//    conv_temp.release();
	//    Imsq_filter.release();
	
    //Should I release this?:
	// D_temp.release();
	// D.release();
	
}

void WindowedPatchDistance(vector<cv::Mat>& Im, vector<cv::Mat>& Patch, cv::Mat& dst) {//using vector version
	
	//Estimate dimensions of Im and Patch:
	// int numDim=Im[0].dims;
	int length = Im[0].rows; //Assuming all matrices have the same size
	int width = Im[0].cols;
	int depth=Im.size();
	// int matriz_sz = length * width;
	
	//cout << "numDim: " << numDim << endl;
	//cout << "depth: " << depth << endl;
	//cout << "matriz_sz: " << matriz_sz << endl;
	
	//int numDim_Patch=Patch[0].dims;
	int length_Patch = Patch[0].size.p[0];
	int width_Patch = Patch[0].size.p[1];
	int depth_Patch=Patch.size();
	//int matriz_sz_Patch = length_Patch * width_Patch;
	
    if (depth != depth_Patch) {
        cout << "The third dimension of IM and Patch must be the same." << endl;//make this an error later
        cv::waitKey(-1);
    }
	
	//Initialize variables:
    cv::Scalar Psqr;
    Psqr[0]=0;
	
    //Maybe check to make sure they are both of the same type.
    int type_data=Im[0].type();
    //cout << "IM type:  " << type_data << endl;
    //cout << "IM size: " << depth << endl;
    //cout << "Patch size: " << depth_Patch << endl;
	
	//Take 3D IM and Patch and save them into Awork and Bwork and the rest is the same.
	
    cv::Mat Awork_cv(length, width, type_data, 0);
    cv::Mat Awork_square(length, width, type_data, 0);
	
    cv::Mat Imsq(length,width, type_data, 0);
    cv::Mat Imsq_temp(length,width, type_data, 0);
    cv::Mat conv_temp(length,width, type_data, 0);
	
    cv::Mat PI(length,width, type_data, 0);
    cv::Mat PI_temp(length,width, type_data, 0);
	
    cv::Mat Bwork_cv(length_Patch, width_Patch,type_data, 0);
    cv::Mat Bwork_square(length_Patch, width_Patch,type_data, 0);
    
    for (int i_chip = 0; i_chip < depth; i_chip++) {
		
		//Read slices of IM:
		Im[i_chip].copyTo(Awork_cv);
		
        //cout << "dimension Awork_cv:  " << Awork_cv.dims << endl;
        //cout << "element 2,3: " << Awork_cv.at<double>(2,3) << endl;
		
		//Read slices of Patch:
		Patch[i_chip].copyTo(Bwork_cv);
		cv::flip(Bwork_cv, Bwork_cv, -1);  // flip prior to convolution
		
		//cout << "dimension Bwork_cv:  " << Bwork_cv.dims << endl;
		
        //Calculate PSQR
        cv::multiply(Bwork_cv, Bwork_cv, Bwork_square);
        Psqr+= cv::sum(Bwork_square);
		
        //cv::multiply(Awork_cv, Awork_cv, Awork_square);
        Awork_square=Awork_cv.mul(Awork_cv);
        if (i_chip == 0){
            //Imsq=Awork_square.clone();
            Awork_square.copyTo(Imsq);
        } else {
            Imsq_temp=Imsq;
            cv::add(Imsq_temp, Awork_square, Imsq);
        }
		
		
        //Create PI (conv2 stuff) 
        conv2(&Awork_cv, &Bwork_cv, &conv_temp);
        
        if (i_chip == 0){
            PI=conv_temp;
        } else {
            PI_temp=PI;
            cv::add(PI_temp, conv_temp, PI);
        }        
		
    } // end for loop in 3rd dimension
	
	cv::Scalar sum_support;
    sum_support[0]=(int)((double)length_Patch/2.0+0.5)-1;
    sum_support[1]=(int)((double)width_Patch/2.0+0.5)-1;
    sum_support[2]=(int)((double)length_Patch/2.0);
    sum_support[3]=(int)((double)width_Patch/2.0);
	
    //if prune==0
    cv::Mat Imsq_filter;
    sumfilter(&Imsq, sum_support, &Imsq_filter);
    Imsq=Imsq_filter;
	
	
	
    //else
    //    II = zeros(size(Im, 1), size(Im, 2));
    //PatchOrientations = +(Patch~=0);
    //for i = 1:dIm
    //        II = II + conv2(Imsq(:, :, i), PatchOrientations(:, :, i), 'same');
    //        end
    //                Imsq = II;
    //end
	
    cv::Mat D_temp(length, width, type_data, 0);
    cv::Mat D(length, width, type_data, 0);
    
    //Calculating:  Imsq - 2 * PI + Psqr
    D=Imsq-2*PI;
    D+=Psqr;
	
	cv::Mat tempmask;
	double thresh = 1e-14;
	cv::compare(D,thresh,tempmask,cv::CMP_GE);
	cv::divide(tempmask,255.0,tempmask);
	tempmask.convertTo(tempmask,TYPE);
	cv::multiply(D,tempmask,D);
	
    sqrt(D, D_temp);
    // D(abs(D) < 10e-15) = 0;
    D=D_temp;
	
    //cout << "dimension D:  " << D.dims << endl;
	
    
	//     //double thresh = 10e-10;
	//     cv::Mat A = cv::abs(D) > thresh;    //(elements we want to keep are one, zeros otherwise)
	//     A.convertTo(A , CV_64FC1);          //(the > operator gives an 8-bit array, I believe
	//                                         //we have to convert it back to the data type of array D
	//                                         //for the next line to work)
	//     D = A.mul(D);
	
    dst=D; //OUTPUT
	
    //Release mat data
	//    PI.release();
	//    Imsq.release();
	//    
	//    Awork_cv.release();
	//    Bwork_cv.release();
	//    
	//    Awork_square.release();
	//    Bwork_square.release();
	//    
	//    Imsq_temp.release();
	//    PI_temp.release();
	//    conv_temp.release();
	//    Imsq_filter.release();
	
    //Should I release this?:
	// D_temp.release();
	// D.release();
	
	
}


void mymaxfilter(cv::Mat& Im, cv::Scalar poolrange, cv::Mat& dst) {
	
    //int type_data=Im.type();
	
    int numrows = Im.rows;
    int numcols = Im.cols;
	
    double dnumrows=(double)Im.rows;
    double dnumcols=(double)Im.cols;
    
    int fullpool=poolrange[3];
    double halfpool = poolrange[3]/2.0;
    
    int length=(int)ceil(dnumrows/(halfpool));
    int width =(int)ceil(dnumcols/(halfpool));
	
    int therowindices[length];
    int thecolindices[width];
	
    //cout << "rowindiceslength: " << length << endl;
    //cout << "colindiceslength: " << width << endl;
	
	
    double i_idx=0;
    for (int ct = 0; ct < length; ct++) {
		therowindices[ct] =(int)ceil(i_idx);
		i_idx=i_idx+halfpool;
    }
    //cout << "therowindices: " << therowindices[0] << endl; 
	
    i_idx=0;
    for (int ct = 0; ct < width; ct++) {
		thecolindices[ct] =(int)ceil(i_idx);
		i_idx=i_idx+halfpool;
    }
	// cout << "thecolindices: " << thecolindices[1] << endl;
	
    cv::Mat max_val(length,width,CV_64FC1,0.0);
	
    int xcount = 0;
    int ycount = 0;
    int i,j;
    cv::Mat Awork_cv;
    for (int ct1 = 0; ct1 < length; ct1++) {
        for (int ct2= 0; ct2 < width; ct2++) {
			
            i=therowindices[ct1];
            j=thecolindices[ct2];
            cv::Range range1(i,min(i+fullpool,numrows)); // JSR corrected off-by-1: Aug 2014
			cv::Range range2(j,min(j+fullpool,numcols)); // JSR corrected off-by-1: Aug 2014
			//cout << i << "," << min(i+fullpool+1,numcols) << endl;          
			//cout << j << "," << min(j+fullpool+1,numcols) << endl;
			cv::Range ranges[] = {range1,range2};
            
            Im(ranges).copyTo(Awork_cv);
			
            //cout << "copied" << endl;
            //Awork_cv=Im(ranges);
            
			//cout << "rows Awork_cv: " << Awork_cv.rows << endl;
			//cout << "cols Awork_cv: " << Awork_cv.cols << endl; 
            
			// maxval(xcount,ycount) = max(max(image(i:min(i+poolrange,numrows),j:min(j+poolrange,numcols))));
			double minVal, maxVal;	
			cv::Point minLoc,maxLoc;
			cv::minMaxLoc(Awork_cv,&minVal,&maxVal,&minLoc,&maxLoc);
			//cout << "maxval: " << maxVal << endl;
			max_val.at<double>(xcount,ycount)=maxVal;
			
            ycount++;
        }
        xcount++;
        ycount=0;
    }
	
	dst=max_val;
}


void WindowedPatchDistance(cv::Mat& Im, cv::Mat& Patch, cv::Mat& dst, int doSkipSquareRoot) 
//dfd: version with doSkipSquareRoot flag.
{    
    //Estimate dimensions of Im and Patch:
    int numDim=Im.dims;
    int length = Im.size.p[0];
    int width = Im.size.p[1];
    int depth=Im.size.p[2];
    int matriz_sz = length * width;
    
    cout << "numDim: " << numDim << endl;
    cout << "depth: " << depth << endl;
    cout << "matriz_sz: " << matriz_sz << endl;
    
    //    int numDim_Patch=Patch.dims;
    int length_Patch = Patch.size.p[0];
    int width_Patch = Patch.size.p[1];
    int depth_Patch=Patch.size.p[2];
    //    int matriz_sz_Patch = length_Patch * width_Patch;
    
    if (depth != depth_Patch) {
        cout << "The third dimension of IM and Patch must be the same." << endl;//make this an error later
    }
    
    //Initialize variables:
    cv::Scalar Psqr;
    Psqr[0]=0;
    
    //Maybe check to make sure they are both of the same type.
    int type_data=Im.type();
    cout << "IM type:  " << type_data << endl;
    
    //Take 3D IM and Patch and save them into Awork and Bwork and the rest is the same.
    
    cv::Mat Awork_cv(length, width, type_data, 0);
    cv::Mat Awork_square(length, width, type_data, 0);
    
    cv::Range range1(0,length);
    cv::Range range2(0,width);
    
    cv::Mat Imsq(length,width, type_data, 0);
    cv::Mat Imsq_temp(length,width, type_data, 0);
    cv::Mat conv_temp(length,width, type_data, 0);
    
    cv::Mat PI(length,width, type_data, 0);
    cv::Mat PI_temp(length,width, type_data, 0);
    
    cv::Mat Bwork_cv(length_Patch, width_Patch,type_data, 0);
    cv::Mat Bwork_square(length_Patch, width_Patch,type_data, 0);
    
    cv::Range range1_Patch(0,length_Patch);
    cv::Range range2_Patch(0,width_Patch);//Reverse direction
    
    
    for (int i_chip = 0; i_chip < depth; i_chip++) {
        
        //Read slices of IM:
        cv::Range range3(i_chip,i_chip+1);
        cv::Range ranges[] = {range1,range2,range3};
        Im(ranges).copyTo(Awork_cv);//OR Aword_cv=Im(ranges)
        
        
        //Awork_cv.create(length,width,type_data);
        cv::Mat miniCubeA;
        miniCubeA.create(length,width,type_data);
        cv::resize(Awork_cv,miniCubeA,miniCubeA.size(),0,0);
        Awork_cv=miniCubeA;
        //        miniCubeA.release();
        
        //cout << "dimension Awork_cv:  " << Awork_cv.dims << endl;
        //cout << "element 2,3: " << Awork_cv.at<double>(2,3) << endl;
        
        //Read slices of Patch:
        cv::Range range3_Patch(i_chip,i_chip+1);
        cv::Range ranges_Patch[] = {range1_Patch,range2_Patch,range3_Patch};
        Patch(ranges_Patch).copyTo(Bwork_cv);//OR Aword_cv=Im(ranges)
        
        //Bwork_cv.create(length_Patch,width_Patch,type_data);
        cv::Mat miniCubeB;
        miniCubeB.create(length_Patch,width_Patch,type_data);
        cv::resize(Bwork_cv,miniCubeB,miniCubeB.size(),0,0);
        Bwork_cv=miniCubeB;
        //        miniCubeB.release();
        
        
        cv::flip(Bwork_cv, Bwork_cv, -1);  // flip prior to convolution
        
        //cout << "dimension Bwork_cv:  " << Bwork_cv.dims << endl;
        
        //Calculate PSQR
        cv::multiply(Bwork_cv, Bwork_cv, Bwork_square);
        Psqr+= cv::sum(Bwork_square);
        
        //cv::multiply(Awork_cv, Awork_cv, Awork_square);
        Awork_square=Awork_cv.mul(Awork_cv);
        if (i_chip == 0){
            //Imsq=Awork_square.clone();
            Awork_square.copyTo(Imsq);
        } else {
            Imsq_temp=Imsq;
            cv::add(Imsq_temp, Awork_square, Imsq);
        }
        
        
        //Create PI (conv2 stuff) 
        conv2(&Awork_cv, &Bwork_cv, &conv_temp);
        
        if (i_chip == 0){
            PI=conv_temp;
        } else {
            PI_temp=PI;
            cv::add(PI_temp, conv_temp, PI);
        }        
        
    } // end for loop in 3rd dimension
    
    cv::Scalar sum_support;
    sum_support[0]=(int)((double)length_Patch/2.0+0.5)-1;
    sum_support[1]=(int)((double)width_Patch/2.0+0.5)-1;
    sum_support[2]=(int)((double)length_Patch/2.0);
    sum_support[3]=(int)((double)width_Patch/2.0);
    
    //if prune==0
    cv::Mat Imsq_filter;
    sumfilter(&Imsq, sum_support, &Imsq_filter);
    Imsq=Imsq_filter;
    
    //else
    //    II = zeros(size(Im, 1), size(Im, 2));
    //PatchOrientations = +(Patch~=0);
    //for i = 1:dIm
    //        II = II + conv2(Imsq(:, :, i), PatchOrientations(:, :, i), 'same');
    //        end
    //                Imsq = II;
    //end
    
    cv::Mat D_temp(length, width, type_data, 0);
    cv::Mat D(length, width, type_data, 0);
    
    //Calculating:  Imsq - 2 * PI + Psqr
    D=Imsq-2*PI;
    D+=Psqr;
    
    cv::Mat tempmask;
    double thresh = 1e-14;
    cv::compare(D,thresh,tempmask,cv::CMP_GE);
    cv::divide(tempmask,255.0,tempmask);
    tempmask.convertTo(tempmask,TYPE);
    cv::multiply(D,tempmask,D);
    
    if (!doSkipSquareRoot) {
        sqrt(D, D_temp);
        // D(abs(D) < 10e-15) = 0;
        D=D_temp;
    }    
    //cout << "dimension D:  " << D.dims << endl;
    
    
    //     //double thresh = 10e-10;
    //     cv::Mat A = cv::abs(D) > thresh;    //(elements we want to keep are one, zeros otherwise)
    //     A.convertTo(A , CV_64FC1);          //(the > operator gives an 8-bit array, I believe
    //                                         //we have to convert it back to the data type of array D
    //                                         //for the next line to work)
    //     D = A.mul(D);
    
    dst=D; //OUTPUT
    
#if 0
    //Release mat data
    PI.release();
    Imsq.release();
    
    Awork_cv.release();
    Bwork_cv.release();
    
    Awork_square.release();
    Bwork_square.release();
    
    Imsq_temp.release();
    PI_temp.release();
    conv_temp.release();
    Imsq_filter.release();
#endif    
    //Should I release this?:
    // D_temp.release();
    // D.release();
    
}

//========================================================================================================

void WindowedPatchDistance(vector<cv::Mat>& Im, vector<cv::Mat>& Patch, cv::Mat& dst, int doSkipSquareRoot) 
//dfd: Vector version with doSkipSquareRoot flag.
{
    
    //Estimate dimensions of Im and Patch:
    int length = Im[0].rows; //Assuming all matrices have the same size
    int width = Im[0].cols;
    int depth=Im.size();
    
    //cout << "numDim: " << numDim << endl;
    //cout << "depth: " << depth << endl;
    //cout << "matriz_sz: " << matriz_sz << endl;
    
    int length_Patch = Patch[0].size.p[0];
    int width_Patch = Patch[0].size.p[1];
    int depth_Patch=Patch.size();
    
#if 1 // dfd: probably faster than if statement
    assert(depth == depth_Patch);
#else
    if (depth != depth_Patch) {
        cout << "The third dimension of IM and Patch must be the same." << endl;//make this an error later
        cv::waitKey(-1);
    }
#endif
    
    //Initialize variables:
    cv::Scalar Psqr;
    Psqr[0]=0;
    
    //Maybe check to make sure they are both of the same type.
    int type_data=Im[0].type();
    //cout << "IM type:  " << type_data << endl;
    //cout << "IM size: " << depth << endl;
    //cout << "Patch size: " << depth_Patch << endl;
    
    //Take 3D IM and Patch and save them into Awork and Bwork and the rest is the same.
    
    cv::Mat Awork_cv(length, width, type_data, 0);
    cv::Mat Awork_square(length, width, type_data, 0);
    
    
    cv::Mat Imsq(length,width, type_data, 0);
    cv::Mat Imsq_temp(length,width, type_data, 0);
    cv::Mat conv_temp(length,width, type_data, 0);
    
    cv::Mat PI(length,width, type_data, 0);
    cv::Mat PI_temp(length,width, type_data, 0);
    
    cv::Mat Bwork_cv(length_Patch, width_Patch,type_data, 0);
    cv::Mat Bwork_square(length_Patch, width_Patch,type_data, 0);
    
    for (int i_chip = 0; i_chip < depth; i_chip++) {
        
        //Read slices of IM:
        Im[i_chip].copyTo(Awork_cv);
        
        //cout << "dimension Awork_cv:  " << Awork_cv.dims << endl;
        //cout << "element 2,3: " << Awork_cv.at<double>(2,3) << endl;
        
        //Read slices of Patch:
        Patch[i_chip].copyTo(Bwork_cv);
        cv::flip(Bwork_cv, Bwork_cv, -1);  // flip prior to convolution
        
        //cout << "dimension Bwork_cv:  " << Bwork_cv.dims << endl;
        
        //Calculate PSQR
        cv::multiply(Bwork_cv, Bwork_cv, Bwork_square);
        Psqr+= cv::sum(Bwork_square);
        
        //cv::multiply(Awork_cv, Awork_cv, Awork_square);
        Awork_square=Awork_cv.mul(Awork_cv);
        if (i_chip == 0){
            //Imsq=Awork_square.clone();
            Awork_square.copyTo(Imsq);
        } else {
            Imsq_temp=Imsq;
            cv::add(Imsq_temp, Awork_square, Imsq);
        }
        
        
        //Create PI (conv2 stuff) 
        conv2(&Awork_cv, &Bwork_cv, &conv_temp);
        
        if (i_chip == 0){
            PI=conv_temp;
        } else {
            PI_temp=PI;
            cv::add(PI_temp, conv_temp, PI);
        }        
        
    } // end for loop in 3rd dimension
    
    cv::Scalar sum_support;
    sum_support[0]=(int)((double)length_Patch/2.0+0.5)-1;
    sum_support[1]=(int)((double)width_Patch/2.0+0.5)-1;
    sum_support[2]=(int)((double)length_Patch/2.0);
    sum_support[3]=(int)((double)width_Patch/2.0);
    
    //if prune==0
    cv::Mat Imsq_filter;
    sumfilter(&Imsq, sum_support, &Imsq_filter);
    Imsq=Imsq_filter;
    
    
    
    //else
    //    II = zeros(size(Im, 1), size(Im, 2));
    //PatchOrientations = +(Patch~=0);
    //for i = 1:dIm
    //        II = II + conv2(Imsq(:, :, i), PatchOrientations(:, :, i), 'same');
    //        end
    //                Imsq = II;
    //end
    
    cv::Mat D_temp(length, width, type_data, 0);
    cv::Mat D(length, width, type_data, 0);
    
    //Calculating:  Imsq - 2 * PI + Psqr
    D=Imsq-2*PI;
    D+=Psqr;
    
    cv::Mat tempmask;
    double thresh = 1e-14;
    cv::compare(D,thresh,tempmask,cv::CMP_GE);
    cv::divide(tempmask,255.0,tempmask);
    tempmask.convertTo(tempmask,TYPE);
    cv::multiply(D,tempmask,D);
    
	if (!doSkipSquareRoot) {
		sqrt(D, D_temp);
		// D(abs(D) < 10e-15) = 0;
		D=D_temp;
	}    
    //cout << "dimension D:  " << D.dims << endl;
    
    
    //     //double thresh = 10e-10;
    //     cv::Mat A = cv::abs(D) > thresh;    //(elements we want to keep are one, zeros otherwise)
    //     A.convertTo(A , CV_64FC1);          //(the > operator gives an 8-bit array, I believe
    //                                         //we have to convert it back to the data type of array D
    //                                         //for the next line to work)
    //     D = A.mul(D);
    
    dst=D; //OUTPUT
    
#if 0
    //Release mat data
    PI.release();
    Imsq.release();
    
    Awork_cv.release();
    Bwork_cv.release();
    
    Awork_square.release();
    Bwork_square.release();
    
    Imsq_temp.release();
    PI_temp.release();
    conv_temp.release();
    Imsq_filter.release();
    
    //Should I release this?:
    // D_temp.release();
    // D.release();
    
#endif    
}

//========================================================================================================
