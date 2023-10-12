#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>
#include <vector>
#include <math.h>
#include <complex>
#include <pcl/common/common.h>
#include <pcl/io/pcd_io.h>
#include <pcl/conversions.h>
#include <pcl/point_cloud.h>
#include <pcl/point_types.h>
#include <pcl/pcl_macros.h>
#include <pcl/kdtree/kdtree_flann.h>
#include <pcl/common/transforms.h>


using namespace std;

struct PointXYZIO
{
    PCL_ADD_POINT4D
    PCL_ADD_INTENSITY;
    int order;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
} EIGEN_ALIGN16;
POINT_CLOUD_REGISTER_POINT_STRUCT (PointXYZIO,
    (float, x, x) (float, y, y) (float, z, z) (float, intensity, intensity)
    (int, order, order)
)

typedef PointXYZIO PointOrder;




void TrajectoryFilterAndZoom(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut1, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut2, float StepLenght, float FactorXY, float FactoryZ) {
	int CloudSize = CloudIn->size();
	int IndexCloudOut1 = 0;
    	float PointDistance = 0;
    	float PointDistanceSquare = 0;
    	CloudOut1->resize(CloudSize);
    	CloudOut2->resize(CloudSize);
      
    	CloudOut1->points[IndexCloudOut1].x = CloudIn->points[0].x;
    	CloudOut1->points[IndexCloudOut1].y = CloudIn->points[0].y;
    	CloudOut1->points[IndexCloudOut1].z = CloudIn->points[0].z;
    	CloudOut1->points[IndexCloudOut1].intensity = CloudIn->points[0].intensity;
    	IndexCloudOut1 = IndexCloudOut1 + 1;
    
    	CloudOut2->points[0].x = CloudIn->points[0].x * FactorXY;
    	CloudOut2->points[0].y = CloudIn->points[0].y * FactorXY;
    	CloudOut2->points[0].z = CloudIn->points[0].z * FactoryZ;
    	CloudOut2->points[0].intensity = CloudIn->points[0].intensity;
    
    
    	for (int i = 1; i < CloudSize; ++i) {
    		CloudOut2->points[i].x = CloudIn->points[i].x * FactorXY;
    		CloudOut2->points[i].y = CloudIn->points[i].y * FactorXY;
    		CloudOut2->points[i].z = CloudIn->points[i].z * FactoryZ;
    		CloudOut2->points[i].intensity = CloudIn->points[i].intensity;
    	
    		PointDistanceSquare = ((CloudIn->points[i].x - CloudOut1->points[IndexCloudOut1 - 1].x) * (CloudIn->points[i].x - CloudOut1->points[IndexCloudOut1 - 1].x)) + ((CloudIn->points[i].y - CloudOut1->points[IndexCloudOut1 - 1].y) * (CloudIn->points[i].y - CloudOut1->points[IndexCloudOut1 - 1].y)) + ((CloudIn->points[i].z - CloudOut1->points[IndexCloudOut1 - 1].z) * (CloudIn->points[i].z - CloudOut1->points[IndexCloudOut1 - 1].z));
    		PointDistance = sqrt(PointDistanceSquare);
    		if (PointDistance >= StepLenght){
    			CloudOut1->points[IndexCloudOut1].x = CloudIn->points[i].x;
    			CloudOut1->points[IndexCloudOut1].y = CloudIn->points[i].y;
    			CloudOut1->points[IndexCloudOut1].z = CloudIn->points[i].z;
    			CloudOut1->points[IndexCloudOut1].intensity = CloudIn->points[i].intensity;
    			IndexCloudOut1 = IndexCloudOut1 + 1;
		}		
    	}
    CloudOut1->resize(IndexCloudOut1);    
}

void TrajectorySmooth(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut, int MaxSmoothTimes) {
	pcl::PointCloud<pcl::PointXYZI>::Ptr Cloud(new pcl::PointCloud<pcl::PointXYZI>);
	int CloudSize = CloudIn->size();
    	CloudOut->resize(CloudSize);
    	*Cloud = *CloudIn;
    	
    	pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdT;
    	std::vector<int> index;
    	std::vector<float> distance;
    	pcl::PointXYZI SmoothPoint;
    	float msx, msy, msz;
    
    	for (int SmoothTimes = 0; SmoothTimes < MaxSmoothTimes; ++SmoothTimes){
    		kdT.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());  
    		kdT->setInputCloud(Cloud);
    		for (int i = 0; i < CloudSize; ++i) {
    			SmoothPoint.x = Cloud->points[i].x;
    			SmoothPoint.y = Cloud->points[i].y;
    			SmoothPoint.z = Cloud->points[i].z;
    	
			kdT->nearestKSearch(SmoothPoint, 5, index, distance);
			msx = 0;
			msy = 0;
			msz = 0;
			for (int j = 0; j < 5; ++j){
    				msx = msx + Cloud->points[index[j]].x;
    				msy = msy + Cloud->points[index[j]].y;
    				msz = msz + Cloud->points[index[j]].z;
    			}
    			CloudOut->points[i].x = msx / 5.00;
    			CloudOut->points[i].y = msy / 5.00;
    			CloudOut->points[i].z = msz / 5.00;
    			CloudOut->points[i].intensity = Cloud->points[i].intensity;	
    		}
    
    		for (int i = 0; i < CloudSize; ++i) {
    			CloudOut->points[i].x = CloudOut->points[i].x - (CloudOut->points[0].x - Cloud->points[0].x);
    			CloudOut->points[i].y = CloudOut->points[i].y - (CloudOut->points[0].y - Cloud->points[0].y);
    			CloudOut->points[i].z = CloudOut->points[i].z - (CloudOut->points[0].z - Cloud->points[0].z);
    			CloudOut->points[i].intensity = CloudOut->points[i].intensity;		
    		}
    
    		//reset trajectoryfiltered, smooth again
    		for (int i = 0; i < CloudSize; ++i) {
    			Cloud->points[i].x = CloudOut->points[i].x;
    			Cloud->points[i].y = CloudOut->points[i].y;
    			Cloud->points[i].z = CloudOut->points[i].z;
    			Cloud->points[i].intensity = CloudOut->points[i].intensity;		
    		}
    	}
	    
}

void DensifyAndProjectTrajectory(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut1, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut2, float MinDistance) {
	
	int CloudInSize = CloudIn->size();
	int CloudOutSize = CloudInSize * (50 / MinDistance);
    	CloudOut1->resize(CloudOutSize);
    	CloudOut2->resize(CloudOutSize);
    
    	int IndexCloudOut = 0;
    	int InterNum = 0;
    	float s1 = 0;
    	float s1Square = 0;
    	float sx = 0;
    	float sy = 0;
    	float sz = 0;
    		
    	for (int i = 0; i < CloudInSize - 1; ++i) {
		s1Square = (CloudIn->points[i + 1].x - CloudIn->points[i].x) * (CloudIn->points[i + 1].x - CloudIn->points[i].x) + (CloudIn->points[i + 1].y - CloudIn->points[i].y) * (CloudIn->points[i + 1].y - CloudIn->points[i].y) + (CloudIn->points[i + 1].z - CloudIn->points[i].z) * (CloudIn->points[i + 1].z - CloudIn->points[i].z);
    		s1 = sqrt(s1Square);
    		InterNum = s1 / MinDistance;
		if (InterNum == 0) {
    			CloudOut1->points[IndexCloudOut].x = CloudIn->points[i].x;
    			CloudOut1->points[IndexCloudOut].y = CloudIn->points[i].y;
    			CloudOut1->points[IndexCloudOut].z = CloudIn->points[i].z;
    			CloudOut1->points[IndexCloudOut].intensity = CloudIn->points[i].intensity;
    			
    			CloudOut2->points[IndexCloudOut].x = CloudIn->points[i].x;
    			CloudOut2->points[IndexCloudOut].y = CloudIn->points[i].y;
    			CloudOut2->points[IndexCloudOut].z = 0;
    			CloudOut2->points[IndexCloudOut].intensity = IndexCloudOut;
    			IndexCloudOut = IndexCloudOut + 1;
    		}
    		if (InterNum > 0) {
    	   		sx = (CloudIn->points[i + 1].x - CloudIn->points[i].x) / InterNum;
    	   		sy = (CloudIn->points[i + 1].y - CloudIn->points[i].y) / InterNum;
    	   		sz = (CloudIn->points[i + 1].z - CloudIn->points[i].z) / InterNum;
    	   		for (int k = 0; k < InterNum; ++k) {
    	   			CloudOut1->points[IndexCloudOut].x = CloudIn->points[i].x + (k * sx);
    				CloudOut1->points[IndexCloudOut].y = CloudIn->points[i].y + (k * sy);
    				CloudOut1->points[IndexCloudOut].z = CloudIn->points[i].z + (k * sz);
    				CloudOut1->points[IndexCloudOut].intensity = CloudIn->points[i].intensity;
    				
    				CloudOut2->points[IndexCloudOut].x = CloudIn->points[i].x + (k * sx);
    				CloudOut2->points[IndexCloudOut].y = CloudIn->points[i].y + (k * sy);
    				CloudOut2->points[IndexCloudOut].z = 0;
    				CloudOut2->points[IndexCloudOut].intensity = IndexCloudOut;
    				IndexCloudOut = IndexCloudOut + 1;
    	   		}   	
    		} 			
    		
    	}
    	CloudOut1->resize(IndexCloudOut);
    	CloudOut2->resize(IndexCloudOut);	    
}

void SortMap(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn1, const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn2, const pcl::PointCloud<PointOrder>::Ptr &CloudOut) {
	
	int CloudInSize = CloudIn1->size();
	CloudOut->resize(CloudInSize);
     	pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdT;
    	kdT.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());  
    	kdT->setInputCloud(CloudIn2);
    	pcl::PointXYZI ReferPoint;
    	std::vector<int> index;
    	std::vector<float> distance;
    
   	for (int i = 0; i < CloudInSize; ++i) {
    		ReferPoint.x = CloudIn1->points[i].x;
    		ReferPoint.y = CloudIn1->points[i].y;
    		ReferPoint.z = CloudIn1->points[i].z;
    	
		kdT->nearestKSearch(ReferPoint, 1, index, distance);
	
		CloudOut->points[i].x = CloudIn1->points[i].x;
    		CloudOut->points[i].y = CloudIn1->points[i].y;
    		CloudOut->points[i].z = CloudIn1->points[i].z;
    		CloudOut->points[i].intensity = CloudIn1->points[i].intensity;
    		CloudOut->points[i].order = CloudIn2->points[index[0]].intensity;
	
    }
    std::sort(CloudOut->points.begin(), CloudOut->points.end(), [](const PointOrder& a, const PointOrder& b) { return a.order < b.order; });
}

void CalculateDensity(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut, int CalculatePointNumber) {
	int CloudSize = CloudIn->size();
    	pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdT;
    	kdT.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());  
    	kdT->setInputCloud(CloudIn);
    	pcl::PointXYZI ReferPoint;
    	std::vector<int> index;
    	std::vector<float> distance;
    	float SumDensity = 0;
    	float AverageDensity = 0;
    	CloudOut->resize(CloudSize);
    	
    	for (int i = 0; i < CloudSize; ++i) {
    		ReferPoint.x = CloudIn->points[i].x;
    		ReferPoint.y = CloudIn->points[i].y;
    		ReferPoint.z = CloudIn->points[i].z;
    		kdT->nearestKSearch(ReferPoint, CalculatePointNumber, index, distance);
    		for (int j = 0; j < CalculatePointNumber; ++j){
    			SumDensity = SumDensity + abs((CloudIn->points[index[j]].x - CloudIn->points[i].x) * (CloudIn->points[index[j]].y - CloudIn->points[i].y) * (CloudIn->points[index[j]].z - CloudIn->points[i].z));	
    		}
    		AverageDensity = SumDensity / (CalculatePointNumber - 1);
    		CloudOut->points[i].x = CloudIn->points[i].x;
		CloudOut->points[i].y = CloudIn->points[i].y;
    		CloudOut->points[i].z = CloudIn->points[i].z;
    		CloudOut->points[i].intensity = AverageDensity;
    		
    		SumDensity = 0;
    		AverageDensity = 0;		
    	}

}

void SmoothCloud(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut, int CalculatePointNumber) {
	int CloudSize = CloudIn->size();
    	pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdT;
    	kdT.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());  
    	kdT->setInputCloud(CloudIn);
    	pcl::PointXYZI ReferPoint;
    	std::vector<int> index;
    	std::vector<float> distance;
    	CloudOut->resize(CloudSize);
    	float SumDensity = 0;
    	float dx = 0;
    	float dy = 0;
    	float dz = 0;
    	for (int i = 0; i < CloudSize; ++i) {
    		ReferPoint.x = CloudIn->points[i].x;
    		ReferPoint.y = CloudIn->points[i].y;
    		ReferPoint.z = CloudIn->points[i].z;
    		kdT->nearestKSearch(ReferPoint, CalculatePointNumber, index, distance);
    		
		for (int j = 0; j < CalculatePointNumber; ++j){
    			dx = dx + (CloudIn->points[index[j]].x / CloudIn->points[index[j]].intensity);
    			dy = dy + (CloudIn->points[index[j]].y / CloudIn->points[index[j]].intensity);
    			dz = dz + (CloudIn->points[index[j]].z / CloudIn->points[index[j]].intensity);
    			SumDensity = SumDensity + (1.0 / CloudIn->points[index[j]].intensity);
    		}
    		CloudOut->points[i].x = dx / SumDensity;
    		CloudOut->points[i].y = dy / SumDensity;
    		CloudOut->points[i].z = dz / SumDensity;
    		CloudOut->points[i].intensity = CloudIn->points[i].intensity;
    		dx = 0;
    		dy = 0;
    		dz = 0;
    		SumDensity = 0;
    	}
}

void DeleteRepeatedPoint(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut, int CalculatePointNumber, float MinDistance) {
	int CloudSize = CloudIn->size();
    	pcl::KdTreeFLANN<pcl::PointXYZI>::Ptr kdT;
    	kdT.reset(new pcl::KdTreeFLANN<pcl::PointXYZI>());  
    	kdT->setInputCloud(CloudIn);
    	pcl::PointXYZI ReferPoint;
    	std::vector<int> index;
    	std::vector<float> distance;
    	CloudOut->resize(CloudSize);
    	float MaxDensity = 0;
    	int IndexCloudOut = 0;

    	for (int i = 0; i < CloudSize; ++i) {
    		ReferPoint.x = CloudIn->points[i].x;
    		ReferPoint.y = CloudIn->points[i].y;
    		ReferPoint.z = CloudIn->points[i].z;
    		kdT->nearestKSearch(ReferPoint, CalculatePointNumber, index, distance);
    		MaxDensity = CloudIn->points[i].intensity;
    		for( int j = 0; j < CalculatePointNumber; ++j){
    			if (MaxDensity < CloudIn->points[index[j]].intensity){
    				MaxDensity = CloudIn->points[index[j]].intensity;
    			}
    		}
    		if (distance[1] >= MinDistance){
    			CloudOut->points[IndexCloudOut].x = CloudIn->points[i].x;
    			CloudOut->points[IndexCloudOut].y = CloudIn->points[i].y;
    			CloudOut->points[IndexCloudOut].z = CloudIn->points[i].z;
    			CloudOut->points[IndexCloudOut].intensity = CloudIn->points[i].intensity;
    			IndexCloudOut = IndexCloudOut + 1;
    			}
    		if (distance[1] < MinDistance && MaxDensity == CloudIn->points[i].intensity){
    			CloudOut->points[IndexCloudOut].x = CloudIn->points[i].x;
    			CloudOut->points[IndexCloudOut].y = CloudIn->points[i].y;
    			CloudOut->points[IndexCloudOut].z = CloudIn->points[i].z;
    			CloudOut->points[IndexCloudOut].intensity = CloudIn->points[i].intensity;
    			IndexCloudOut = IndexCloudOut + 1;
    		}
    	}
    	CloudOut->resize(IndexCloudOut);
}

void SortCrossSection(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut) {
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud1(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud2(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointCloud<pcl::PointXYZI>::Ptr cloud3(new pcl::PointCloud<pcl::PointXYZI>);
	int CloudSize = CloudIn->size();
    	cloud1->resize(CloudSize);
     	cloud2->resize(CloudSize);
     	cloud3->resize(CloudSize);
    	CloudOut->resize(CloudSize);
    	pcl::PointXYZI MinXYZ;
    	pcl::PointXYZI MaxXYZ;
    	pcl::getMinMax3D(*CloudIn, MinXYZ, MaxXYZ);
    	float LeftPart = MinXYZ.y + (MaxXYZ.y -MinXYZ.y) / 4;
     	float RightPart = MaxXYZ.y - (MaxXYZ.y -MinXYZ.y) / 4;
     	int IndexLeft = 0;
     	int IndexMid = 0;
     	int IndexRight = 0;
     	for (int i = 0; i < CloudSize; ++i) {
		if (CloudIn->points[i].y <= LeftPart) {
			cloud1->points[IndexLeft].x = 0;
			cloud1->points[IndexLeft].y = CloudIn->points[i].y;
			cloud1->points[IndexLeft].z = CloudIn->points[i].z;
			cloud1->points[IndexLeft].intensity = CloudIn->points[i].intensity;
			IndexLeft += 1;
		}
		if (CloudIn->points[i].y > LeftPart && CloudIn->points[i].y <= RightPart) {
			cloud2->points[IndexMid].x = 0;
			cloud2->points[IndexMid].y = CloudIn->points[i].y;
			cloud2->points[IndexMid].z = CloudIn->points[i].z;
			cloud2->points[IndexMid].intensity = CloudIn->points[i].intensity;
			IndexMid += 1;
		}
		if (CloudIn->points[i].y > RightPart) {
			cloud3->points[IndexRight].x = 0;
			cloud3->points[IndexRight].y = CloudIn->points[i].y;
			cloud3->points[IndexRight].z = CloudIn->points[i].z;
			cloud3->points[IndexRight].intensity = CloudIn->points[i].intensity;
			IndexRight += 1;
		}	
     	}
     	cloud1->resize(IndexLeft);
     	cloud2->resize(IndexMid);
     	cloud3->resize(IndexRight);
     	std::sort(cloud1->points.begin(), cloud1->points.end(), [](const pcl::PointXYZI& a, const pcl::PointXYZI& b) { return a.z > b.z; });
     	std::sort(cloud2->points.begin(), cloud2->points.end(), [](const pcl::PointXYZI& a, const pcl::PointXYZI& b) { return a.y < b.y; });
     	std::sort(cloud3->points.begin(), cloud3->points.end(), [](const pcl::PointXYZI& a, const pcl::PointXYZI& b) { return a.z < b.z; });
     	*CloudOut = *cloud1 + *cloud2 + *cloud3;
}

void DensifyCrossSection(const pcl::PointCloud<pcl::PointXYZI>::Ptr& CloudIn, const pcl::PointCloud<pcl::PointXYZI>::Ptr &CloudOut, float MinDistance, float Hight) {
	pcl::PointCloud<pcl::PointXYZI>::Ptr CloudFiltered(new pcl::PointCloud<pcl::PointXYZI>);
	pcl::PointXYZI MinXYZ;
    	pcl::PointXYZI MaxXYZ;
    	pcl::getMinMax3D(*CloudIn, MinXYZ, MaxXYZ);
    	float MaxHight = Hight + MinXYZ.z;
	int CloudSize = CloudIn->size();
	CloudFiltered->resize(CloudSize);
	int IndexCloudFiltered = 0;
	for (int i = 0; i < CloudSize; ++i) {
		if (CloudIn->points[i].z <= MaxHight){
    			CloudFiltered->points[IndexCloudFiltered].x = CloudIn->points[i].x;
    			CloudFiltered->points[IndexCloudFiltered].y = CloudIn->points[i].y;
    			CloudFiltered->points[IndexCloudFiltered].z = CloudIn->points[i].z;
    			CloudFiltered->points[IndexCloudFiltered].intensity = CloudIn->points[i].intensity;
    			IndexCloudFiltered = IndexCloudFiltered + 1;
    		}
	}
	CloudFiltered->resize(IndexCloudFiltered);
	
	
    	int PreSize = (300 / MinDistance) + IndexCloudFiltered;
    	CloudOut->resize(PreSize);
    	int IndexCloudOut = 0;
    	int InterNum = 0;
    	float s1 = 0;
    	float s1Square = 0;
    	float sx = 0;
    	float sy = 0;
    	float sz = 0;
    		
    	for (int i = 0; i < IndexCloudFiltered - 1; ++i) {
		s1Square = (CloudFiltered->points[i + 1].y - CloudFiltered->points[i].y) * (CloudFiltered->points[i + 1].y - CloudFiltered->points[i].y) + (CloudFiltered->points[i + 1].z - CloudFiltered->points[i].z) * (CloudFiltered->points[i + 1].z - CloudFiltered->points[i].z);
    		s1 = sqrt(s1Square);
    		InterNum = s1 / MinDistance;
		if (InterNum == 0) {
    			CloudOut->points[IndexCloudOut].x = CloudFiltered->points[i].x;
    			CloudOut->points[IndexCloudOut].y = CloudFiltered->points[i].y;
    			CloudOut->points[IndexCloudOut].z = CloudFiltered->points[i].z;
    			CloudOut->points[IndexCloudOut].intensity = CloudFiltered->points[i].intensity;
    			IndexCloudOut = IndexCloudOut + 1;
    		}
    		if (InterNum > 0) {
    	   		sx = (CloudFiltered->points[i + 1].x - CloudFiltered->points[i].x) / InterNum;
    	   		sy = (CloudFiltered->points[i + 1].y - CloudFiltered->points[i].y) / InterNum;
    	   		sz = (CloudFiltered->points[i + 1].z - CloudFiltered->points[i].z) / InterNum;
    	   		for (int k = 0; k < InterNum; ++k) {
    	   			CloudOut->points[IndexCloudOut].x = CloudFiltered->points[i].x + (k * sx);
    				CloudOut->points[IndexCloudOut].y = CloudFiltered->points[i].y + (k * sy);
    				CloudOut->points[IndexCloudOut].z = CloudFiltered->points[i].z + (k * sz);
    				CloudOut->points[IndexCloudOut].intensity = CloudFiltered->points[i].intensity;
    				IndexCloudOut = IndexCloudOut + 1;
    	   		}   	
    		} 			
    		
    	}
    	CloudOut->resize(IndexCloudOut);
}


int main(int argc,char *argv[])
{
    std::cout << "PCD ZOOM ELEVATION START" << std::endl;
    
    std::cout << "0 CONFIGURE PARAMETERS" << std::endl;
    float referXY = 1000.00;
    float referZ = 100.00;
    std::cout << "INPUT REFERED XY (m)" << std::endl;
    std::cin >> referXY;
    std::cout << "INPUT REFERED Z (m)" << std::endl;
    std::cin >> referZ;
    float dtra0 = 10;
    std::cout << "INPUT MINI DISRANCE OF TWO INITIAL TRAJECTORY POINTS (m)" << std::endl;
    std::cin >> dtra0;
    float dtra1 = 1;
    std::cout << "INPUT MINI DISRANCE OF TWO DENSIFIED TRAJECTORY POINTS (cm)" << std::endl;
    std::cin >> dtra1;
    dtra1 = dtra1 * 0.01;
    int smoothnumber = 1;
    std::cout << "INPUT SMOOTH TRAJECTORY TIME" << std::endl;
    std::cin >> smoothnumber;
    int intensitynum = 100;
    std::cout << "INPUT DENSITY CALCULATE POINT NUMBER" << std::endl;
    std::cin >> intensitynum;
    int smNUM = 10;
    std::cout << "INPUT NEAREST POINT NUMBER" << std::endl;
    std::cin >> smNUM;
    int dividelenght0 = 20;
    std::cout << "INPUT DIVISION LENGHT (cm) AT LEAST 10cm" << std::endl;
    std::cin >> dividelenght0;
    float zmax = 30;
    std::cout << "INPUT CROSS-SECTION HIGHT (m)" << std::endl;
    std::cin >> zmax;
    float denrate = 10;
    std::cout << "INPUT CROSS-SECTION DENSIFY RATE (cm)" << std::endl;
    std::cin >> denrate;
    denrate = denrate * 0.01;
    
    
    std::string savepath = "/home/thierry/Downloads/collected-data/reconstruction/workflow3/data/dividedgroup/";
    
    
    std::cout << "1 FILTER AND ZOOM TRAJECTORY" << std::endl;
    pcl::PointCloud<pcl::PointXYZI>::Ptr trajectory(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr trajectoryzoom(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr trajectoryfiltered(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::io::loadPCDFile(savepath + "trajectory.pcd", *trajectory);  //read trajectory pcd
    
    int trajectorysize = trajectory->size();
    
    float zoomfactorXY = 1.00;
    float zoomfactor = 1.00;
    
    zoomfactorXY = referXY / (sqrt((((trajectory->points[trajectorysize - 1].x + trajectory->points[trajectorysize - 2].x + trajectory->points[trajectorysize - 3].x) / 3) * ((trajectory->points[trajectorysize - 1].x + trajectory->points[trajectorysize - 2].x + trajectory->points[trajectorysize - 3].x) / 3)) + (((trajectory->points[trajectorysize - 1].y + trajectory->points[trajectorysize - 2].y + trajectory->points[trajectorysize - 3].y) / 3) * ((trajectory->points[trajectorysize - 1].y + trajectory->points[trajectorysize - 2].y + trajectory->points[trajectorysize - 3].y) / 3))));
    zoomfactor = referZ / ((trajectory->points[trajectorysize - 1].z + trajectory->points[trajectorysize - 2].z + trajectory->points[trajectorysize - 3].z) / 3);
    
    TrajectoryFilterAndZoom(trajectory, trajectoryfiltered, trajectoryzoom, dtra0, zoomfactorXY, zoomfactor);
    pcl::io::savePCDFileBinary(savepath + "trajectoryfiltered.pcd", *trajectoryfiltered); //save trajectoryfiltered pcd
    pcl::io::savePCDFileBinary(savepath + "trajectoryzoom.pcd", *trajectoryzoom); //save trajectoryzoom pcd
    std::cout << "FILTER AND ZOOM TRAJECTORY FINISHED" << std::endl;
    
    
    std::cout << "2 SMOOTH FILTERED TRAJECTORY START" << std::endl;
    pcl::PointCloud<pcl::PointXYZI>::Ptr smoothtrajectoryfiltered(new pcl::PointCloud<pcl::PointXYZI>);
    
    TrajectorySmooth(trajectoryfiltered, smoothtrajectoryfiltered, smoothnumber);
    
    pcl::io::savePCDFileBinary(savepath + "smoothtrajectoryfiltered.pcd", *smoothtrajectoryfiltered); //save smoothtrajectoryfiltered pcd
    std::cout << "SMOOTH FILTERED TRAJECTORY FINISHED" << std::endl;
    

    std::cout << "3 DENSIFY SMOOTH FILTERED TRAJECTORY START" << std::endl;
    pcl::PointCloud<pcl::PointXYZI>::Ptr densifysmoothtrajectoryfiltered(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<pcl::PointXYZI>::Ptr densifysmoothtrajectoryfilteredprojection(new pcl::PointCloud<pcl::PointXYZI>);
    
    DensifyAndProjectTrajectory(smoothtrajectoryfiltered, densifysmoothtrajectoryfiltered, densifysmoothtrajectoryfilteredprojection, dtra1);
    
    pcl::io::savePCDFileBinary(savepath + "densifysmoothtrajectoryfiltered.pcd", *densifysmoothtrajectoryfiltered); //save densifysmoothtrajectoryfiltered pcd
    pcl::io::savePCDFileBinary(savepath + "densifysmoothtrajectoryfilteredprojection.pcd", *densifysmoothtrajectoryfilteredprojection); //save densifysmoothtrajectoryfilteredprojection pcd
    std::cout << "DENSIFY SMOOTH FILTERED TRAJECTORY FINISHED" << std::endl;
    
    
    
    std::cout << "4 GlobalMap SORT START" << std::endl;
    pcl::PointCloud<pcl::PointXYZI>::Ptr GlobalMap(new pcl::PointCloud<pcl::PointXYZI>);
    pcl::PointCloud<PointOrder>::Ptr GlobalMaporder(new pcl::PointCloud<PointOrder>);
    pcl::io::loadPCDFile(savepath + "GlobalMap.pcd", *GlobalMap);  //read GlobalMap pcd
    
    SortMap(GlobalMap, densifysmoothtrajectoryfilteredprojection, GlobalMaporder);
    
    pcl::io::savePCDFileBinary(savepath + "GlobalMaporder.pcd", *GlobalMaporder); //save GlobalMaporder pcd
    std::cout << "GlobalMap SORT FINISHED" << std::endl;
    
    
    std::cout << "5 GlobalMap ZOOM START" << std::endl;
    pcl::PointCloud<pcl::PointXYZI>::Ptr GlobalMapzoom(new pcl::PointCloud<pcl::PointXYZI>);
    int GlobalMapzoomsize = GlobalMaporder->size();
    GlobalMapzoom->resize(GlobalMapzoomsize);
    
    int zoomindex = 0;
    for (int i = 0; i < GlobalMapzoomsize; ++i) {
    	    zoomindex = GlobalMaporder->points[i].order;
    	    
    	    GlobalMapzoom->points[i].x = GlobalMaporder->points[i].x - (densifysmoothtrajectoryfiltered->points[zoomindex].x * (1.00 - zoomfactorXY));
    	    GlobalMapzoom->points[i].y = GlobalMaporder->points[i].y - (densifysmoothtrajectoryfiltered->points[zoomindex].y * (1.00 - zoomfactorXY));
    	    GlobalMapzoom->points[i].z = GlobalMaporder->points[i].z - (densifysmoothtrajectoryfiltered->points[zoomindex].z * (1.00 - zoomfactor));
    	    
    	    GlobalMapzoom->points[i].intensity = GlobalMaporder->points[i].intensity;
    	        	    
    }
    pcl::io::savePCDFileBinary(savepath + "GlobalMapzoom.pcd", *GlobalMapzoom); //save GlobalMapzoom pcd
    std::cout << "GlobalMap ZOOM FINISHED" << std::endl;
    
    
    
    std::cout << "6 CALCULATE POINT DENSITY START" << std::endl;  
    pcl::PointCloud<pcl::PointXYZI>::Ptr intensitycloudOut(new pcl::PointCloud<pcl::PointXYZI>);
    
    CalculateDensity(GlobalMapzoom, intensitycloudOut, intensitynum);
	
    pcl::io::savePCDFileBinary(savepath + "intensitycloudOut.pcd", *intensitycloudOut); //save intensitycloudOut pcd
    std::cout << "CALCULATE POINT DENSITY FINISHED" << std::endl;
    
    std::cout << "7 SMOOTH PCD START" << std::endl;
    pcl::PointCloud<pcl::PointXYZI>::Ptr GlobalMapzoomsmooth(new pcl::PointCloud<pcl::PointXYZI>);
    
    SmoothCloud(intensitycloudOut, GlobalMapzoomsmooth, smNUM);
    
    pcl::io::savePCDFileBinary(savepath + "GlobalMapzoomsmooth.pcd", *GlobalMapzoomsmooth); //save GlobalMapzoomsmooth pcd
    std::cout << "PCD SMOOTH FINISHED" << std::endl;
    
    
    std::cout << "8 GlobalMapzoomsmooth DIVISION START" << std::endl;
    
    pcl::PointCloud<pcl::PointXYZI>::Ptr dividedcloud(new pcl::PointCloud<pcl::PointXYZI>);
    int dividedcloudsize = GlobalMapzoomsize;
    dividedcloud->resize(dividedcloudsize);
    
    
    int dividelenght = 0;
    
    
    
    int indextrajectoryi;
    std::string indextrajectory;
    std::string indextrajectoryx;
    std::string indextrajectoryy;
    std::string indextrajectoryz;
    std::string indexsection;
    int i0 = 0;
    int CrosssectionNumber = 1;
    float paraa = 0;
    float parab = 0;
    float parac = 0;
    float parat = 0;
    float paraab = 0;
    float pointix = 0;
    float pointiy = 0;
    
    for (int i = 0; i < GlobalMapzoomsize - 1; ++i) {
    	    dividedcloud->points[i0].x = GlobalMapzoomsmooth->points[i].x;
    	    dividedcloud->points[i0].y = GlobalMapzoomsmooth->points[i].y;
    	    dividedcloud->points[i0].z = GlobalMapzoomsmooth->points[i].z;
    	    dividedcloud->points[i0].intensity = GlobalMapzoomsmooth->points[i].intensity;
    	    i0 = i0 + 1;
    	    if (GlobalMaporder->points[i].order != GlobalMaporder->points[i + 1].order || i + 2 == GlobalMapzoomsize) {
    	       dividelenght = dividelenght + 1;
    	       if (dividelenght == dividelenght0){
    	       dividedcloudsize = i0;
    	       dividedcloud->resize(dividedcloudsize);
    	       indextrajectoryi = GlobalMaporder->points[i].order - (dividelenght0 / 2);
    	       indextrajectory = std::to_string(indextrajectoryi);
    	       indextrajectoryx = std::to_string(densifysmoothtrajectoryfiltered->points[indextrajectoryi].x * zoomfactorXY);
    	       indextrajectoryy = std::to_string(densifysmoothtrajectoryfiltered->points[indextrajectoryi].y * zoomfactorXY);
    	       indextrajectoryz = std::to_string(densifysmoothtrajectoryfiltered->points[indextrajectoryi].z * zoomfactor);
    	       pcl::io::savePCDFileBinary(savepath + "dividedcloud/Order" + indextrajectory + "X" + indextrajectoryx + "Y" + indextrajectoryy + "Z" + indextrajectoryz +".pcd", *dividedcloud);
    	       
    	       pcl::PointCloud<pcl::PointXYZI>::Ptr projectiondividedcloud(new pcl::PointCloud<pcl::PointXYZI>);
    	       pcl::PointCloud<pcl::PointXYZI>::Ptr fprojectiondividedcloud(new pcl::PointCloud<pcl::PointXYZI>);
    	       projectiondividedcloud->resize(i0);
    	       fprojectiondividedcloud->resize(i0);
    	       for (int p0 = 0; p0 < i0; ++p0){
    	       	paraa = densifysmoothtrajectoryfiltered->points[indextrajectoryi + (dividelenght0 / 2)].x - densifysmoothtrajectoryfiltered->points[indextrajectoryi - (dividelenght0 / 2)].x;
    	       	parab = densifysmoothtrajectoryfiltered->points[indextrajectoryi + (dividelenght0 / 2)].y - densifysmoothtrajectoryfiltered->points[indextrajectoryi - (dividelenght0 / 2)].y;
    	       	// parac = densifysmoothtrajectoryfiltered->points[indextrajectoryi + (dividelenght0 / 2)].z - densifysmoothtrajectoryfiltered->points[indextrajectoryi - (dividelenght0 / 2)].z;
    	       	parac = 0;
    	       	parat = (paraa * (densifysmoothtrajectoryfiltered->points[indextrajectoryi].x - dividedcloud->points[p0].x) +parab * (densifysmoothtrajectoryfiltered->points[indextrajectoryi].y - dividedcloud->points[p0].y) + parac * (densifysmoothtrajectoryfiltered->points[indextrajectoryi].z - dividedcloud->points[p0].z)) / ((paraa *paraa) + (parab *parab) +(parac *parac));
    	       	projectiondividedcloud->points[p0].x = dividedcloud->points[p0].x + (paraa * parat);
    	       	projectiondividedcloud->points[p0].y = dividedcloud->points[p0].y + (parab * parat);
    	       	projectiondividedcloud->points[p0].z = dividedcloud->points[p0].z + (parac * parat);
    	       	
    	       	paraab = sqrt(paraa * paraa + parab * parab);
    	       	pointix = (paraa / paraab) * densifysmoothtrajectoryfiltered->points[indextrajectoryi].x + (parab / paraab) * densifysmoothtrajectoryfiltered->points[indextrajectoryi].y;
    	       	pointiy = (paraa / paraab) * densifysmoothtrajectoryfiltered->points[indextrajectoryi].y - (parab / paraab) * densifysmoothtrajectoryfiltered->points[indextrajectoryi].x;
    	       	fprojectiondividedcloud->points[p0].x = (paraa / paraab) * projectiondividedcloud->points[p0].x + (parab / paraab) * projectiondividedcloud->points[p0].y - pointix;
    	       	fprojectiondividedcloud->points[p0].y = (paraa / paraab) * projectiondividedcloud->points[p0].y - (parab / paraab) * projectiondividedcloud->points[p0].x - pointiy;;
    	       	fprojectiondividedcloud->points[p0].z = projectiondividedcloud->points[p0].z;
    	       }
    	       pcl::io::savePCDFileBinary(savepath + "projectiondividedcloud/Order" + indextrajectory + "X" + indextrajectoryx + "Y" + indextrajectoryy + "Z" + indextrajectoryz +".pcd", *projectiondividedcloud);
    	       pcl::io::savePCDFileBinary(savepath + "fprojectiondividedcloud/Order" + indextrajectory + "X" + indextrajectoryx + "Y" + indextrajectoryy + "Z" + indextrajectoryz +".pcd", *fprojectiondividedcloud);
    	       
    	       pcl::PointCloud<pcl::PointXYZI>::Ptr fprojectiondividedcloudfiltered(new pcl::PointCloud<pcl::PointXYZI>);
    		pcl::PointCloud<pcl::PointXYZI>::Ptr fprojectiondividedclouddensity(new pcl::PointCloud<pcl::PointXYZI>);
    		pcl::PointCloud<pcl::PointXYZI>::Ptr fprojectiondividedcloudsmooth(new pcl::PointCloud<pcl::PointXYZI>);
    		pcl::PointCloud<pcl::PointXYZI>::Ptr fprojectiondividedcloudsmoothfiltered(new pcl::PointCloud<pcl::PointXYZI>);
    		pcl::PointCloud<pcl::PointXYZI>::Ptr crosssection0(new pcl::PointCloud<pcl::PointXYZI>);
    		pcl::PointCloud<pcl::PointXYZI>::Ptr crosssection(new pcl::PointCloud<pcl::PointXYZI>);
    		
    		indexsection = std::to_string(CrosssectionNumber);
    	
    		DeleteRepeatedPoint(fprojectiondividedcloud, fprojectiondividedcloudfiltered, 5, 0.01);
    	
    		CalculateDensity(fprojectiondividedcloudfiltered, fprojectiondividedclouddensity, intensitynum);
    	
    		SmoothCloud(fprojectiondividedclouddensity, fprojectiondividedcloudsmooth, smNUM);
    	
    		DeleteRepeatedPoint(fprojectiondividedcloudsmooth, fprojectiondividedcloudsmoothfiltered, 5, 0.01);
    	
    		SortCrossSection(fprojectiondividedcloudsmoothfiltered, crosssection0);
    	
    		DensifyCrossSection(crosssection0, crosssection, denrate, zmax);
    	
    		pcl::io::savePCDFileBinary(savepath + "crosssection/Order" + indexsection + ".pcd", *crosssection);
    		CrosssectionNumber = CrosssectionNumber + 1;
    	       
    	       dividedcloudsize = GlobalMapzoomsize;
    	       dividedcloud->resize(dividedcloudsize);
    	       i0 = 0;
    	       dividelenght = 0;
    	       }
    	       }
    	    
    }
    std::cout << "GlobalMapzoomsmooth DIVISION FINISHED" << std::endl;
    
   
    std::cout << "PCD ZOOM ELEVATION FINISHED" << std::endl;
    
    
    
    return(0);
}
