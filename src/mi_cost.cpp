#include "mi_cost.h"
#include <boost/unordered_map.hpp>


long getKeyForPoint (float x, float y, float z, float resolution) {
    
    long key = 0;
    
    int keyx = x / resolution;
    int keyy = y / resolution;
    int keyz = z / resolution;
    
    if (x < 0)
        keyx -= 1;
    
    if (y < 0)
        keyy -= 1;
    
    if (z < 0)
        keyz -= 1;
    
    if (keyx > P1)
    	std::cerr << "Warning: Max limit for leaves along x" << std::endl;

    key += keyx * P1;
    key += keyy;
    key *= P2;
    key += keyz;
    
    return key;
}

void createMap (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
                pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
                pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan,
                float& boundsMinX, float& boundsMaxX,
                float& boundsMinY, float& boundsMaxY,
                float& boundsMinZ, float& boundsMaxZ,
                boost::unordered::unordered_map<long, std::vector<int> >& leafMap,
                float resolution) {
    
    
    double minX = std::numeric_limits<float>::max (), minY = std::numeric_limits<float>::max (), minZ = std::numeric_limits<float>::max ();
    double maxX = -std::numeric_limits<float>::max(), maxY = -std::numeric_limits<float>::max(), maxZ = -std::numeric_limits<float>::max();
    
    for (int i = 0; i < A->size(); i++) {
        pcl::PointXYZRGBA temp (A->points[i]);
        
        if (!pcl::isFinite (temp)) //Check to make sure transform didn't make point not finite
            continue;
        
        if (temp.x < minX)
            minX = temp.x;
        if (temp.y < minY)
            minY = temp.y;
        if (temp.z < minZ)
            minZ = temp.z;
        if (temp.x > maxX)
            maxX = temp.x;
        if (temp.y > maxY)
            maxY = temp.y;
        if (temp.z > maxZ)
            maxZ = temp.z;
        
        pcl::PointXYZRGBL p;
        p.x = temp.x;
        p.y = temp.y;
        p.z = temp.z;
        p.r = temp.r;
        p.g = temp.g;
        p.b = temp.b;
        p.label = 0;
        
        combinedScan->push_back(p);
    }
    
    double minX_A = minX, minY_A = minY, minZ_A = minZ;
    double maxX_A = maxX, maxY_A = maxY, maxZ_A = maxZ;
    
    minX = std::numeric_limits<float>::max (), minY = std::numeric_limits<float>::max (), minZ = std::numeric_limits<float>::max ();
    maxX = -std::numeric_limits<float>::max(), maxY = -std::numeric_limits<float>::max(), maxZ = -std::numeric_limits<float>::max();
    
    for (int i = 0; i < B->size(); i++) {
        
        pcl::PointXYZRGBA temp (B->points[i]);
        
        if (!pcl::isFinite (temp)) //Check to make sure transform didn't make point not finite
            continue;
        
        if (temp.x < minX)
            minX = temp.x;
        if (temp.y < minY)
            minY = temp.y;
        if (temp.z < minZ)
            minZ = temp.z;
        if (temp.x > maxX)
            maxX = temp.x;
        if (temp.y > maxY)
            maxY = temp.y;
        if (temp.z > maxZ)
            maxZ = temp.z;
        
        pcl::PointXYZRGBL p;
        p.x = temp.x;
        p.y = temp.y;
        p.z = temp.z;
        p.r = temp.r;
        p.g = temp.g;
        p.b = temp.b;
        p.label = 1;
        
        combinedScan->push_back(p);
    }
    
    // B min/max
    double minX_B = minX, minY_B = minY, minZ_B = minZ;
    double maxX_B = maxX, maxY_B = maxY, maxZ_B = maxZ;
    
    minX = std::max(minX_A, minX_B);
    maxX = std::min(maxX_A, maxX_B);
    minY = std::max(minY_A, minY_B);
    maxY = std::min(maxY_A, maxY_B);
    minZ = std::max(minZ_A, minZ_B);
    maxZ = std::min(maxZ_A, maxZ_B);
    
    boundsMinX = resolution * floor(minX / resolution) ;
    boundsMaxX = resolution * ceil(maxX / resolution) ;
    
    boundsMinY = resolution * floor(minY / resolution) ;
    boundsMaxY = resolution * ceil(maxY / resolution) ;
    
    boundsMinZ = resolution * floor(minZ / resolution) ;
    boundsMaxZ = resolution * ceil(maxZ / resolution) ;
    
    // create map for only leaves that need to be used
    // Create bounding box from minb maxb
    
    for (int i = 0; i < combinedScan->size(); i++) {
        
        pcl::PointXYZRGBL p (combinedScan->points[i]);
        
        if (boundsMinX < p.x && boundsMaxX > p.x &&
            boundsMinY < p.y && boundsMaxY > p.y &&
            boundsMinZ < p.z && boundsMaxZ > p.z) {
            
            long key = getKeyForPoint(p.x, p.y, p.z, resolution);
            
            if (leafMap.find(key) != leafMap.end()) {
                leafMap[key].push_back(i);
            } else {
                std::vector<int> indices;
                indices.push_back(i);
                leafMap.insert(std::pair<long, std::vector<int> > (key, indices));
            }
        }
        
    }
    
}


double
calculateMIFromMap (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA,
                    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB,
                    float resolution) {
    
    long maxPointsInVoxelForEachScan = 10000;
    float minXo, minYo, minZo;
    float maxXo, maxYo, maxZo;
    
    std::vector<std::pair<double, double> > data;
    pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBL> > (new pcl::PointCloud<pcl::PointXYZRGBL> ());
    
    boost::unordered::unordered_map<long, std::vector<int> > leafMap;
    createMap(scanA, scanB, combinedScan, minXo, maxXo, minYo, maxYo, minZo, maxZo, leafMap, resolution);

    /*
    std::cout << "XMin: " << minXo << std::endl;
    std::cout << "XMax: " << maxXo << std::endl;
    std::cout << "YMin: " << minYo << std::endl;
    std::cout << "YMax: " << maxYo << std::endl;
    std::cout << "ZMin: " << minZo << std::endl;
    std::cout << "ZMax: " << maxZo << std::endl;
    */
    
    boost::unordered::unordered_map<int, int> umapA;
    boost::unordered::unordered_map<int, int> umapB;
    boost::unordered::unordered_map<long, int> umapAB;
    
    boost::unordered::unordered_map<long, std::vector<int> >::iterator iterator = leafMap.begin();
    
    while (iterator != leafMap.end()) {
        
//        std::cout << "Key: " << iterator->first << std::endl;
        std::vector<int> indexVector = iterator->second;
        
        int pointsFromA(0), pointsFromB(0);
        
        pcl::PointXYZRGBL p;
        for (int i = 0; i < indexVector.size(); i++) {
            
            int scanIndex = indexVector[i];
            p = combinedScan->at(scanIndex);

//            std::cout << "x: " << p.x << std::endl;
//            std::cout << "y: " << p.y << std::endl;
//            std::cout << "z: " << p.z << std::endl;
            
            if (p.label == 0)
                pointsFromA++;
            else
                pointsFromB++;
        }
        
        
        if (umapA.find(pointsFromA) == umapA.end())
            umapA.insert(std::pair<int, int> (pointsFromA, 1));
        else
            umapA[pointsFromA]++;
        
        if (umapB.find(pointsFromB) == umapB.end())
            umapB.insert(std::pair<int, int> (pointsFromB, 1));
        else
            umapB[pointsFromB]++;
        
        long key = pointsFromA * maxPointsInVoxelForEachScan + pointsFromB;
        
        if (umapAB.find(key) == umapAB.end()) {
            umapAB.insert(std::pair<long, int> (key, 1));
        } else {
            umapAB[key]++;
        }
        
        data.push_back(std::pair<int, int> (pointsFromA, pointsFromB));
        iterator++;
        
    }
    
    unsigned long total = (maxXo - minXo) * (maxYo - minYo) * (maxZo - minZo) / (resolution * resolution * resolution);
    
    unsigned long leavesNotInOctree = total-leafMap.size();
    
    if (umapA.find(0) == umapA.end()) {
        umapA.insert(std::pair<int, int> (0, leavesNotInOctree));
    } else {
        umapA[0] += leavesNotInOctree;
    }
    
    if (umapB.find(0) == umapB.end()) {
        umapB.insert(std::pair<int, int> (0, leavesNotInOctree));
    } else {
        umapB[0] += leavesNotInOctree;
    }
    
    umapAB.insert(std::pair<long, int> (0, leavesNotInOctree));
    
    
    double jointEntropy(0), marginalEntropyA(0), marginalEntropyB(0);
    
    boost::unordered::unordered_map<int, int>::iterator itr = umapA.begin();
    while (itr != umapA.end()) {
        
        double probability = ( (double) itr->second) / total;

        marginalEntropyA += -probability * log(probability);
        
        itr++;
    }
    
    itr = umapB.begin();
    while (itr != umapB.end()) {
        
        double probability = ( (double) itr->second) / total;
        marginalEntropyB += -probability * log(probability);
        
        itr++;
    }
    
    boost::unordered::unordered_map<long, int>::iterator itrAB = umapAB.begin();
    while (itrAB != umapAB.end()) {
        
        double probability = ( (double) itrAB->second) / total;
        jointEntropy += -probability * log(probability);
        
        itrAB++;
    }
    
    double mi = marginalEntropyA + marginalEntropyB - jointEntropy;
    return mi;
}

void getCombinedScan (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr A,
                      pcl::PointCloud<pcl::PointXYZRGBA>::Ptr B,
                      pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan,
                      pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::Ptr octree,
                      int& minXo, int& minYo,
                      int& minZo, int& maxXo,
                      int& maxYo, int& maxZo,
                      float resolution) {
    
    octree->setResolution(resolution);
    double minX = std::numeric_limits<float>::max (), minY = std::numeric_limits<float>::max (), minZ = std::numeric_limits<float>::max ();
    double maxX = -std::numeric_limits<float>::max(), maxY = -std::numeric_limits<float>::max(), maxZ = -std::numeric_limits<float>::max();
    
    for (int i = 0; i < A->size(); i++) {
        pcl::PointXYZRGBA temp (A->points[i]);
        
        if (!pcl::isFinite (temp)) //Check to make sure transform didn't make point not finite
            continue;
        
        if (temp.x < minX)
            minX = temp.x;
        if (temp.y < minY)
            minY = temp.y;
        if (temp.z < minZ)
            minZ = temp.z;
        if (temp.x > maxX)
            maxX = temp.x;
        if (temp.y > maxY)
            maxY = temp.y;
        if (temp.z > maxZ)
            maxZ = temp.z;
        
        pcl::PointXYZRGBL p;
        p.x = temp.x;
        p.y = temp.y;
        p.z = temp.z;
        p.r = temp.r;
        p.g = temp.g;
        p.b = temp.b;
        p.label = 0;
        
        combinedScan->push_back(p);
    }
    
    double minX_A = minX, minY_A = minY, minZ_A = minZ;
    double maxX_A = maxX, maxY_A = maxY, maxZ_A = maxZ;
    
    minX = std::numeric_limits<float>::max (), minY = std::numeric_limits<float>::max (), minZ = std::numeric_limits<float>::max ();
    maxX = -std::numeric_limits<float>::max(), maxY = -std::numeric_limits<float>::max(), maxZ = -std::numeric_limits<float>::max();
    
    for (int i = 0; i < B->size(); i++) {
        
        pcl::PointXYZRGBA temp (B->points[i]);
        
        if (!pcl::isFinite (temp)) //Check to make sure transform didn't make point not finite
            continue;
        
        if (temp.x < minX)
            minX = temp.x;
        if (temp.y < minY)
            minY = temp.y;
        if (temp.z < minZ)
            minZ = temp.z;
        if (temp.x > maxX)
            maxX = temp.x;
        if (temp.y > maxY)
            maxY = temp.y;
        if (temp.z > maxZ)
            maxZ = temp.z;
        
        pcl::PointXYZRGBL p;
        p.x = temp.x;
        p.y = temp.y;
        p.z = temp.z;
        p.r = temp.r;
        p.g = temp.g;
        p.b = temp.b;
        p.label = 1;
        
        combinedScan->push_back(p);
    }
    
    double minX_B = minX, minY_B = minY, minZ_B = minZ;
    double maxX_B = maxX, maxY_B = maxY, maxZ_B = maxZ;
    
    minX = std::min(minX_A, minX_B);
    maxX = std::max(maxX_A, maxX_B);
    minY = std::min(minY_A, minY_B);
    maxY = std::max(maxY_A, maxY_B);
    minZ = std::min(minZ_A, minZ_B);
    maxZ = std::max(maxZ_A, maxZ_B);
    
    int minXb = minX/resolution;
    minXb = resolution * minXb;
    if (minXb < 0)
        minXb -= resolution;
    
    int minYb = minY/resolution;
    minYb = resolution * minYb;
    if  (minYb < 0)
        minYb -= resolution;
    
    int minZb = minZ/resolution;
    minZb = resolution * minZb;
    if (minZb < 0)
        minZb -= resolution;
    
    int maxXb = maxX/resolution;
    maxXb = resolution * maxXb;
    if (maxXb > 0)
        maxXb += resolution;
    
    int maxYb = maxY / resolution;
    maxYb = resolution * maxYb;
    if (maxYb > 0)
        maxYb += resolution;
    
    int  maxZb = maxZ / resolution;
    maxZb = resolution*maxZb;
    if (maxZb > 0)
        maxZb += resolution;
    
    octree->defineBoundingBox (minXb, minYb, minZb, maxXb, maxYb, maxZb);
    octree->setInputCloud(combinedScan);
    octree->addPointsFromInputCloud();
    //octree->getBoundingBox(minX, minY, minZ, maxX, maxY, maxZ);
    
    cout << "MinX: " << minX << std::endl;
    cout << "MinY: " << minY << std::endl;
    cout << "MinZ: " << minZ << std::endl;
    cout << "MaxX: " << maxX << std::endl;
    cout << "MaxY: " << maxY << std::endl;
    cout << "MaxZ: " << maxZ << std::endl;
    cout << std::endl;
    cout << "MinX: " << minXb << std::endl;
    cout << "MinY: " << minYb << std::endl;
    cout << "MinZ: " << minZb << std::endl;
    cout << "MaxX: " << maxXb << std::endl;
    cout << "MaxY: " << maxYb << std::endl;
    cout << "MaxZ: " << maxZb << std::endl;
    
    minXo = floor(std::max(minX_A, minX_B));
    minYo = floor(std::max(minY_A, minY_B));
    minZo = floor(std::max(minZ_A, minZ_B));
    
    maxXo = ceil(std::min(maxX_A, maxX_B));
    maxYo = ceil(std::min(maxY_A, maxY_B));
    maxZo = ceil(std::min(maxZ_A, maxZ_B));
    
    //	cout << "MinXo: " << minXo << std::endl;
    //	cout << "MinYo: " << minYo << std::endl;
    //	cout << "MinZo: " << minZo << std::endl;
    //	cout << "MaxXo: " << maxXo << std::endl;
    //	cout << "MaxYo: " << maxYo << std::endl;
    //	cout << "MaxZo: " << maxZo << std::endl;
    
}


double calculateMI (pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanA,
                    pcl::PointCloud<pcl::PointXYZRGBA>::Ptr scanB,
                    float resolution) {
    
    long maxPointsInVoxelForEachScan = 10000;
    int minXo, minYo, minZo;
    int maxXo, maxYo, maxZo;
    
    std::vector<std::pair<double, double> > data;
    
    pcl::PointCloud<pcl::PointXYZRGBL>::Ptr combinedScan = boost::shared_ptr <pcl::PointCloud<pcl::PointXYZRGBL> > (new pcl::PointCloud<pcl::PointXYZRGBL> ());
    pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::Ptr octree = boost::shared_ptr <pcl::octree::OctreePointCloud<pcl::PointXYZRGBL> > (new pcl::octree::OctreePointCloud<pcl::PointXYZRGBL> (resolution));
    
    getCombinedScan(scanA, scanB, combinedScan, octree, minXo, minYo, minZo, maxXo, maxYo, maxZo, resolution);
    
    typename pcl::octree::OctreePointCloud<pcl::PointXYZRGBL>::LeafNodeIterator iterator = octree->leaf_begin();
    std::vector<int> indexVector;
    
    boost::unordered::unordered_map<int, int> umapA;
    boost::unordered::unordered_map<int, int> umapB;
    boost::unordered::unordered_map<long, int> umapAB;
    
    int overlappingLeavesInOctree = 0;
    while (iterator != octree->leaf_end()) {
        
        typename pcl::octree::OctreePointCloud<pcl::PointXYZRGBA>::LeafNode *leafNode = (typename pcl::octree::OctreePointCloud<pcl::PointXYZRGBA>::LeafNode*) *iterator;
        
        indexVector.clear();
        leafNode->getContainer().getPointIndices(indexVector);
        
        int pointsFromA(0), pointsFromB(0);
        
        pcl::PointXYZRGBL p;
        for (int i = 0; i < indexVector.size(); i++) {
            
            int scanIndex = indexVector[i];
            p = combinedScan->at(scanIndex);
            
            if (p.label == 0)
                pointsFromA++;
            else
                pointsFromB++;
        }
        
        int xleafMin, xLeafMax, yLeafMin, yLeafMax, zLeafMin, zLeafMax;
        xleafMin = resolution * p.x/resolution;
        yLeafMin = resolution * p.y/resolution;
        zLeafMin = resolution * p.z/resolution;
        
        if (p.x < 0)
            xleafMin -= resolution;
        if (p.y < 0)
            yLeafMin -= resolution;
        if (p.z < 0)
            zLeafMin -= resolution;;
        
        xLeafMax = xleafMin + resolution;
        yLeafMax = yLeafMin + resolution;
        zLeafMax = zLeafMin + resolution;
        
        // check leaf in overlapping region
        
        if (	xleafMin > minXo && xLeafMax < maxXo &&
            yLeafMin > minYo && yLeafMax < maxYo &&
            zLeafMin > minZo && zLeafMax < maxZo) {
            
            // leaf in overlapping region
            overlappingLeavesInOctree++;
            
            if (umapA.find(pointsFromA) == umapA.end())
                umapA.insert(std::pair<int, int> (pointsFromA, 1));
            else
                umapA[pointsFromA]++;
            
            if (umapB.find(pointsFromB) == umapB.end())
                umapB.insert(std::pair<int, int> (pointsFromB, 1));
            else
                umapB[pointsFromB]++;
            
            long key = pointsFromA * maxPointsInVoxelForEachScan + pointsFromB;
            
            if (umapAB.find(key) == umapAB.end()) {
                umapAB.insert(std::pair<long, int> (key, 1));
            } else {
                umapAB[key]++;
            }
            
            data.push_back(std::pair<int, int> (pointsFromA, pointsFromB));
            
        }
        
        iterator++;
    }
    
    
    long total = (maxXo - minXo) * (maxYo - minYo) * (maxZo - minZo);
    //	cout << "Total voxels including empty space in overlapping cuboid: " << total << std::endl;
    
    long leavesNotInOctree = total-overlappingLeavesInOctree;
    
    if (umapA.find(0) == umapA.end()) {
        umapA.insert(std::pair<int, int> (0, leavesNotInOctree));
    } else {
        umapA[0] += leavesNotInOctree;
    }
    
    if (umapB.find(0) == umapB.end()) {
        umapB.insert(std::pair<int, int> (0, leavesNotInOctree));
    } else {
        umapB[0] += leavesNotInOctree;
    }
    
    umapAB.insert(std::pair<long, int> (0, leavesNotInOctree));
    
    double jointEntropy(0), marginalEntropyA(0), marginalEntropyB(0);
    
    boost::unordered::unordered_map<int, int>::iterator itr = umapA.begin();
    while (itr != umapA.end()) {
        
        double probability = ( (double) itr->second) / total;
        marginalEntropyA += -probability * log(probability);
        
        itr++;
    }
    
    itr = umapB.begin();
    while (itr != umapB.end()) {
        
        double probability = ( (double) itr->second) / total;
        marginalEntropyB += -probability * log(probability);
        
        itr++;
    }
    
    boost::unordered::unordered_map<long, int>::iterator itrAB = umapAB.begin();
    while (itrAB != umapAB.end()) {
        
        double probability = ( (double) itrAB->second) / total;
        jointEntropy += -probability * log(probability);
        
        itrAB++;
    }
    
    //	cout << "Marginal Entropy A: " << marginalEntropyA << std::endl;
    //	cout << "Marginal Entropy B: " << marginalEntropyB << std::endl;
    //	cout << "Joint Entropy: " << jointEntropy << std::endl;
    
    double mi = marginalEntropyA + marginalEntropyB - jointEntropy;    
    return mi;
}
