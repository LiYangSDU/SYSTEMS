#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <string>
#include <string.h>
#include <stdlib.h>
#include <dirent.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <algorithm>
#include <math.h>
#include <limits.h>
#include <list>
#include <queue>

using namespace std;

#define MIN(a,b) ((a)<(b)?(a):(b))
#define MAX(a,b) ((a)>(b)?(a):(b))
#define DIR_LENGTH 200
#define PARTITION_CUTOFF 0.382
#define EDGE_CUTOFF 0
#define EDGE_CUTOFF2 0
#define WEIGHT_CUTOFF 0.5
#define BASE_WEIGHT 0.5
#define JACCARD_CUTOFF 0.618
#define MAX_NODE 3
#define SUM_WEIGHT_CUTOFF 0.2
#define TOP 10

typedef struct
{
	vector<int> node;
	vector<pair<int, int> > edge;
	vector<int> occurr;
	int freq;
	float weight;
} mod;

typedef pair<int, int> pairInt;
typedef vector<float> vecFloat;
typedef vector<int> vecInt;
typedef vector<mod> vecMod;

typedef struct
{
	map<int, int> nodeList;
	vector<vecFloat> edgeList;
	map<int, list<int> > neighborList;
	vector<list<int> > compList;
} graph;

// Function: partition
// Used to be called by quickSort

template <class T>
int partition(T a[], int start, int stop, int id[])
{
        int temp_id, up = start, down = stop - 1;
        T temp_value, part = a[stop];
        if (stop <= start) return start;

        while (true)
        {
                while (a[up ] <part) up ++;
                while (part < a[down] && (up < down)) down --;

                if (up >= down) break;

                temp_value = a[up];  a[up] = a[down];
		a[down] = temp_value;
                temp_id = id[up]; id[up] = id[down];
		id[down] = temp_id;

                up ++; down --;
        }

        temp_value = a[up]; a[up] = a[stop]; a[stop] = temp_value;
        temp_id = id[up]; id[up] = id[stop]; id[stop] = temp_id;
        return up;
}

// Function: quickSort
// Used to sort the values in array

template <class T>
void quickSort(T a[], int start, int stop, int id[])
{
        int i;
        if (stop <= start) return;

        i = partition(a, start, stop, id);
        quickSort(a, start, i - 1, id);
        quickSort(a, i + 1, stop, id);
}

// Function: getDim
// Used to get dimensions of data

void getDim(string& origPath, int& rowNum, 
		int& colNum)
{
	ifstream inFileHd(origPath.c_str());

	if (!inFileHd)
	{
		cerr << "Error: Cannot open file: " 
			<< origPath << ".\n";
		exit(1);
	}

	string strLn, strSep; colNum = 0;
	getline(inFileHd, strLn);

	if (strLn.empty())
	{
		cerr << "Error: The first line of file: " 
			<< origPath << " is empty.\n";
		exit(1);
	}

	istringstream sepHd(strLn);

	while (sepHd >> strSep)
	{
		colNum ++;
//		cout << colNum << "\n";
	}

//	cout << colNum << "\n";

	rowNum = 1;

	while (getline(inFileHd, strLn))
	{
		if (strLn.empty())
		{
			continue;
		}

		rowNum ++;
	}

	cout << "The size of the original data is: " 
		<< rowNum << " X " << colNum - 2 << ".\n";
}

// Function: readData
// Used to read original data into program

void readData(string& origPath, float** origMat, 
		int& rowNum, int& colNum, 
		graph& sumGraph, int& freqCutoff)
{
	ifstream inFileHd(origPath.c_str());
	sumGraph.edgeList.reserve(rowNum);

	for (int i = 0; i < rowNum; i ++)
	{
		for (int j = 0; j < colNum; j ++)
		{
			inFileHd >> origMat[i][j];
		}
	}

	inFileHd.close();
	
	cout << "Finished reading file: " << 
		origPath << " into program.\n";

	for (int i = 0; i < rowNum; i ++)
	{
		vecFloat vecEdge; vecEdge.reserve(3);
		vecEdge.push_back(origMat[i][0]);
		vecEdge.push_back(origMat[i][1]);

		sumGraph.nodeList[int(origMat[i][0])] = 0;
		sumGraph.nodeList[int(origMat[i][1])] = 0;

		float weight = 0; int occurr = 0;

		for (int j = 2; j < colNum; j ++)
		{
			if (origMat[i][j] > 0)
			{
				occurr ++;
				weight += origMat[i][j];
			}
		}

		if (occurr < freqCutoff)
		{
			continue;
		}

		vecEdge.push_back(weight);

		sumGraph.edgeList.push_back(vecEdge);
	}

	cout << "Finished constructing the summary graph.\n" 
		<< "There're altogether " << sumGraph.edgeList.size() 
		<< " edges.\n";
}

// Function: normGraph
// Used to normalize the ege weights

void normGraph(graph& G)
{
	for (int i = 0; i < G.edgeList.size(); i ++)
	{
//		G.edgeList[i][2] = G.edgeList[i][2] * 
//			G.edgeList[i][2];
//		G.edgeList[i][2] = pow(2, G.edgeList[i][2]);
	}

	double maxVal = 0;

	for (int i = 0; i < G.edgeList.size(); i ++)
	{
		if (maxVal < G.edgeList[i][2])
		{
			maxVal = G.edgeList[i][2];
		}
	}

	for (int i = 0; i < G.edgeList.size(); i ++)
	{
		G.edgeList[i][2] /= maxVal;
	}

	cout << "Finished normalizing the edge weights.\n";
}

// Function: getNeighbor
// Used to get hte neighbor list of nodes
// in graph

void getNeighbor(graph& G, int startId)
{
	for (map<int, list<int> >::iterator it = G.neighborList.begin(); 
		it != G.neighborList.end(); it ++)
	{
		(it -> second).clear();
	}

	for (int i = startId; i < G.edgeList.size(); i ++)
	{
		if (G.neighborList.count(G.edgeList[i][0]) == 0)
		{
			list<int> oneList;
			oneList.push_back(G.edgeList[i][1]);
			G.neighborList[G.edgeList[i][0]] = oneList;
		}
		else
		{
			G.neighborList[G.edgeList[i][0]].push_back
				(G.edgeList[i][1]);
		}

		if (G.neighborList.count(G.edgeList[i][1]) == 0)
		{
			list<int> oneList;
			oneList.push_back(G.edgeList[i][0]);
			G.neighborList[G.edgeList[i][1]] = oneList;
		}
		else
		{
			G.neighborList[G.edgeList[i][1]].push_back
				(G.edgeList[i][0]);
		}

//		G.nodeList[G.edgeList[i][0]] = 1;
//		G.nodeList[G.edgeList[i][1]] = 1;

	}

	for (map<int, int>::iterator it = G.nodeList.begin(); it != 
		G.nodeList.end(); it ++)
	{
		if (G.neighborList.count(it -> first) > 0)
		{
			continue;
		}

		list<int> oneList;
		G.neighborList[it -> first] = oneList;
	}

//	for (map<int, list<int> >::iterator it = G.neighborList.begin(); 
//		it != G.neighborList.end(); it ++)
//	{
//		cout << it -> first << "\t";
//
//		for (list<int>::iterator iit = (it -> second).begin(); 
//			iit != (it -> second).end(); iit ++)
//		{
//			cout << *iit << " ";
//		}
//
//		cout << "\n";
//	}

}

// Function: rankEdge
// Used to rank the edges

void rankEdge(graph& G)
{
	float* weightArray; int *idArray;

	weightArray = new float [G.edgeList.size()];
	idArray = new int [G.edgeList.size()];

	for (int i = 0; i < G.edgeList.size(); i ++)
	{
		weightArray[i] = G.edgeList[i][2];
		idArray[i] = i;
	}

	quickSort(weightArray, 0, G.edgeList.size() - 1, 
			idArray);

	{
		vector<vecFloat> vecSwap;

		for (int i = 0; i < G.edgeList.size(); i ++)
		{
			vecSwap.push_back(G.edgeList[idArray[i]]);
		}

		G.edgeList.swap(vecSwap);
	}


	cout << "Finished ranking the edges according to weights.\n";
}

// Function: findRot
// Used to get the root to perform BFS

int findRoot(map<int, int>& mapFlag)
{
	int root = -1;

	for (map<int, int>::iterator it = mapFlag.begin(); 
		it != mapFlag.end(); it ++)
	{
		if ((it -> second) == 0)
		{
			root = (it -> first);
			break;
		}
	}

	return root;
}

// Function: addNeighbor
// Used to add neighbors of one node to component

void addNeighbor(graph& G, int& thisNode, map<int, int>& 
		mapFlag, list<int>& nodeInComp, queue<int>& 
		nodeBatch)
{
	for (list<int>::iterator it = G.neighborList[thisNode].begin(); 
		it != G.neighborList[thisNode].end(); it ++)
	{
//		cout << *it << "\tFlag\t" << mapFlag[*it] << "\n";

		if (mapFlag[*it] == 1)
		{
			continue;
		}

//		cout << *it << "\tFlag\t" << mapFlag[*it] << "\n";
//		cout << "2 flag: " << mapFlag[2] << "\n";

		nodeBatch.push(*it);
//		nodeInComp.push_back(*it);
		mapFlag[*it] = 1;
	}
}

// Function: getComp
// Used to get the components of one graph
// with part of the edges retained

void getComp(graph& G)
{
	G.compList.clear();
	map<int, int> mapFlag;

	for (map<int, list<int> >::iterator it = G.neighborList.begin(); 
		it != G.neighborList.end(); it ++)
	{
//		cout << it -> first << "\n";

		mapFlag[it -> first] = 0;
	}

	while (1)
	{
		int root = findRoot(mapFlag);

//		cout << " Root: " << root << "\n";

		if (root == -1)
		{
			break;
		}

		list<int> nodeInComp;
		queue<int> nodeBatch;
		mapFlag[root] = 1;
//		nodeInComp.push_back(root);
		nodeBatch.push(root);

//		cout << "Before while: " << root << "\n";

		while (nodeBatch.size() != 0)
		{
			int thisNode = nodeBatch.front();
			
//			if (thisNode == 2)
//			{
//				cout << thisNode << "\t";
//				nodeBatch.pop();
//				cout <<nodeBatch.front() << "\n";
//				exit(0);
//			}

//			cout << "In while: " << thisNode << "\n";
//			cout << "Queue size 1: " << nodeBatch.size() << "\n";

			nodeBatch.pop();

//			cout << nodeBatch.size() << "\n";

			nodeInComp.push_back(thisNode);
//			mapFlag[thisNode] = 1;
//			cout << thisNode << "\t";

//			cout << "Queue size 2: " << nodeBatch.size() << "\n";

			addNeighbor(G, thisNode, mapFlag, 
				nodeInComp, nodeBatch);
			
//			cout << "Queue size 3: " << nodeBatch.size() << "\n";
		}

//		cout << "\n";

//		cout << "After while: " << root << "\n";

		G.compList.push_back(nodeInComp);
	}

//	for (int i = 0; i < G.compList.size(); i ++)
//	{
//		for (list<int>::iterator it = G.compList[i].begin(); 
//				it != G.compList[i].end(); it ++)
//		{
//			cout << *it << " ";
//		}
//
//		cout << "\n";
//
//		cout << G.compList[i].size() << "\n";
//	}

//	exit(0);
}

// Function: getMaxCompSize
// Used to get the maximum component size

int getMaxCompSize(graph& G)
{
	int maxSize = 0;

	for (int i = 0; i < G.compList.size(); i ++)
	{
		if (G.compList[i].size() > maxSize)
		{
			maxSize = G.compList[i].size();
		}
	}

	return maxSize;
}

// Function: optimalComp
// Used to get the components of a graph

void optimalComp(graph& G)
{
	rankEdge(G);

//	for (int i = 0; i < G.edgeList.size(); i ++)
//	{
//		cout << "optimalComp: " << i << "\n";
//
//		getNeighbor(G, i);
//
//		cout <<"Before getComp: " << i << "\n";
//
//		getComp(G);
//		
//		cout <<"After getComp: " << i << "\n";
//		
//		if (getMaxCompSize(G) <= MAX_NODE)
//		{
//			break;
//		}
//	}

	int maxId = G.edgeList.size() - 1;
	int minId = 0;
	int optimalId = maxId / 2;

	while (minId <= maxId)
	{
		optimalId = (maxId + minId) / 2;
//		getNeighbor(G, 0);
		getNeighbor(G, optimalId);
		getComp(G);
		int maxCompSize = getMaxCompSize(G);

		if (maxCompSize == MAX_NODE)
		{
			break;
		}
		else if (maxCompSize < MAX_NODE)
		{
			maxId = optimalId - 1;
		}
		else
		{
			minId = optimalId + 1;
		}
	}

	for (int i = 0; i < G.compList.size(); i ++)
	{
//		for (list<int>::iterator it = G.compList[i].begin(); 
//				it != G.compList[i].end(); it ++)
//		{
//			cout << *it << " ";
//		}
//
//		cout << "\n";
//
		cout << G.compList[i].size() << "\n";
	}

	exit(0);

	cout << "Finished partitioning the node sets.\n" << 
		"There're altogether " << G.compList.size() << 
		" connected components.\n";
}

// Function: Array2Vec
// Used to transform array into vector

template <class T>

void Array2Vec(T* lnArray, int start, int end, 
		vector<T>& lnVec)
{
	lnVec.clear();
	lnVec.reserve(end - start + 1);

	for (int i = start; i <= end; i ++)
	{
		lnVec.push_back(lnArray[i]);
	}
}

// Function: getFreq
// Used to get frequency of an edge

template <class T>

int getFreq(T* lnArray, int begin, int end)
{
	int freq = 0;

	for (int i = begin; i <= end; i ++)
	{
		if (lnArray[i] > 0)
		{
			freq ++;
		}
	}

	return freq;
}

// Function: getLocalMap
// Used to build map for one graph

void getLocalMap(map<pairInt, vecFloat>& localMap, 
		vector<int>& nodeSet, float** 
		origMat, int& rowNum, int& colNum, 
		int& freqCutoff)
{
	map<int, int> nodeMap;

	for (int i = 0; i < nodeSet.size(); i ++)
	{
		nodeMap[nodeSet[i]] = 0;
	}

	for (int i = 0; i < rowNum; i ++)
	{
		if (nodeMap.count(int(origMat[i][0])) 
			== 0 || nodeMap.count
			(int(origMat[i][1])) 
			== 0)
		{
			continue;
		}
		
		if (getFreq(origMat[i], 2, colNum - 1) 
				< freqCutoff)
		{
			continue;
		}

		vector<float> weightVec;
		weightVec.reserve(colNum - 2);
		Array2Vec(origMat[i], 2, colNum - 1, 
				weightVec);
		localMap[make_pair(int(origMat[i][0]), 
			int(origMat[i][1]))] = 
			weightVec;
	}
}

// Function: getBiMap
// Used to build map for one bipartite graph

//void getBiMap(map<pairInt, vecFloat>& biMap, 
//		vector<int>& nodeSet1, vector<int>& 
//		nodeSet2, float** origMat, int& 
//		rowNum, int& colNum, int& freqCutoff)
//{
//	map<int, int> nodeMap1, nodeMap2;
//
//	for (int i = 0; i < nodeSet1.size(); i ++)
//	{
//		nodeMap1[nodeSet1[i]] = 0;
//	}
//
//	for (int i = 0; i < nodeSet2.size(); i ++)
//	{
//		nodeMap2[nodeSet2[i]] = 0;
//	}
//
//	for (int i = 0; i < rowNum; i ++)
//	{
//		if (getFreq(origMat[i], 2, colNum - 1) 
//				< freqCutoff)
//		{
//			continue;
//		}
//
//		if (nodeMap1.count(int(origMat[i][0])) 
//			== 0 && nodeMap2.count(int(
//			origMat[i][1])) == 0 || nodeMap1.count
//			(int(origMat[i][1])) == 0 && 
//			nodeMap2.count(int(origMat[i][0])) == 0)
//		{
//			continue;
//		}
//
//		vector<float> weightVec;
//		weightVec.reserve(colNum - 2);
//		Array2Vec(origMat[i], 2, colNum - 1, 
//				weightVec);
//		biMap[make_pair(int(origMat[i][0]), 
//			int(origMat[i][1]))] = 
//			weightVec;
//	}
//}

// Function: getNbrList
// Used to get neighbor list

void getNbrList(map<int, vecInt>& nbrListMap, 
		map<pairInt, vecFloat>& localMap)
{
	for (map<pairInt, vecFloat>::iterator it = 
		localMap.begin(); it != localMap.end(); 
		it ++)
	{
		nbrListMap[(it -> first).first].
			push_back((it -> first).second);
		nbrListMap[(it -> first).second].
			push_back((it -> first).first);
	}
}

// Function: getCommonNbr
// Used to get the common neighbors 
// between two nodes

void getCommonNbr(vector<int>& commonNbr, vector<int>& 
		vec1, vector<int>& vec2)
{
	commonNbr.clear(); commonNbr.reserve(
		MIN(vec1.size(), vec2.size()));

	for (int i = 0; i < vec1.size(); i ++)
	{
		if (find(vec2.begin(), 
			vec2.end(), 
			vec1[i]) != 
			vec2.end())
		{
			commonNbr.push_back
			(vec1[i]);
		}
	}


//	for (int i = 0; i < nbrListMap[6].size(); i ++)
//	{
//		cout << nbrListMap[6][i] << " ";
//	}
//	
//	cout << "\n";
//	exit(0);

//	for (int i = 0; i < commonNbr.size(); i ++)
//	{
//		cout << commonNbr[i] << " ";
//	}

//	cout << "\n";

}

// Function: innerProd
// Used to get inner product of three vectors

float innerProd(vector<float>& centerVec, vector<float>& 
		vec1, vector<float>& vec2)
{
	float weight = 0;

	for (int i = 0; i < centerVec.size(); i ++)
	{
		weight += centerVec[i] * vec1[i] * 
			vec2[i];
//		cout << centerVec[i] << "\t" << vec1[i] << "\t" << vec2[i] << "\n";
	}

	return weight;
}

// Function: getWeight
// Used to get weight of an edge

void getWeight(float& weight, pairInt centerPair, 
		vecFloat& centerVec, vector<int>& 
		commonNbr, map<pairInt, vecFloat>& 
		localMap)
{
	weight = 0; float maxWeight = 0;

	for (int i = 0; i < centerVec.size(); i ++)
	{
		weight += centerVec[i];
	}

	for (int i = 0; i < commonNbr.size(); i ++)
	{
		pairInt edge1, edge2;

		if (commonNbr[i] > centerPair.first)
		{
			edge1 = make_pair(commonNbr[i], 
				centerPair.first);
		}
		else
		{
			edge1 = make_pair(centerPair.first, 
					commonNbr[i]);
		}

		if (commonNbr[i] > centerPair.second)
		{
			edge2 = make_pair(commonNbr[i], 
				centerPair.second);
		}
		else
		{
			edge2 = make_pair(centerPair.second, 
					commonNbr[i]);
		}

		weight += innerProd(centerVec, localMap[edge1], 
				localMap[edge2]);
	}
}

// Function: getCompGraph
// Used to get compatible graph

void getCompGraph(vector<vector<float> >& compGraph, 
		map<pairInt, vecFloat>& localMap, 
		float& actDenCutoff)
{
	map<int, vecInt> nbrListMap; float maxWeight = 0;
	getNbrList(nbrListMap, localMap); float minWeight = INT_MAX;
	compGraph.reserve(localMap.size());

	for (map<pairInt, vecFloat>::iterator it = 
		localMap.begin(); it != localMap.end(); 
		it ++)
	{
		vector<int> commonNbr; float weight;
		getCommonNbr(commonNbr, nbrListMap[(it -> first).
		first], nbrListMap[(it -> first).second]);
		getWeight(weight, it -> first, it -> 
			second, commonNbr, localMap);

//		cout << (it -> first).first << " " << 
//			(it -> first).second << "\n";
//		weight = 1;

		if (weight == 0)
		{
			continue;
		}

		vector<float> lnVec; lnVec.push_back((it -> 
					first).first);
		lnVec.push_back((it -> first).second);
		lnVec.push_back(weight);
		compGraph.push_back(lnVec);

		if (maxWeight < weight)
		{
			maxWeight = weight;
		}

		if (minWeight > weight)
		{
			minWeight = weight;
		}
	}

	if (maxWeight == 0)
	{
		return;
	}

	for (int i = 0; i < compGraph.size(); i ++)
	{
		if (maxWeight == minWeight)
		{
			compGraph[i][2] = BASE_WEIGHT;
		}
		else
		{
			compGraph[i][2] = BASE_WEIGHT + 
				(compGraph[i][2] - minWeight) * 
				(1 - BASE_WEIGHT) / 
				(maxWeight - minWeight);
		}

//		cout << compGraph[i][0] << " " << 
//			compGraph[i][1] << " " << 
//			compGraph[i][2] << "\n";
//		compGraph[i][2] = 1;
	}
	
//	cout << "\n\n";

//	float minWeight = INT_MAX;
//
//	float avgWeight = 0; float stdDev = 0;
//
//	for (int i = 0; i < compGraph.size(); i ++)
//	{
//		avgWeight += compGraph[i][2];
//		stdDev += compGraph[i][2] * compGraph[i][2];
//
//		if (compGraph[i][2] == 0)
//		{
//			continue;
//		}
//
//		if (minWeight > compGraph[i][2])
//		{
//			minWeight = compGraph[i][2];
//		}
//	}
//
//	avgWeight /= compGraph.size();
//	stdDev /= compGraph.size();
//	stdDev -= (avgWeight * avgWeight);
//	stdDev = sqrt(stdDev);
//	actDenCutoff *= minWeight;
//	actDenCutoff *= avgWeight;
//	actDenCutoff *= SHRINK;
//	cout << "Density cutoff: " << actDenCutoff << "\n";
//	float weightCutoff = avgWeight + EDGE_CUTOFF * stdDev;
//	cout << "Edge cutoff: " << weightCutoff << "\n";
//	cout << "Minimum weight: " << minWeight << "\n";

	{
		vector<vecFloat> swapVec; swapVec.reserve
			(compGraph.size());

		for (int i = 0; i < compGraph.size(); i ++)
		{
			if (compGraph[i][2] < WEIGHT_CUTOFF)
			{
				continue;
			}

			swapVec.push_back(compGraph[i]);
		}

		swapVec.swap(compGraph);
	}
}

// Function: writeGraph
// Used to write one weighted graph into file

void writeGraph(vector<vector<float> >& compGraph, 
		string& compGraphPath)
{
	if (compGraph.size() <= 0)
	{
		cerr << "Error: The size of " 
		<< "compatible graph is wrong.\n";
		exit(1);
	}

	ofstream outputHd(compGraphPath.c_str());

	if (!outputHd)
	{
		cerr << "Error: Cannot open file: " << 
			compGraphPath << " to write.\n";
		exit(1);
	}

	for (int i = 0; i < compGraph.size(); i ++)
	{
		if (compGraph[i][2] < SUM_WEIGHT_CUTOFF)
		{
			continue;
		}

		outputHd << compGraph[i][0] << 
			"\t" << compGraph[i][1] << 
			"\t" << compGraph[i][2] << "\n";
//			"\t" << "1" << "\n";
	}

	outputHd.close();
}

// Function: writeFullGraph
// Used to write one weighted graph into file

void writeFullGraph(vector<vector<float> >& compGraph, 
		string& compGraphPath)
{
	if (compGraph.size() <= 0)
	{
		cerr << "Error: The size of " 
		<< "compatible graph is wrong.\n";
		exit(1);
	}

	ofstream outputHd(compGraphPath.c_str());

	if (!outputHd)
	{
		cerr << "Error: Cannot open file: " << 
			compGraphPath << " to write.\n";
		exit(1);
	}

	for (int i = 0; i < compGraph.size(); i ++)
	{
		outputHd << compGraph[i][0] << 
			"\t" << compGraph[i][1] << 
			"\t1\n";
//			"\t" << "1" << "\n";
	}

	outputHd.close();
}

// Function: runClusterONE
// Used to run clusterONE

void runClusterONE(string& compGraphPath, string& 
		psModPath, int sizeCutoff, float 
		denCutoff)
{
//	cout << "Graph: " << compGraphPath << "\n";
	string shellCmd = 
		"java -jar ../cluster_one-1.0.jar ";
	shellCmd += compGraphPath;
	shellCmd += " -d ";
	ostringstream sepHd;
	sepHd << denCutoff;
	shellCmd += sepHd.str();
	shellCmd += " -s ";
	sepHd.str("");
	sepHd << sizeCutoff;
	shellCmd += sepHd.str();
	shellCmd += " > ";
	shellCmd += psModPath;
	system(shellCmd.c_str());
}

// Function: runClusterONEBySeed
// Used to run clusterONE

void runClusterONEBySeed(string& compGraphPath, string& 
		psModPath, string& seedPath, int sizeCutoff, 
		float 
		denCutoff)
{
//	cout << "Graph: " << compGraphPath << "\n";
	string shellCmd = 
		"java -jar ../cluster_one-1.0.jar ";
	shellCmd += compGraphPath;
	shellCmd += " -d ";
	ostringstream sepHd;
	sepHd << denCutoff;
	shellCmd += sepHd.str();
	shellCmd += " -s ";
	sepHd.str("");
	sepHd << sizeCutoff;
	shellCmd += sepHd.str();
	shellCmd += " --seed-method 'file(";
	shellCmd += seedPath;
	shellCmd += ")' > ";
	shellCmd += psModPath;
	system(shellCmd.c_str());
}

// Function: getNameList
// Used to get the name list of files with the given
// suffix

void getNameList(vector<string>& nameList, string& tail)
{
	char pathName[] = "./";
        DIR* dir = opendir(pathName);
        struct dirent* ptr;

        while ((ptr = readdir(dir)) != NULL)
        {
                if (strcmp(ptr -> d_name,".") == 0 ||
                        strcmp(ptr -> d_name,"..") == 0)
                {
                        continue;
                }

                string fileName(ptr -> d_name);
                string suffixName =
                	fileName.substr(fileName.substr
			(0, fileName.find_last_of('.')).
			find_last_of('.') + 1);
                
		if (suffixName.compare(tail) == 0)
                {
                        nameList.push_back(fileName);
                }
        }

        closedir(dir);
}

// Function: readMod
// Ysed to read the modules into program

void readMod(vector<mod>& psModSet, string& 
		name)
{
	ifstream inputHd(name.c_str());

	if (!inputHd)
	{
		cerr << "Cannot open file " << 
			name << ".\n";
		exit(1);
	}

	string strLn;

	while (getline(inputHd, strLn))
	{
		istringstream sepHd(strLn);
		mod lnMod; string strUnit;

		while (sepHd >> strUnit)
		{
			lnMod.node.push_back
			(atoi(strUnit.c_str()));
		}

		psModSet.push_back(lnMod);
	}

	inputHd.close();

	cout << "Finished reading pseudo modules into program.\n" 
		<< "There're altogether " << psModSet.size() << 
		" groups of pseudo modules discovered.\n";
}

// Function: checkInc
// Used to check the inclusive relationship between 
// two modules

bool checkInc(mod& mod1, mod& mod2)
{
	if (mod2.node.size() > mod1.node.size())
	{
		cerr << "Error: The size of modules " 
			<< "is wrong.\n";
		exit(1);
	}

	bool res = true;

	for (int i = 0; i < mod2.node.size(); i ++)
	{
		if (find(mod1.node.begin(), mod1.
			node.end(), mod2.node[i]) 
				== mod1.node.end())
		{
			res = false;
			break;
		}
	}

	return res;
}

// Function: getJaccardSim
// Used to get Jaccard similarity

float getJaccardSim(vector<int>& vec1, vector<int>& vec2)
{
	int smallSize = MIN(vec1.size(), vec2.size());

	if (smallSize == 0)
	{
		cerr << "Error: There's empty modules!.\n";
		exit(1);
	}

	int commonSize = 0;

	for (int i = 0; i < vec1.size(); i ++)
	{
		if (find(vec2.begin(), vec2.end(), vec1[i]) 
				!= vec2.end())
		{
			commonSize ++;
		}
	}

	return float(commonSize) / float(smallSize);
}

// Function: unionSet
// Used to merge multiple sets into one

void unionSet(mod newMod, vector<int>& modIdSet, 
	vector<mod>& modSet, vector<int>& modFlag)
{
	map<int, int> nodeFlagMap; newMod.node.clear();

	for (int i = 0; i < modIdSet.size(); i ++)
	{
		for (int j = 0; j < modSet[modIdSet[i]].node.size(); 
				j ++)
		{
			nodeFlagMap[modSet[modIdSet[i]].node[j]] 
				= 0;
		}

		modFlag[modIdSet[i]] = 1;
	}

	for (map<int, int>::iterator it = nodeFlagMap.begin(); 
			it != nodeFlagMap.end(); it ++)
	{
		newMod.node.push_back(it -> first);
	}
}

// Function: constructMod
// Used to reconstruct the modules according to 
// their similarities.

void constructMod(vector<mod>& modSet, vector<vector<int> >& 
		mergeModSet)
{
	vector<mod> swapVec; swapVec.reserve(modSet.size());
	vector<int> modFlag(modSet.size(), 0);

	for (int i = 0; i < mergeModSet.size(); i ++)
	{
		mod newMod;

		unionSet(newMod, mergeModSet[i], modSet, 
			modFlag);
		swapVec.push_back(newMod);
	}

	for (int i = 0; i < modFlag.size(); i ++)
	{
		if (modFlag[i] == 1)
		{
			continue;
		}

		swapVec.push_back(modSet[i]);
	}

	swapVec.swap(modSet);

	cout << "Finished re-constructing modules.\n" 
		<< "There're altogether " << modSet.size() 
		<< " modules.\n";
}

// Function: getBiWeight
// Used to get the weight of an edge in hyper graph

float getBiWeight(vector<int>& nodeSet1, vector<int>& nodeSet2, 
		float** origMat, int& rowNum, int& colNum, 
		int& freqCutoff)
{
	map<pairInt, int> edgeMap;

	for (int i = 0; i < rowNum; i ++)
	{
		if (getFreq(origMat[i], 2, colNum - 1) < freqCutoff)
		{
			continue;
		}

		if (find (nodeSet1.begin(), nodeSet1.end(), 
			int(origMat[i][0])) != nodeSet1.end() && 
			find (nodeSet2.begin(), nodeSet2.end(), 
			int(origMat[i][1])) != nodeSet2.end())
		{
			edgeMap[make_pair(int(origMat[i][0]), 
					int(origMat[i][1]))] = i;
		}
		else if (find (nodeSet1.begin(), nodeSet1.end(), 
			int(origMat[i][1])) != nodeSet1.end() && 
			find (nodeSet2.begin(), nodeSet2.end(), 
			int(origMat[i][0])) != nodeSet2.end())
		{
			edgeMap[make_pair(int(origMat[i][1]), 
					int(origMat[i][0]))] = i;
		}
	}

	map<int, int> nodeMap1, nodeMap2;

	for (map<pairInt, int>::iterator it = edgeMap.begin(); 
		it != edgeMap.end(); it ++)
	{
		nodeMap1[(it -> first).first] = 0;
		nodeMap2[(it -> first).second] = 0;
	}

	if (nodeMap1.size() * nodeMap2.size() == 0)
	{
		return 0;
	}

//	cout << nodeMap1.size() << " " << nodeMap2.size() << "\n";
	float den = float(edgeMap.size()) / 
		float(nodeMap1.size() * nodeMap2.size());

	vecFloat weightVec(colNum - 2, 1);

	for (map<pairInt, int>::iterator it = edgeMap.begin(); 
		it != edgeMap.end(); it ++)
	{
		for (int i = 0; i < colNum - 2; i ++)
		{
			weightVec[i] *= origMat[it -> second][i + 2];
			
//			if (origMat[it -> second][i + 2] > 0)
//			{
//				weightVec[i] *= 1;
//			}
//			else
//			{
//				weightVec[i] *= 0;
//			}
		}
	}

	float weight = 0;

	for (int i = 0; i < colNum - 2; i ++)
	{
		weight += weightVec[i];
//		cout << weightVec[i] << "\n";
	}

	weight *= den;
//	cout << den << "\n";
//	cout << weight << "\n";
	return weight;
}

// Function: mergeMod
// Used to merge similar modules

void mergeMod(vector<mod>& mergedModSet, vector<vecMod>& 
	modSet, string& mergePath, 
	string& mergeModPath, float& jaccardSimCutoff, 
	float** origMat, int& rowNum, int& colNum, int& 
	freqCutoff, float& denCutoff, string& origPath)
{
	vector<vecFloat> simGraph;
	map<pairInt, int> modIdMap;
	map<int, pairInt> idModMap;

	int iter = 0; float maxWeight = 0;
	float minWeight = INT_MAX;

	for (int i = 0; i < modSet.size(); i ++)
	{
		for (int k = 0; k < modSet[i].size(); k ++)
		{
			pairInt tmpPair = make_pair(i, k);
			idModMap[iter] = tmpPair;
			modIdMap[tmpPair] = iter;

//			cout << tmpPair.first << " " << 
//				tmpPair.second << " " << 
//				iter << "\n";

			iter ++;
		}
	}

//	for (map<pairInt, int>::iterator it = modIdMap.begin(); 
//			it != modIdMap.end(); it ++)
//	{
//		cout << (it -> first).first << " " << 
//			(it -> first).second << "\t" << 
//			it -> second << "\n";
//	}

	for (int i = 1; i < modSet.size(); i ++)
	{
		for (int j = 0; j < i; j ++)
		{
			for (int k = 0; k < modSet[i].size(); k ++)
			{
				for (int l = 0; l < modSet[j].size(); l ++)
				{
					vecFloat sglEdge;
					float weight = getBiWeight
						(modSet[i][k].node, 
						 modSet[j][l].node, 
						 origMat, rowNum, 
						 colNum, freqCutoff);

					if (weight == 0)
					{
						continue;
					}

//					idModMap[iter] = make_pair(i, k);
//					sglEdge.push_back(iter);
//					iter ++;
//					idModMap[iter] = make_pair(j, l);
//					sglEdge.push_back(iter);
//					iter ++;
					pairInt pair1 = make_pair(i, k);
					pairInt pair2 = make_pair(j, l);
					int id1 = modIdMap[pair1];
					int id2 = modIdMap[pair2];

//					cout << pair1.first << " " << 
//						pair1.second << "\t" << 
//						pair2.first << " " << 
//						pair2.second << "\n";
//					cout << id1 << "\t" << id2 << 
//						"\n";

					sglEdge.push_back(id1);
					sglEdge.push_back(id2);
					sglEdge.push_back(weight);
					simGraph.push_back(sglEdge);

					if (maxWeight < weight)
					{
						maxWeight = weight;
					}

					if (minWeight > weight)
					{
						minWeight = weight;
					}
				}
			}
		}
	}

	cout << "Finished building the module similarity graph.\n" 
		<< "There're " << simGraph.size() << " edges.\n";

	mergedModSet.clear();

	if (simGraph.size() == 0)
	{
		for (int i = 0; i < modSet.size(); i ++)
		{
			for (int j = 0; j < modSet[i].size(); j ++)
			{
				mergedModSet.push_back(modSet[i][j]);
			}
		}
	}
	else
	{
		map<int, int> modFlagMap;

		for (map<pairInt, int>::iterator it = modIdMap.begin(); 
				it != modIdMap.end(); it ++)
		{
			modFlagMap[it -> second] = 0;
		}

//		float avgWeight = 0; float stdDev = 0;

		for (int i = 0; i < simGraph.size(); i ++)
		{
			simGraph[i][2] = BASE_WEIGHT + 
				(simGraph[i][2] - minWeight) * 
				(1 - BASE_WEIGHT) / 
				(maxWeight - minWeight);
//			modFlagMap[simGraph[i][0]] = 0;
//			modFlagMap[simGraph[i][1]] = 0;
//			simGraph[i][2] /= maxWeight;
//			avgWeight += simGraph[i][2];
//			stdDev += simGraph[i][2] * 
//				simGraph[i][2];
		}

//		avgWeight /= simGraph.size();
//		stdDev /= simGraph.size();
//		stdDev -= avgWeight * avgWeight;
//		stdDev = sqrt(stdDev);
//		float weightCutoff = avgWeight + stdDev 
//			* EDGE_CUTOFF2;
		
//		{
//			vector<vecFloat> swapVec;
//			swapVec.reserve(simGraph.size());
//
//			for (int i = 0; i < simGraph.size(); 
//					i ++)
//			{
//				if (simGraph[i][2] < weightCutoff)
//				{
//					continue;
//				}
//
//				swapVec.push_back(simGraph[i]);
//			}
//
//			swapVec.swap(simGraph);
//		}

		string simGraphPath = origPath;
		simGraphPath += ".simGraph.tmp";
		string simModPath = origPath;
		simModPath += ".simMod.tmp";
		ofstream simGraphHd(simGraphPath.c_str());

		if (!simGraphHd)
		{
			cerr << "Error: Cannot open file " 
				<< simGraphPath << ".\n";
			exit(0);
		}

		for (int i = 0; i < simGraph.size(); i ++)
		{
			simGraphHd << simGraph[i][0] << "\t" 
				<< simGraph[i][1] << "\t" << 
				simGraph[i][2] << "\n";
		}

		simGraphHd.close();

		float actDenCutoff = denCutoff;
		actDenCutoff *= BASE_WEIGHT;
		runClusterONE(simGraphPath, simModPath, 2, 
				actDenCutoff);

		ifstream simModHd(simModPath.c_str());

		if (!simModHd)
		{
			cerr << "Error: Cannot open file " 
				<< simModPath << ".\n";
			exit(0);
		}

		vector<vecInt> mergedModVec; string strLn;

		while (getline(simModHd, strLn))
		{
			if (strLn.empty())
			{
				continue;
			}

			string strSep; vecInt sglVec;
			istringstream sepHd(strLn);

			while (sepHd >> strSep)
			{
				sglVec.push_back
				(atoi(strSep.c_str()));
			}

			mergedModVec.push_back(sglVec);
		}

		simModHd.close();

		for (int i = 0; i < mergedModVec.size(); i ++)
		{
			for (int j = 0; j < mergedModVec[i].size(); 
					j ++)
			{
				modFlagMap[mergedModVec[i][j]] = 1;
			}
		}

		for (map<int, int>::iterator it = modFlagMap.begin(); 
				it != modFlagMap.end(); it ++)
		{
			if (it -> second == 0)
			{
				vecInt sglVec;
//				cout << it -> first << "\n";
				sglVec.push_back(it -> first);
				mergedModVec.push_back(sglVec);
			}
		}

//		for (int i = 0; i < mergedModVec.size(); i ++)
//		{
//			for (int j = 0; j < mergedModVec[i].size(); j ++)
//			{
//				cout << mergedModVec[i][j] << " ";
//			}
//
//			cout << "\n";
//		}

		for (int i = 0; i < mergedModVec.size(); i ++)
		{
			map<int, int> modNodeMap;

			for (int j = 0; j < mergedModVec[i].size(); j ++)
			{
				int k = idModMap[mergedModVec[i][j]].first;
				int l = idModMap[mergedModVec[i][j]].second;

				for (int m = 0; m < modSet[k][l].node.size(); 
						m ++)
				{
					modNodeMap[modSet[k][l].node[m]] 
						= 0;
				}
			}

			mod sglMod; sglMod.node.reserve(modNodeMap.size());

			for (map<int, int>::iterator it = modNodeMap.begin(); 
				it != modNodeMap.end(); it ++)
			{
				sglMod.node.push_back(it -> first);
			}

			mergedModSet.push_back(sglMod);
		}
	}

	cout << "Finished merging the pseudo modules.\n" 
		<< "There're altogether " << mergedModSet.size() 
		<< " pseudo modules.\n";
}


// Function: copyVec
// Used to copy a number of values from one vector
// to another vector

void copyVec(vector<int>& smallVec, vector<int>& bigVec, 
		int start, int end)
{
	smallVec.clear();

	if (end > bigVec.size() - 1)
	{
		end = bigVec.size() - 1;
	}

	smallVec.reserve(end - start + 1);

	for (int i = start; i <= end; i ++)
	{
		smallVec.push_back(bigVec[i]);
	}
}

// Function: localCL
// Used to clustering in a local manner

void localCL(vector<mod>& mergedModSet, vector<mod>& modSet, 
		float** origMat, int& rowNum, 
	int& colNum, int& freqCutoff, float& 
	denCutoff, string& origPath, int& sizeCutoff)
{
	mergedModSet.clear();

	for (int i = 0; i < modSet.size(); i ++)
	{
		map<pairInt, vecFloat> localMap;
		getLocalMap(localMap, modSet[i].node, origMat, 
				rowNum, colNum, freqCutoff);

		if (localMap.size() == 0)
		{
				continue;
		}

		vector<vector<float> > compGraph;
		compGraph.reserve(localMap.size());
		float actDenCutoff = denCutoff;
		actDenCutoff *= BASE_WEIGHT;
		getCompGraph(compGraph, localMap, actDenCutoff);
                string compGraphPath = origPath;
                ostringstream sepHd; sepHd << (i + 1);
                compGraphPath += ".";
                compGraphPath += sepHd.str();
                compGraphPath += ".compGraph.tmp";
                writeGraph(compGraph, compGraphPath);
                string psModPath = origPath;
                psModPath += ".";
                psModPath += sepHd.str();
                psModPath += ".psMod.tmp";
                runClusterONE(compGraphPath, psModPath,
                        sizeCutoff, actDenCutoff);
		readMod(mergedModSet, psModPath);
	}

	cout << "Finished clustering on the pseudo graph.\n";
}

// Function: getBicMat
// Used to get the matrix for biclustering

void getBicMat(vector<vecInt>& unitMat, vector<int>& 
	listOfNodes, float** origMat, int& rowNum, int& 
	colNum, map<pairInt, int>& edgeIdMap, 
	map<int, pairInt>& idEdgeMap, int& sizeCutoff)
{
//	if (listOfNodes.size() < sizeCutoff)
//	{
//		cerr << "Warning: Small pseudo module.\n";
//		return;
//	}

	map<int, int> nodeMap;

	for (vector<int>::iterator it = listOfNodes.begin(); 
		it != listOfNodes.end(); it ++)
	{
		nodeMap[*it] = 0;
	}

	int unitId = 0;

	for (int i = 0; i < rowNum; i ++)
	{
		if (nodeMap.count(int(origMat[i][0])) > 0 && 
			nodeMap.count(int(origMat[i][1])) > 0)
		{
			pairInt unitPair;
			unitPair = make_pair(origMat[i][0], 
					origMat[i][1]);
			unitId ++;
//			cout << unitId << "\n";
//			cout << unitId << " " << 
//				unitPair.first << "_" 
//				<< unitPair.second << "\n";
			edgeIdMap[unitPair] = unitId;
			idEdgeMap[unitId] = unitPair;
			vector<int> unitLn;
			unitLn.reserve(colNum - 2);

			for (int j = 2; j < colNum; j ++)
			{
				if (origMat[i][j] > 0)
				{
					unitLn.push_back
						(1);
				}
				else
				{
					unitLn.push_back
						(0);
				}
			}

			unitMat.push_back(unitLn);
		}
	}
}

// Function: addEdge
// Used to add edges to one module

void addEdge(mod& thisMod, vector<pairInt>& edgeSet)
{
	map<int, int> mapNode;

	for (int i = 0; i < thisMod.node.size(); i ++)
	{
		mapNode[thisMod.node[i]] = 0;
	}

	for (int i = 0; i < edgeSet.size(); i ++)
	{
		if (mapNode.count(edgeSet[i].first) == 0 || 
			mapNode.count(edgeSet[i].second) == 0)
		{
			continue;
		}

		thisMod.edge.push_back(edgeSet[i]);
	}
}

// Function: getCommonNode
// Used to get the common node between one list and one 
// vector

int getCommonNode(list<int>& thisList, vector<int>& thisVec)
{
	map<int, int> thisMap; int number = 0;

	for (list<int>::iterator it = thisList.begin(); it != 
			thisList.end(); it ++)
	{
		thisMap[*it] = 0;
	}

	for (int i = 0; i < thisVec.size(); i ++)
	{
		if (thisMap.count(thisVec[i]) > 0)
		{
			number ++;
		}
	}

	return number;
}

// Function: getFullMod
// Used to get the full modules

void getFullMod(graph& G, vector<mod>& psModSet)
{
	getNeighbor(G, 0);

	vector<int> nodeToAdd; map<int, int> nodeIncluded;

	for (int i = 0; i < psModSet.size(); i ++)
	{
		for (int j = 0; j < psModSet[i].node.size(); 
				j ++)
		{
			nodeIncluded[psModSet[i].node[j]] 
				= 0;
		}
	}

	cout << "There're altogether " << nodeIncluded.size() 
		<< " nodes covered by the preliminary modules.\n";

	for (map<int, int>::iterator it = G.nodeList.begin(); 
			it != G.nodeList.end(); it ++)
	{
		if (nodeIncluded.count(it -> first) == 0)
		{
			nodeToAdd.push_back(it -> first);
		}
	}

	cout << "There're altogether " << nodeToAdd.size() << 
		" nodes to be added to the modules.\n";

	vector<int> vecLn;
	vector<vecInt> patchVec(psModSet.size(), vecLn);

	for (int i = 0; i < nodeToAdd.size(); i ++)
	{
		int *scoreArray, *idArray;
		scoreArray = new int [psModSet.size()];
		idArray = new int [psModSet.size()];

		for (int j = 0; j < psModSet.size(); j ++)
		{
			idArray[j] = j;

			scoreArray[j] = getCommonNode
				(G.neighborList[nodeToAdd[i]], 
					psModSet[j].node);
		}

		quickSort(scoreArray, 0, psModSet.size() - 1, 
				idArray);

		for (int j = 0; j < psModSet.size(); j ++)
		{
			if (j >= TOP)
			{
				break;
			}

			patchVec[idArray[psModSet.size() - 1 - j]]
				.push_back(nodeToAdd[i]);
		}

		delete [] scoreArray; delete [] idArray;
	}

	for (int i = 0; i < psModSet.size(); i ++)
	{
		psModSet[i].node.insert(psModSet[i].node.end(), 
			patchVec[i].begin(), patchVec[i].end());
	}
}

// Function: bicCL
// Used to perform biclustering

void bicCL(vector<mod>& psModSet, 
		float** origMat, int& rowNum, 
		int& colNum, int& freqCutoff, int& sizeCutoff, 
		float& denCutoff, vector<mod>& realModSet)
{
	if (psModSet.size() == 0)
	{
		cerr << "There's no pseudo modules discovered.\n";
		exit(0);
	}

	for (int i = 0; i < psModSet.size(); i ++)
	{
		vector<vecInt> unitBicMat;
		map<pairInt, int> edgeIdMap;
		map<int, pairInt> idEdgeMap;

		if (psModSet[i].node.size() < sizeCutoff)
		{
			continue;
		}

		getBicMat(unitBicMat, psModSet[i].node, 
			origMat, rowNum, colNum, edgeIdMap, 
			idEdgeMap, sizeCutoff);

		if (unitBicMat.size() == 0)
		{
			continue;
		}

		string idStr; ostringstream idHd;
		idHd << (i + 1);
		string psMatPath = "PS." + idHd.str();
		psMatPath += ".psMat.tmp";
		string psBlkPath = "PS." + idHd.str();
		psBlkPath += ".psBlk.tmp";
		string realModPrefix = "REAL." + idHd.str();
		realModPrefix += '.';
		float height = sizeCutoff * (sizeCutoff - 1) / 2;
		height *= denCutoff;
		int width = freqCutoff;

		if (height < sizeCutoff)
		{
			height = sizeCutoff;
		}

		ofstream psMatHd(psMatPath.c_str());

		psMatHd << unitBicMat.size() << " " 
			<< unitBicMat[0].size() << " " 
			<< ceil(height) << " " << width << "\n";

		for (int j = 0; j < unitBicMat.size(); j ++)
		{
			for (int k = 0; k < unitBicMat[j].size() - 1; 
					k ++)
			{
				psMatHd << unitBicMat[j][k] << " ";
			}

			psMatHd << unitBicMat[j][unitBicMat[j].size() - 1] 
				<< "\n";
		}

		psMatHd.close();

		string bimaxCmd = "../bimax " + psMatPath;
		bimaxCmd += " > ";
		bimaxCmd += psBlkPath;

		system(bimaxCmd.c_str());

		ifstream psBlkHd(psBlkPath.c_str());

		if (!psBlkHd)
		{
			cerr << "Cannot open file " << 
				psBlkPath << ".\n";
			exit(1);
		}

		string blkLn; int secId = 0;

		while(getline(psBlkHd, blkLn))
		{
			if (blkLn.empty())
			{
				getline(psBlkHd, blkLn);
				getline(psBlkHd, blkLn);
				istringstream sepHd(blkLn);
				string strSep; secId ++;
				ostringstream secSep;
				secSep << secId;
				string realGraphPath = 
				realModPrefix + secSep.str();
				realGraphPath += ".realGraph.tmp";
				string realModPath = 
				realModPrefix + secSep.str();
				realModPath += ".realMod.tmp";
				vector<pairInt> edgeSet;

				while(sepHd >> strSep)
				{
					edgeSet.push_back
					(idEdgeMap[atoi
					 (strSep.c_str())]);
				}

				ofstream realGraphHd(realGraphPath.c_str());

				for (int j = 0; j < edgeSet.size(); j ++)
				{
					realGraphHd << edgeSet[j].first << 
						"\t" << edgeSet[j].second << 
						"\t1\n";
				}

				realGraphHd.close();

				runClusterONE(realGraphPath, realModPath, 
						sizeCutoff, denCutoff);

				ifstream realModHd(realModPath.c_str());
				string modLn;
				
				while (getline(realModHd, modLn))
				{
					istringstream modSepHd(modLn);
					mod sglMod; string sepStr;

					while(modSepHd >> sepStr)
					{
						sglMod.node.push_back
						(atoi(sepStr.c_str()));
					}

					addEdge(sglMod, edgeSet);

					realModSet.push_back(sglMod);
//					cout << sglMod.edge.size() << "\n";
				}

				realModHd.close();
			}
		}
	}

	cout << "Finished biclustering for the edge sets.\n" 
		<< "There're altogether " << realModSet.size() 
		<< " modules discovered.\n";
}

// Function: getModWeight
// Used to get the weight of module

void getModWeight(mod& thisMod, float** origMat, 
	int& rowNum, int& colNum)
{
	map<pairInt, int> mapEdge;

//	cout << "Edge: " << thisMod.edge.size() << "\n";

	for (int i = 0; i < thisMod.edge.size(); i ++)
	{
		mapEdge[thisMod.edge[i]] = 0;
	}

	thisMod.weight = 0; thisMod.occurr.clear();

	for (int i = 2; i < colNum; i ++)
	{
		thisMod.occurr.push_back(1);
	}

	for (int i = 0; i < rowNum; i ++)
	{
		pairInt thisPair = make_pair
			(origMat[i][0], 
				origMat[i][1]);

		if (mapEdge.count(thisPair) == 0)
		{
			continue;
		}

//		cout << thisPair.first << "\t" << 
//			thisPair.second << "\n";

		for (int j = 2; j < colNum; j ++)
		{
			if (origMat[i][j] > 0)
			{
				thisMod.occurr[j - 2] *= 1;
			}
			else
			{
				thisMod.occurr[j - 2] *= 0;
			}

		}

	}

	for (int i = 0; i < rowNum; i ++)
	{
		pairInt thisPair = make_pair
			(origMat[i][0], 
				origMat[i][1]);

		if (mapEdge.count(thisPair) == 0)
		{
			continue;
		}

		for (int j = 0; j < thisMod.occurr.size(); j ++)
		{
			if (thisMod.occurr[j] == 0)
			{
				continue;
			}

			thisMod.weight += origMat[i][2 + j];
		}
	}

	thisMod.freq = 0;

	for (int i = 0; i < colNum - 2; i ++)
	{
		thisMod.freq += thisMod.occurr[i];
	}

//	cout << thisMod.occurr.size() << "\n";
}

// Fnction: node2edge
// Used to transform nodes into edges Idiscarded)
// get the module weight

//void node2edge(mod& sglMod, float** origMat, int& rowNum, 
//		int& colNum, int& freqCutoff)
//{
//	map<int, int> sglMap;
//
//	for (int i = 0; i < sglMod.node.size(); i ++)
//	{
//		sglMap[sglMod.node[i]] = 0;
//	}
//
//	for (int i = 0; i < rowNum; i ++)
//	{
//		if (sglMap.count(int(origMat[i][0])) == 0 
//			|| sglMap.count(int(origMat[i][1])) 
//			== 0)
//		{
//			continue;
//		}
//
//		sglMod.occurr.clear();
//		sglMod.weight = 0;
//
//		for (int j = 2; j < colNum; j ++)
//		{
//			if (origMat[i][j] != 0)
//			{
//				sglMod.occurr.push_back(1);
//				sglMod.weight += 
//				origMat[i][j];
//			}
//			else
//			{
//				sglMod.occurr.push_back(0);
//				sglMod.weight += 
//				origMat[i][j];
//			}
//		}
//
//		sglMod.freq = 0;
//
//		for (int j = 0; j < sglMod.occurr.size(); j ++)
//		{
//			sglMod.freq += sglMod.occurr[j];
//		}
//
//		if (sglMod.freq < freqCutoff)
//		{
//			continue;
//		}
//
//		pairInt sglEdge = make_pair
//			(int(origMat[i][0]), 
//			 int(origMat[i][1]));
//		sglMod.edge.push_back(sglEdge);
//	}
//}

// Function: refineMod
// Used to refine modules

void refineMod(vector<mod>& modSet, float** origMat, 
		int& rowNum, int& colNum, int& freqCutoff, 
		float& denCutoff)
{
	if (modSet.size() == 1)
	{
		getModWeight(modSet[0], origMat, rowNum, 
				colNum);
		return;
	}

	float* weightArray; int* idArray;
	
	weightArray = new float [modSet.size()];
	idArray = new int [modSet.size()];

	for (int i = 0; i < modSet.size(); i ++)
	{
//		node2edge(modSet[i], origMat, rowNum, 
//				colNum, freqCutoff);
		getModWeight(modSet[i], origMat, rowNum, 
				colNum);
		weightArray[i] = modSet[i].weight;
		idArray[i] = i;
	}

	quickSort(weightArray, 0, modSet.size() - 1,
			idArray);

	{
		vector<mod> swapVec;
		swapVec.reserve(modSet.size());

		for (int i = modSet.size() - 1; i >= 0; 
				i --)
		{
			swapVec.push_back(modSet[idArray[i]]);
		}

		swapVec.swap(modSet);
	}

	{
		vector<mod> swapVec; vector<int> flagVec
			(modSet.size(), 1);
		swapVec.reserve(modSet.size());

		for (int i = 0; i < modSet.size() - 1; i ++)
		{
			if (flagVec[i] == 0)
			{
				continue;
			}

			for (int j = i + 1; j < modSet.size(); 
					j ++)
			{
				if (flagVec[j] == 0)
				{
					continue;
				}

				if (getJaccardSim(modSet[i].node, 
					modSet[j].node) >= 
						JACCARD_CUTOFF)
				{
					flagVec[j] = 0;
				}
			}
		}

		for (int i = 0; i < flagVec.size(); i ++)
		{
			if (flagVec[i] == 0)
			{
				continue;
			}

			swapVec.push_back(modSet[i]);
		}

		swapVec.swap(modSet);
	}

	delete [] idArray;
	delete [] weightArray;
	cout << "Finished refining the modules.\n" << 
		"There're altogether " << modSet.size() << 
		" modules.\n";
}

// Function: writeMod
// Used to write modules into file

void writeMod(vector<mod>& modSet, string& modPath)
{
	ofstream modHd(modPath.c_str());

	modHd << "Frequency;Weight;Node;Edge;Occurrence\n";

	for (int i = 0; i < modSet.size(); i ++)
	{
//		cout << "Id: " << i << "\n";

		modHd << modSet[i].freq << ";" << 
			modSet[i].weight << ";"; 
		
		for (int j = 0; j < modSet[i].node.
				size() - 1; j ++)
		{
			modHd << modSet[i].node[j] 
				<< ",";
		}

//		cout << "Line: 2058\n";

		modHd << modSet[i].node[modSet[i].node.
			size() - 1] << ";";

		for (int j = 0; j < modSet[i].edge.
				size() - 1; j ++)
		{
			modHd << modSet[i].edge[j].first << 
				"-" << modSet[i].edge[j].
				second << ",";
		}

//		cout << "Line: 2070\n";

		modHd << modSet[i].edge[modSet[i].edge.size() 
			- 1].first << "-" << modSet[i].edge
			[modSet[i].edge.size() - 1].second << 
			";";

//		cout << "Line: 2077\n";
//		cout << modSet[i].occurr.size() << "\n";

		for (int j = 0; j < modSet[i].occurr.size() 
				- 1; j ++)
		{
			modHd << modSet[i].occurr[j] << ",";
		}

//		cout << "Line: 2084\n";

		modHd << modSet[i].occurr[modSet[i].occurr.
			size() - 1] << "\n";
	}

	modHd.close();

	cout << "Finished writing the " << modSet.size() << 
		" modules into file path: " 
		<< modPath << ".\n";
}

// Function: miMod
// Used to integrate all functions

void miMod(string& origPath, int& freqCutoff, 
	float& denCutoff, int& sizeCutoff, float& 
	jaccardSimCutoff)
{
	float** origMat; int rowNum, colNum;
	getDim(origPath, rowNum, colNum);
//	cout << colNum << "\n";
	
	origMat = new float* [rowNum];
	
	for (int i = 0; i < rowNum; i ++)
	{
		origMat[i] = new float [colNum];
	}

	graph sumGraph;
	readData(origPath, origMat, rowNum, colNum, 
			sumGraph, freqCutoff);
	normGraph(sumGraph);

//	for (int i = 0; i < sumGraph.edgeList.size(); i ++)
//	{
//		cout << sumGraph.edgeList[i][0] << " " 
//			<< sumGraph.edgeList[i][1] << "\n";
//	}

//	optimalComp(sumGraph);
	
	string workDir = origPath;
	workDir += "_Work_tmp";
	mkdir(workDir.c_str(), 00700);
        int chdirStatus = chdir(workDir.c_str());

	if (chdirStatus == -1)
	{
		cerr << "Error: Cannot open directory: " 
			<< workDir << ".\n";
		exit(1);
	}

        char currentDir[DIR_LENGTH];
        getcwd(currentDir, DIR_LENGTH);
        cout << "Current directory: " << currentDir
                << ".\n";

	string sumGraphPath = origPath;
	sumGraphPath += ".sumGraph.tmp";
	string sumModPath = origPath;
	sumModPath += ".sumMod.tmp";
	writeGraph(sumGraph.edgeList, sumGraphPath);
	runClusterONE(sumGraphPath, sumModPath, sizeCutoff, 
			denCutoff * SUM_WEIGHT_CUTOFF);
//	string fullGraphPath = origPath;
//	fullGraphPath += ".fullGraph.tmp";
//	string fullModPath = origPath;
//	fullModPath += ".fullMod.tmp";
//	writeFullGraph(sumGraph.edgeList, fullGraphPath);
//	runClusterONEBySeed(fullGraphPath, fullModPath, 
//			sumModPath, sizeCutoff, 
//			denCutoff);
	vector<mod> psModSet;
	readMod(psModSet, sumModPath);
	getFullMod(sumGraph, psModSet);
	vector<mod> mergedModSet;
	localCL(mergedModSet, psModSet, origMat, rowNum, 
		colNum, freqCutoff, 
		denCutoff, origPath, sizeCutoff);
	vector<mod> realModSet;
	bicCL(mergedModSet, origMat, rowNum, 
			colNum, freqCutoff, sizeCutoff, 
			denCutoff, realModSet);
	
	if (realModSet.size() == 0)
	{
		cout << "There're 0 modules discovered.\n";
		exit(0);
	}

	refineMod(realModSet, origMat, rowNum, colNum, 
			freqCutoff, denCutoff);
	chdir("..");
        getcwd(currentDir, DIR_LENGTH);
        cout << "Current directory: " << currentDir
                << ".\n";

	string modPath = origPath;
	modPath += ".mod";

//	cout << "Line: 2156\n";

	writeMod(realModSet, modPath);

	for (int i = 0; i < rowNum; i ++)
	{
		delete [] origMat[i];
	}

	delete [] origMat;
}

// Function: main function
// Used to implement miMod

int main(int argc, char** argv)
{
	string origPath = argv[1];
	int freqCutoff = atoi(argv[2]);
	float denCutoff = atof(argv[3]);
	int sizeCutoff = atoi(argv[4]);
	float jaccardSimCutoff = atof(argv[5]);

	miMod(origPath, freqCutoff, denCutoff, 
		sizeCutoff, jaccardSimCutoff);

	return 0;
}
