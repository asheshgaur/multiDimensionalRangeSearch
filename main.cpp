#include <iostream>
using namespace std;
#include <vector>
#include <queue>
#include <stack>
#include <set>
#include <cstdlib>
#include <algorithm>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <functional>
#include <regex>
#include <ctime>
#include <random>
#include <chrono>
#include "classbench_reader.h"
#include "rules_and_packets.h"
#include "trace_tools.h"

// default constructor
rule::rule(int dimensions, int ruleCounter)
{
	this->ruleNumber = ruleCounter + 1;
    this->numDimensions = dimensions;
	fields.resize(5);
	prefix_length.resize(2);
}

// non-default constructor
rule::rule(vector<pair<unsigned int, unsigned int>>field, int num, int dimensions)
{
	this->ruleNumber = num;
	this->fields = field;
	this->numDimensions = dimensions;
}
rule::~rule() {
 //Don't need it cause no dynamic allocation.
}

//To check if a certain rule satisfies a packet or not.
bool rule::satisfiesPacket(const Packet& p) {
	for (int i = 0; i < MAX_DIMENSIONS; i++) {
		if ((fields[i].first <= p[i]) && (fields[i].second >= p[i]))
			continue;
		else return false;
	}
	return true;
}

bool rule::satisfiesPacketDim(const Packet& p, int dim) {
	if ((fields[dim].first <= p[dim]) && (fields[dim].second >= p[dim]))
		return true;
	else return false;
	
}

// display method for rule
void rule::printRule()
{
	cout << "Rule " << ruleNumber << ": " << endl;
	for(int i = 0; i < numDimensions; i++)
		cout << "Field" << i+1 << ": " << fields[i].first << " " << fields[i].second << endl;
}

int rule::getLR(int dimension){
	return fields[dimension - 1].first;
}

int rule::getRR(int dimension){
	return fields[dimension - 1].second;
}

int rule::getDimensionality(){
	return this->numDimensions;
}
unsigned int rule::getPrefixLength(int dim){
    return this->prefix_length[dim];
}

unsigned int rule::getDimRangeMin(int dimension){
    return fields[dimension].first;
}
unsigned int rule::getDimRangeMax(int dimension){
    return fields[dimension].second;
}
unsigned int rule::getRuleNumber(){
	return this->ruleNumber;
}
void rule::setDimRangeMin(int dimension, unsigned int value){
		this->fields[dimension].first = value;

}
void rule::setDimRangeMax(int dimension, unsigned int value){
		this->fields[dimension].second = value;

}

void rule::setPrefixLength(int dim, unsigned int value){
	this->prefix_length[dim] = value;
}

int maxDimensions = 5;
unsigned int maxRulesPerNode = 8;

//Global declarations
vector<rule*> rules; // vector for storing the rules from the classbench files.
vector<Packet> packets; // vector for storing randomly generated packets.
vector<pair<Packet*, set<int>>> classificationResult; // vector for storing the result of classification.

vector<unsigned int>* getEndPoints (vector<int>* ruleIDs, int d){
	vector<unsigned int>* eP = new vector<unsigned int> ();
	int numRules = (*ruleIDs).size();
	for (int i = 0; i < numRules; i++) {
        (*eP).push_back(rules[(*ruleIDs)[i]]->getDimRangeMin(d));
        (*eP).push_back(rules[(*ruleIDs)[i]]->getDimRangeMax(d));
	}

	return eP;
}

/* we may not need this program */
int countEndPoints (vector<unsigned int>* eP) {
	int count = 0;
	int numEP = (*eP).size();
	if (numEP == 0) return count;
	count = 1;
	for (int i = 1; i < numEP; i++) {
		if ((*eP)[i-1] != (*eP)[i]) count++;
	}

	return count;
}

class RangeTree {
public:
	unsigned int low;
	unsigned int high;
	int dimension;
	RangeTree* left;
	RangeTree* right;
	RangeTree* nextDimension;
	vector<int>* ruleIDs;
	RangeTree ();
	~RangeTree();
	RangeTree (unsigned int l, unsigned int h, int dim, RangeTree* lChild, RangeTree* rChild);
	bool insertARule (int ruleNum, int dim);
	bool satisfiesPacket(unsigned int point);
	bool satisfiesRule( int ruleNum, int dim);
	bool isleaf();
	int getDimension();
	void preOrder();
};
RangeTree::RangeTree () {
	low = 0;
	high = 0;
	dimension = 0;
	left = NULL;
	right = NULL;
	ruleIDs = NULL;
	nextDimension = NULL;
}

RangeTree::RangeTree (unsigned int l, unsigned int h, int dim, RangeTree* lChild, RangeTree* rChild) {
	low = l;
	high = h;
	dimension = dim;
	left = lChild;
	right = rChild;
	ruleIDs = NULL;
	nextDimension = NULL;
}

RangeTree::~RangeTree() {
	delete ruleIDs;
	delete left;
	delete right;
	delete nextDimension;
}
bool RangeTree::isleaf() {
	if ((this->left == NULL) && (this->right == NULL)) {
		return true;
	}
	else return false;
}

int RangeTree::getDimension(){
	return dimension;
}

bool RangeTree:: insertARule (int ruleNum, int dim) {

	bool inserted = false;
    if (satisfiesRule(ruleNum, dim)) {
		if (ruleIDs == NULL) {
			ruleIDs = new vector<int>();
		}
		(*ruleIDs).push_back(ruleNum);
		inserted = true;
	}
	return inserted;
}

bool RangeTree::satisfiesPacket(unsigned int point) {
	if ((low <= point) && (high >= point))
		return true;
	else return false;
}

bool RangeTree::satisfiesRule(int ruleNum, int dim) {
	if ((rules[ruleNum]->getDimRangeMin(dim) <= low) && (rules[ruleNum]->getDimRangeMax(dim) >= high))
		return true;
	else return false;
}


void RangeTree::preOrder () {
	stack<RangeTree*> s;
	RangeTree* t;
	s.push(this);
	while (!s.empty()) {
		t = s.top();
		s.pop();
		if ((*t).ruleIDs != NULL) {
			//cout << "(" << (*t).low << ": " << (*t).high << "); Size = " << (*(*t).ruleIDs).size() << "==";
			int numRules = (*(*t).ruleIDs).size();
			for (int k=0; k < numRules; k++)
				//cout << " " << (*(*t).ruleIDs)[k];
			//cout << endl;
			if ((*t).nextDimension != NULL) {
				//cout << "******************************** (1)" << endl;
				(*(*t).nextDimension).preOrder();
				//cout << "******************************** (2)" << endl;
			}
		}
		else
			//cout << "(" << (*t).low << ": " << (*t).high << "); Size = 0" << endl;
		if ((*t).right != NULL) s.push((*t).right);
		if ((*t).left != NULL) s.push((*t).left);
	}
}

void insertRules (RangeTree* node, vector<int>* ruleIDs, int dim) {
    // cout << endl << "Inside insert rules with node ("<<node->low
    //     << " : " << node->high << ") and dim: " << dim << endl;
	stack<RangeTree*> s;
	RangeTree* t;
	bool inserted;
	int numRules = (*ruleIDs).size();
	for (int i=0; i < numRules; i++) {
        // cout << "Inside for loop for iteration: " << i << endl;
		s.push(node);
		while (!s.empty()) {
			t = s.top();
			s.pop();
			//cout << "(" << (*t).low << ": " << (*t).high << ")" << endl;
			inserted = 	(*t).insertARule ((*ruleIDs)[i], dim);
			if (!inserted) {
				if ((*t).right != NULL) s.push((*t).right);
				if ((*t).left != NULL) s.push((*t).left);
			}
		}
	}
}
//This function only goes one way when a rule is not inserted into a node.
void insertRulesOnce (RangeTree* node, vector<int>* ruleIDs, int dim) {
    // cout << endl << "Inside insert rules with node ("<<node->low
    //     << " : " << node->high << ") and dim: " << dim << endl;
	stack<RangeTree*> s;
	RangeTree* t;
	bool inserted;

	int numRules = (*ruleIDs).size();
	for (int i=0; i < numRules; i++) {
        // cout << "Inside for loop for iteration: " << i << endl;
		s.push(node);
		while (!s.empty()) {
			t = s.top();
			s.pop();
			//cout << "(" << (*t).low << ": " << (*t).high << ")" << endl;
			inserted = 	(*t).insertARule ((*ruleIDs)[i], dim);
			if (!inserted) {
				if ((*t).right != NULL){
					if ((*t).right->satisfiesRule((*ruleIDs)[i], dim))
						s.push((*t).right);
				} 
				else if ((*t).left != NULL){
					if ((*t).left->satisfiesRule((*ruleIDs)[i], dim))
						s.push((*t).left);
				} 
			}
		}
	}
}

unsigned int inline atoui(const string& in) {
	std::istringstream reader(in);
	unsigned int val;
	reader >> val;
	return val;
}

RangeTree* buildTree (vector<int>* ruleIDs, int dim, const vector<int>& fieldOrder) {
	vector<unsigned int>* endPoints;
	RangeTree* aNode;
	RangeTree* root = NULL;
	stack<RangeTree*> s;
	RangeTree* t;
	queue<RangeTree*> myQ;
	unsigned int k;
    // cout << endl << "Inside buildTree for dimension: " << dim
    //     << endl << "And the number of rules: " << ruleIDs->size() << endl;
	unsigned int numRules = (*ruleIDs).size();
	if (maxRulesPerNode < numRules) {
	
		endPoints = getEndPoints (ruleIDs,fieldOrder[dim]);
		sort((*endPoints).begin(), (*endPoints).end());
		//int numEndPoints = countEndPoints (endPoints);
		//cout << "Number of end points = " << numEndPoints << endl;

		k = (*endPoints)[0];
		int numEP = (*endPoints).size();
		for (int i=1; i < numEP; i++) {
			if ((*endPoints)[i] != k) {
				aNode = new RangeTree(k, (*endPoints)[i], fieldOrder[dim], NULL, NULL);
				k = (*endPoints)[i];
				myQ.push(aNode);
				//cout << dim << ": >>>>> " << (*aNode).low << ": " << (*aNode).high << endl;
			}
		}
		delete endPoints;
		
		while (myQ.size() > 1) {
			RangeTree* l = myQ.front();
			myQ.pop();
			RangeTree* r = myQ.front();
			if ((*l).high == (*r).low) {
				myQ.pop();
				RangeTree* n = new RangeTree ((*l).low, (*r).high, fieldOrder[dim], l, r);
				//cout << dim << ": ++++++++" << (*n).low << ": " << (*n).high << endl;
				myQ.push(n);
			}
			else myQ.push(l);
		}

		root = myQ.front(); //root of the dimension one range tree

		//Insert the rules in the first dimension range tree - root

		//insertRulesOnce(root, ruleIDs, fieldOrder[dim]);
		insertRules(root, ruleIDs, fieldOrder[dim]);
		//Process the next Dimension, if dim < maxDimensions
		if ((dim+1) < MAX_DIMENSIONS) { //dim 0 is already procesed
			s.push(root);
			while (!s.empty()) {
				t = s.top();
				s.pop();
				if ((*t).right != NULL) s.push((*t).right);
				if ((*t).left != NULL) s.push((*t).left);
				if (((*t).ruleIDs != NULL) /*&& ((*(*t).ruleIDs).size() > 1)*/) {
					//cout << "@@@@: (" << (*t).low << ": " << (*t).high << ")" << endl;
					RangeTree* another = buildTree ((*t).ruleIDs, dim+1, fieldOrder);
					(*t).nextDimension = another;
					// if (((*t).low == 15) && ((*t).high == 20)) {
					// 	//cout << "$$$$$$$$$$$$$$$$$$$$$$$$" << endl;
					// 	(*another).preOrder();
					// }
				}
			}
		}

		//(*root).preOrder();
	}
	else {
		endPoints = getEndPoints (ruleIDs,fieldOrder[dim]);
		//sort((*endPoints).begin(), (*endPoints).end());
		int numEP = (*endPoints).size();
		uint32_t minEP = 4294967295;
		uint32_t maxEP = 0;
		//Finding min endpoint
		for (int i = 0; i < numEP; i++) {
			if (minEP > (*endPoints)[i]){
				minEP = (*endPoints)[i];
			}
			if (maxEP < (*endPoints)[i]) {
				maxEP = (*endPoints)[i];
			}
		}

		root = new RangeTree (minEP, maxEP, fieldOrder[dim], NULL,NULL);
		root->ruleIDs = new vector<int>();
		for (unsigned int i = 0; i < ruleIDs->size(); i++) {
			root->ruleIDs->push_back((*(ruleIDs))[i]);
		}

		// for (int i = 0; i < ruleIDs->size(); i++){
		// 	bool ruleInserted = (*root).insertARule ((*ruleIDs)[i], fieldOrder[dim]);
		// 	if (!ruleInserted) {
		// 		cout << std::endl << "Rule was not inserted in the tree with dim " << fieldOrder[dim] ;
		// 		cout << std::endl << "Tree low and high: " << k << " " << m;
		// 		cout << std::endl << rules[(*ruleIDs)[i]]->getDimRangeMin(fieldOrder[dim]);
		// 		cout << std::endl << rules[(*ruleIDs)[i]]->getDimRangeMax(fieldOrder[dim]);
		// 	}
		// }
	}

	return root;
}


void classify (RangeTree* root, Packet& packet, int dim, vector<int>* matches) {

	RangeTree* node = root;
	//cout << "%%%% packet[dim] = " << packet[dim] <<"; dim = " << dim << endl;
	//cout << "%%%% node values (" << (*node).low <<":" << (*node).high << ")" << endl;
	if (node != NULL) {
		if ((packet[dim] >= (*node).low) && (packet[dim] <= (*node).high)) {
			//cout << "Packet Attribute: " << dim << " = " << packet[dim];
			//cout << "; node values (" << (*node).low <<":" << (*node).high << ")" << endl;
			if (((*node).nextDimension == NULL) && ((*node).ruleIDs != NULL)) {
				int numCurrNodeRules = (*(*node).ruleIDs).size();
				for (int k=0; k < numCurrNodeRules; k++)
					(*matches).push_back((*(*node).ruleIDs)[k]);
					//cout << " " << (*(*node).ruleIDs)[k];
			}
			if ((*node).nextDimension != NULL) {
				//cout << "next dimension" << endl;
				//RangeTree* t = node;
				/*if ((*t).ruleIDs != NULL) {
					cout << "(" << (*t).low << ": " << (*t).high << "); Size = " << (*(*t).ruleIDs).size() << "==";
					for (int k=0; k < (*(*t).ruleIDs).size(); k++)
						cout << " " << (*(*t).ruleIDs)[k];
					cout << endl;
				}*/
				classify((*node).nextDimension, packet, dim+1, matches);
			}
		}
		if (!(((*node).low > packet[dim]) || ((*node).high < packet[dim]))) {
			if ((*node).left != NULL)
				classify((*node).left, packet, dim, matches);
			if ((*node).right != NULL) 
				classify((*node).right, packet, dim, matches);
		}
	}
} 
//The one used for numbers in the paper submitted.
void classifyFirst (RangeTree* root, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	RangeTree* node = root;
	bool found = false;

	if (node != NULL) {
		if ((packet[fieldOrder[dim]] >= (*node).low) && (packet[fieldOrder[dim]] <= (*node).high)) {
			if (((*node).nextDimension == NULL) && ((*node).ruleIDs != NULL)) {
				if ((*(*node).ruleIDs).size() > 0) {
					found = true;
					matchRule = (*(*node).ruleIDs)[0];
					//cout << "^^^^^ " << matchRule << endl;
					return;
				}
			}
			if (((*node).nextDimension != NULL) && (!found)) {
				//RangeTree* t = node;
				classifyFirst((*node).nextDimension, packet, dim+1, fieldOrder, matchRule);
			}
		}
		if (!found) {

			if (!(((*node).low > packet[dim]) || ((*node).high < packet[dim]))) {
				if ((*node).left != NULL)
					classifyFirst((*node).left, packet, dim, fieldOrder, matchRule);
				else if ((*node).right != NULL) 
					classifyFirst((*node).right, packet, dim, fieldOrder, matchRule);
			}
		}
	}
}
//This one classifies to a single rule, but still accurate since the rule classified to is the dont care rule.
void classifyFirstTest (RangeTree* root, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	RangeTree* node = root;
	//bool found = false;

	if (node != NULL) {
		if ((packet[fieldOrder[dim]] >= (*node).low) && (packet[fieldOrder[dim]] <= (*node).high)) {
			if ((node->ruleIDs != NULL) && (*(*node).ruleIDs).size() > 0) {
				if ((*node).nextDimension == NULL){
					matchRule = (*(*node).ruleIDs)[0];
					return;
				}
				else
					classifyFirstTest((*node).nextDimension, packet, dim+1, fieldOrder, matchRule);
			}
			else {
				if ((*node).left != NULL){
					if ((packet[fieldOrder[dim]] >= (*(*node).left).low) && (packet[fieldOrder[dim]] <= (*(*node).left).high)) {
						classifyFirstTest((*node).left, packet, dim, fieldOrder, matchRule);
					}
				}
				else if ((*node).right != NULL){
					if ((packet[fieldOrder[dim]] >= (*(*node).right).low) && (packet[fieldOrder[dim]] <= (*(*node).right).high)) {
						classifyFirstTest((*node).right, packet, dim, fieldOrder, matchRule);
					}
				}
				else return;
			}
		}
		else return;
	}
}

//Function for finding the most specific rule applicable to the packet. This one does classify to other rules but not accurate.
void classifyLastRecursive (RangeTree* node, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	//RangeTree* t;
	//bool found = false;

	if (!node) return;

	if ((*node).satisfiesPacket(packet[fieldOrder[dim]])) {
		if (((*node).isleaf()) || 
			(((*node).ruleIDs != NULL) && ((*(*node).ruleIDs).size() == 1))) { // if node is leaf or it only has a single rule.
			if ((*node).ruleIDs != NULL) matchRule = (*((*node).ruleIDs))[0]; 
			return;
		}
		else //Rule ambiguity is there 
			classifyLastRecursive ((*node).nextDimension, packet, dim + 1, fieldOrder, matchRule);
	}
	// Based on packet characteristics, push left or right node to the stack.
	if (((*node).left != NULL) && ((*(*node).left).satisfiesPacket(packet[node->getDimension()])))
		classifyLastRecursive ((*node).left, packet, dim, fieldOrder, matchRule);
	else if (((*node).right != NULL) && ((*(*node).right).satisfiesPacket(packet[node->getDimension()])))
		classifyLastRecursive ((*node).right, packet, dim, fieldOrder, matchRule);
	
	
	//No match was found
	return;

}

//Function for finding the most specific rule applicable to the packet. This one classifies to a single rule.
void classifyLastIterative (RangeTree* node, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	//RangeTree* t;
	stack<RangeTree*> s;
	//bool found = false;

	s.push(node);
	while (!s.empty()){
		//t = s.top();
		s.pop();
		if ((*node).satisfiesPacket(packet[node->getDimension()])) {
			if (((*node).isleaf()) || 
				(((*node).ruleIDs != NULL) && ((*(*node).ruleIDs).size() > 0))) { // if node is leaf or it only has a single rule.
				matchRule = (*((*node).ruleIDs))[0]; 
				return;
			}
			else //Rule ambiguity is there 
				s.push((*node).nextDimension);
		}
		// Based on packet characteristics, push left or right node to the stack.
		if (((*node).left != NULL) && ((*(*node).left).satisfiesPacket(packet[node->getDimension()])))
			s.push((*node).left);
		else if (((*node).right != NULL) && ((*(*node).right).satisfiesPacket(packet[node->getDimension()])))
			s.push((*node).left);
		
	}
	//No match was found
	return;

}

bool classifyLastIterativeBacktrack (RangeTree* node, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	std::stack<RangeTree*> s;

	if ((node != NULL) && (node->satisfiesPacket(packet[node->getDimension()])) /*&& ((node->ruleIDs) != NULL)*/) { //cannot force the root to have rules.
		if (node->ruleIDs != NULL){
			s.push(node);
		}
			
	}
	else return false;
	bool flag = false;

	while (!flag) {
		if ((node->right != NULL) && (node->right->satisfiesPacket(packet[fieldOrder[dim]]))) {
			if ((node->right->ruleIDs != NULL))
				s.push(node->right);
			node = node->right;
		}
		else if ((node->left != NULL) && (node->left->satisfiesPacket(packet[fieldOrder[dim]]))) {
			if ((node->left->ruleIDs != NULL))
				s.push(node->left);
			node = node->left;
		}
		else flag = true; //Leaf node was encountered so all relevant nodes are in the stack.
	} // While loop ends here. The stack 's' should have all nodes satisfying the packet in a single dimension.

	//The next while loop will go through these nodes in the stack and check the next dimension.

	while (!s.empty()) { //Backtracking through the nodes by popping them out now.
		RangeTree* y = s.top();
		s.pop(); 	//backtracking requires removing the node visited.

		if ((y->ruleIDs != NULL) && (y->ruleIDs->size() == 1)) { //just one rule
			//just check for the remaining dimensions (if there are any).
			bool ruleFound = false;
			if (y->nextDimension != NULL) {
				RangeTree* z = y->nextDimension;
				
				while (z != NULL) {
					if (z->satisfiesPacket(packet[z->getDimension()])) {
						if (z->nextDimension != NULL){
							z = z->nextDimension;
						}
							
						else {
							matchRule = (*(z->ruleIDs))[0];
							ruleFound = true;
							return true;
						}
					}
					else {
						ruleFound = false;
						break;
					} 
				}
				if (ruleFound) break;
			}
			//if (ruleFound) return;
		}
		else { //have to do the same stuff with the remaining dimensions.
			bool ruleFound = classifyLastIterativeBacktrack (y->nextDimension, packet, dim+1, fieldOrder, matchRule);
			
			if (ruleFound) return true;
		}
	}
	return false;
}

bool classifyLastIterativeBacktrackLeafSize (RangeTree* node, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	std::stack<RangeTree*> s;

	if ((node != NULL) && (node->satisfiesPacket(packet[node->getDimension()])) /*&& ((node->ruleIDs) != NULL)*/) { //cannot force the root to have rules.
		if (node->ruleIDs != NULL){
			s.push(node);
		}
			
	}
	else return false;
	bool flag = false;

	while (!flag) {
		if ((node->right != NULL) && (node->right->satisfiesPacket(packet[fieldOrder[dim]]))) {
			if ((node->right->ruleIDs != NULL))
				s.push(node->right);
			node = node->right;
		}
		else if ((node->left != NULL) && (node->left->satisfiesPacket(packet[fieldOrder[dim]]))) {
			if ((node->left->ruleIDs != NULL))
				s.push(node->left);
			node = node->left;
		}
		else flag = true; //Leaf node was encountered so all relevant nodes are in the stack.
	} // While loop ends here. The stack 's' should have all nodes satisfying the packet in a single dimension.

	//The next while loop will go through these nodes in the stack and check the next dimension.

	while (!s.empty()) { //Backtracking through the nodes by popping them out now.
		RangeTree* y = s.top();
		s.pop(); 	//backtracking requires removing the node visited.

		if (y->nextDimension == NULL) { 
			//looping through each rule in the node to see which one satisfies MAX_DIMENSIONS - dim fields.
			for (unsigned int i = 0; i < y->ruleIDs->size(); i++) {
				for (int j = dim; j < MAX_DIMENSIONS; j++) {
					if (rules[(*(y->ruleIDs))[i]]->satisfiesPacketDim(packet, fieldOrder[j])) {
						matchRule = (*(y->ruleIDs))[i];
						return true;
					}
				}
			}
			return false;
		}
		else { //have to do the same stuff with the remaining dimensions.
			if (classifyLastIterativeBacktrackLeafSize (y->nextDimension, packet, dim+1, fieldOrder, matchRule)) return true;
		}
	}
	return false;
}

bool classifyLastIterativeBacktrackLeafSizeAllRules (RangeTree* node, Packet packet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	std::stack<RangeTree*> s;

	if ((node != NULL) && (node->satisfiesPacket(packet[node->getDimension()])) /*&& ((node->ruleIDs) != NULL)*/) { //cannot force the root to have rules.
		if (node->ruleIDs != NULL){
			s.push(node);
		}
			
	}
	else return false;
	bool flag = false;

	while (!flag) {
		if ((node->right != NULL) && (node->right->satisfiesPacket(packet[fieldOrder[dim]]))) {
			if ((node->right->ruleIDs != NULL))
				s.push(node->right);
			node = node->right;
		}
		else if ((node->left != NULL) && (node->left->satisfiesPacket(packet[fieldOrder[dim]]))) {
			if ((node->left->ruleIDs != NULL))
				s.push(node->left);
			node = node->left;
		}
		else flag = true; //Leaf node was encountered so all relevant nodes are in the stack.
	} // While loop ends here. The stack 's' should have all nodes satisfying the packet in a single dimension.

	//The next while loop will go through these nodes in the stack and check the next dimension.

	while (!s.empty()) { //Backtracking through the nodes by popping them out now.
		RangeTree* y = s.top();
		s.pop(); 	//backtracking requires removing the node visited.

		if (y->nextDimension == NULL) { 
			//looping through each rule in the node to see which one satisfies MAX_DIMENSIONS - dim fields.
			for (unsigned int i = 0; i < y->ruleIDs->size(); i++) {
				for (int j = dim; j < MAX_DIMENSIONS; j++) {
					if (rules[(*(y->ruleIDs))[i]]->satisfiesPacketDim(packet, fieldOrder[j])) {
						matchRule = (*(y->ruleIDs))[i];
						return true;
					}
				}
			}
			return false;
		}
		else { //have to do the same stuff with the remaining dimensions.
			if (classifyLastIterativeBacktrackLeafSizeAllRules (y->nextDimension, packet, dim+1, fieldOrder, matchRule)) return true;
		}
	}
	return false;
}

//Called from within the classifyBruteForce() function.
void determineMostRelevantRule (const vector<int>& ruleSet, int dim, const vector<int>& fieldOrder, int& matchRule) {
	unsigned int minRangeSpan = 1;
	vector<int> minRules;
	if (dim == MAX_DIMENSIONS)
		return;
	if (ruleSet.size() == 1){
		matchRule = ruleSet[0];
		return;
	}
	int numRules = ruleSet.size();
	for (int i = 0; i < numRules; i++) {
		unsigned int rangeSpan = rules[ruleSet[i]]->getDimRangeMax(fieldOrder[dim]) - 
							rules[ruleSet[i]]->getDimRangeMin(fieldOrder[dim]);
		if (i == 0) {
			minRangeSpan = rangeSpan;
		}
		else {
			if (minRangeSpan > rangeSpan)
				minRangeSpan = rangeSpan;
		}

	}
	// After getting the value of minRangeSpan we will check the rules in ruleSet and push them into a vector to pass to next
	// dimension processing.

	for (int i = 0; i < numRules; i++) {
		unsigned int rangeSpan = rules[ruleSet[i]]->getDimRangeMax(fieldOrder[dim]) - 
							rules[ruleSet[i]]->getDimRangeMin(fieldOrder[dim]);
		if (rangeSpan == minRangeSpan) {
			minRules.push_back(ruleSet[i]);
		}
	}
	if (minRules.size() == 1) {
		matchRule = minRules[0];
		return;
	}
	else determineMostRelevantRule (minRules, dim + 1, fieldOrder, matchRule);
}
vector<int> classifyBruteForce () {
	vector<vector<int>> queryResult;		//Stores all the applicable rules for each packet.
	vector<int> mostRelevantRules;		//Stores the most relevant rule for each packet.
	//Checking which rules are satisfied by each packet.

	//A temporary backdoor to print the values of each relevant rule for each packet to a file.
	std::ofstream allRelevantRules ("AllApplicableRules.txt");
	int numPackets = packets.size();
	for (int i = 0; i < numPackets; i++){
		vector<int> ruleSet;
		queryResult.push_back(ruleSet);
		
		int numRules = rules.size();
		for (int j = 0; j < numRules; j++) {
			if (rules[j]->satisfiesPacket(packets[i])) {              
				queryResult[i].push_back(j);
			}
		}
		int numQueryResult = queryResult[i].size();
		allRelevantRules << numQueryResult;
		for ( int k = 0; k < numQueryResult; k++) {
			allRelevantRules << " "<< queryResult[i][k];
		}
		allRelevantRules << std::endl;
	}
	allRelevantRules.close();
	//Now to check the rule which is the most relevant one
	
	int matchRule = -1;
	int dim = 0;
	int numQueryResult = queryResult.size();
	for (int i = 0; i < numQueryResult; i++) {
		determineMostRelevantRule(queryResult[i], dim, fieldOrder, matchRule);
		mostRelevantRules.push_back(matchRule);
		matchRule = -1;
	}
	return mostRelevantRules;
}

vector<Packet> generatePacketsFromFile(string fileName) {
	ifstream inputFile (fileName);
	if (!inputFile) {
        std::cerr << "Failed to open the file for reading.\n";
        exit (0);
    }
	vector<Packet> Packets;
	std::string line;
    while (getline(inputFile, line)) {
        Packet packet;
		stringstream ss(line);
		string field;
        while (getline(ss, field, ' ')) {
			unsigned int intField = atoui(field);
			packet.push_back(intField);
		}
		Packets.push_back(packet);
    }
	return Packets;
}

vector<rule*> generateRulesFromFile(string fileName) {
	ifstream inputFile (fileName);
	if (!inputFile) {
        std::cerr << "Failed to open the file for reading.\n";
        exit (0);
    }
	vector<rule*> Rules;
	std::string line;
	int ruleCounter = 0;
    while (getline(inputFile, line)) {
		stringstream ss(line);
		string field;
		int pairCounter = 1;
		
		unsigned int intField1 = 0;
		unsigned int intField2 = 0;
		vector <pair<unsigned int, unsigned int>> fields;
        while (getline(ss, field, ' ')) {
			
			if (pairCounter % 2 != 0)
				intField1 = atoui(field);
			else {
				intField2 = atoui(field);
				pair<unsigned int, unsigned int> fieldPair = {intField1, intField2};
				fields.push_back(fieldPair);
			}	
			pairCounter++;
		}
		rule* Rule = new rule (fields, ruleCounter, MAX_DIMENSIONS);
		Rules.push_back(Rule);
		ruleCounter++;
    }
	return Rules;
}
// Function to calculate the mean
double calculateMean(const std::vector<double>& values) {
    return std::accumulate(values.begin(), values.end(), 0.0) / values.size();
}

// Function to calculate the standard deviation
double calculateStdDev(const std::vector<double>& values, double mean) {
    double sum = 0.0;
    for (auto value : values) {
        sum += (value - mean) * (value - mean);
    }
    return std::sqrt(sum / values.size());
}
int main (int argc, char* argv[]) {

	RangeTree* mainRoot;
	set<int> ruleMatches;
	
	int brute_force_choice = 0; //Default to no brute classification.
    int dim = 0;
	if (argc < 4)
	{
		cout << "Please use the following format:" << endl;
		cout << "./simulation {Rules_File} {maximum_nodes_per_leaf} {Brute = 0,1}\n";
		exit(-1);
	}
    //numRules = stoi(argv[1]); //numRules is a global variable.
    srand(time(NULL));
	rules = InputReader::ReadFilterFileClassBench(argv[1]);
	//rules = generateRulesFromFile(argv[1]);
    int numRules = rules.size();
	for (int i = 0; i < numRules; i++) {
		for (int j = 0; j < 5; j++) {
			if (rules[i]->getDimRangeMin(j) == rules[i]->getDimRangeMax(j))
				rules[i]->setDimRangeMax(j,rules[i]->getDimRangeMax(j) + 1);
		}
		
	}
	maxRulesPerNode = std::stoi(argv[2]);
	brute_force_choice = std::stoi(argv[3]);
	packets = GeneratePacketsFromRuleset(rules);
	
	vector<int>* ruleIDs = new vector<int>();
	for (int i =0; i < numRules; i++) {
		(*ruleIDs).push_back(i);
	}

    chrono::time_point<chrono::steady_clock> start = chrono::steady_clock::now();
	
	mainRoot = buildTree (ruleIDs, 0, fieldOrder);
    chrono::time_point<chrono::steady_clock> end = chrono::steady_clock::now();
    chrono::duration<double> elapsed = end - start;
	// cout << endl << "The time taken to create tree: " << elapsed.count() 
    //     << " seconds"<< endl;
    int ruleClassified = -1;
	vector<int> queryResultSpecific;
	vector<int> queryResultBroadest;
	std::chrono::duration<double> sum_time(0);
	//std::vector<double> times;
	for (int j = 0; j < 10; j++) {
		queryResultSpecific.clear();
		
		start = chrono::steady_clock::now();
		for (int i = 0; i < NUM_PACKETS; i++) {
			if (maxRulesPerNode == 0) {
				classifyLastIterativeBacktrack (mainRoot, packets[i], dim, fieldOrder, ruleClassified);
			}
			  	
			else {
			// 	if (j == 9) {
			// 		auto start1 = chrono::steady_clock::now();
			// 		classifyLastIterativeBacktrackLeafSizeAllRules (mainRoot, packets[i], dim, fieldOrder, ruleClassified);
			// 		auto end1 = chrono::steady_clock::now();
			// 		chrono::duration<double> elapsed1 = end1 - start1;
			// 		times.push_back(elapsed1.count());
			// 	}
			// 	else {
			 		classifyLastIterativeBacktrackLeafSizeAllRules (mainRoot, packets[i], dim, fieldOrder, ruleClassified);
			// 	}
				
			}
			queryResultSpecific.push_back(ruleClassified);
			ruleClassified = -1;
		}
		end = chrono::steady_clock::now();
		elapsed = end - start;
		sum_time += elapsed;
	}
	//cout << endl << "The time taken to classify specific rules for " << numPackets << " packets is: " << sum_time.count() / 10 
      //  << " seconds"<< endl;
	cout << sum_time.count() / 10 << ", ";

	// double mean = calculateMean(times);
    // double stdDev = calculateStdDev(times, mean);

    // // Assuming a 95% confidence interval for a normal distribution
    // double confidenceInterval = 1.96 * stdDev / std::sqrt(times.size()); 

    // std::cout << "Mean duration: " << mean << " ms\n";
    // std::cout << "Standard Deviation: " << stdDev << " ms\n";
    // std::cout << "95% Confidence Interval: +/-" << confidenceInterval << " ms\n";

	std::chrono::duration<double> sum_time2(0);

	for (int j = 0; j < 10; j++) {
		queryResultBroadest.clear();
		start = chrono::steady_clock::now();
		for (int i = 0; i < NUM_PACKETS; i++) {
			classifyFirstTest (mainRoot, packets[i], dim, fieldOrder, ruleClassified);
			queryResultBroadest.push_back(ruleClassified);
			ruleClassified = -1;
		}
		end = chrono::steady_clock::now();
		elapsed = end - start;
		sum_time2 += elapsed;
	}
	//cout << endl << "The time taken to classify broadest rules for " << numPackets << " packets is: " << sum_time2.count() / 10 
       // << " seconds"<< endl;
	cout << sum_time2.count() / 10 ;
	
	//Classifying brute force
	if (brute_force_choice){
		cout <<"\n Classifying Brute force: \n";
		start = chrono::steady_clock::now();
		vector<int> bruteForceResult = classifyBruteForce();
		end = chrono::steady_clock::now();
		elapsed = end - start;
		cout << "\n The time taken to brute force classify " << NUM_PACKETS << " packets is:" << elapsed.count() << " seconds\n";
		std::ofstream outfileBrute ("BruteClassification.txt");
		
		if (!outfileBrute) {
			std::cerr << "Failed to open the file for writing.\n";
			exit (0);
		}
		// Write the vector's contents to the file
		for (const auto& value : bruteForceResult) {
			outfileBrute << value << '\n';  // Writes each element on a new line
		}

		// Close the file (optional here, since the file will close when outfile goes out of scope)
		outfileBrute.close();
	}

	//cleanup.
	for (rule* currRule : rules) {
		delete currRule;
	}
	delete ruleIDs;
	delete mainRoot;
	
}




