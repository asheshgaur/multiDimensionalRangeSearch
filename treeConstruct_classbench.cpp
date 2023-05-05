// create a tree with given ranges

#include <iostream>
#include <vector>
#include <algorithm>
#include <set>
#include <stack>
#include <string>
#include <fstream>
#include <sstream>
#include <iterator>
#include <functional>
#include <regex>

using namespace std;

// rule class to keep track of the input rules
class rule
{
	int ruleNumber;
    unsigned int prefixLength; //The CIDR prefix for SIP and DIP.
	array<array<unsigned int,2>, 5> fields;
	int lRange;
	int rRange;
	int numDimensions;

	public:
		// constructors
		rule(int dimensions);
		rule(array<array<unsigned int,2>, 5> fields, int num, int dimensions);

		// display
		void printRule();
		int getLR(int dimension);
		int getRR(int dimension);
		int getDimensionality();
        unsigned int getPrefixLength();
        array<unsigned int,2> getDimRange(int dimension);
};

// default constructor
rule::rule(int dimensions)
{
	ruleNumber = -1;
	lRange = -1;
	rRange = -1;
    this->numDimensions = dimensions;
}

// non-default constructor
rule::rule(array<array<unsigned int,2>,5>field, int num, int dimensions)
{
	ruleNumber = num + 1;
	this->fields = field;
	this->numDimensions = dimensions;
}

// display method for rule
void rule::printRule()
{
	cout << "Rule " << ruleNumber << ": " << endl;
	for(int i = 0; i < numDimensions; i++)
		cout << "Field" << i+1 << ": " << fields[i][0] << " " << fields[i][1] << endl;
}

int rule::getLR(int dimension){
	return fields[dimension - 1][0];
}

int rule::getRR(int dimension){
	return fields[dimension - 1][1];
}

int rule::getDimensionality(){
	return this->numDimensions;
}
unsigned int rule::getPrefixLength(){
    return this->prefixLength;
}

array<unsigned int,2> rule::getDimRange(int dimension){
    return fields[dimension];
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// node class to keep track of the intermediate nodes in the tree
class node
{
	int lRange;
	int rRange;
	vector<rule*> rulesAssigned;
	node* leftChild;
	node* rightChild;
	node* midChild;

	public: 
		// constructors
		node();
		node(int lR, int rR);

		// display
		void printNode(); 
		
		// setters
		void setRule(rule &newRule);
		void setLeftChild(node* l);
		void setRightChild(node* r);
		void setMidChild(node* m);

		//getters
		int getLR();
		int getRR();
		node* getLC();
		node* getRC();
		node* getMC();
		vector<rule*> getRulesAssigned();
};

// default constructor
node::node()
{
		lRange = -1;
		rRange = -1;
}

// non-default constructor
node::node(int lR, int rR)
{
		lRange = lR;
		rRange = rR;
		leftChild = (node*)NULL;
		rightChild = (node*)NULL;
		midChild = (node*)NULL;
}

// getter for left range
int node::getLR()
{
	return lRange;
}

// getter for right range
int node::getRR()
{
	return rRange;
}

// getter for left child
node* node::getLC()
{
	return leftChild;
}

// getter for right range
node* node::getRC()
{
	return rightChild;
}

node* node::getMC(){
	return midChild;
}

vector<rule*> node::getRulesAssigned() {
	return (this->rulesAssigned);
}

// display method for node
void node::printNode()
{
	cout << "(" << lRange << " - " << rRange << ")" << endl;
	if (this->rulesAssigned.size()!= 0)
	{
		cout << "Rules for this node:" << endl;
		for (int i = 0; i < this->rulesAssigned.size(); i++){
			(this->rulesAssigned[i])->printRule();
			cout << " ";
		}
	}
	cout << endl;
}

// add rule into a specific node
void node::setRule(rule& newRule)
{
	rulesAssigned.push_back(&newRule); // push pointer to the rule
}

// set left child
void node::setLeftChild(node* l)
{
	leftChild = l;
}

// set right child
void node::setRightChild(node* r)
{
	rightChild = r;
}

void node::setMidChild(node* m) {
	this->midChild = m;
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FUNCTIONS not part of class

// create the tree and return the root node of the tree
node* createTree(vector<node*>& allNodes, int l, int r)
{
	// Base Case
	if (r == l)
		return allNodes[r];

	// left recursion
	node* left = createTree(allNodes, l, l + ((r - l) / 2));

	// create the root node
	// sending the left and the right range of this subset
	node *root = new node(allNodes[l]->getLR(), allNodes[r]->getRR()); 
	root->setLeftChild(left); // set left child

	// right recursion - set right child
	root->setRightChild(createTree(allNodes, l + ((r - l) / 2) + 1, r));

	return root; // return root
}

bool isPartiallyOverlapping(rule* r, node* n, int dim){
	if (r->getLR(dim) == n->getLR() && r->getRR(dim) == n->getRR()){
		return false;
	}
	else if (r->getLR(dim) >= n->getRR()){
		return false;
	}
	else if (r->getRR(dim) <= n->getLR()){
		return false;
	}
	else return true;
}
// populate one rule in the tree - left subtree
void populateSubTree(node* rNode, rule* rCurrent, int dim)
{
	if (rNode == NULL)
		return;

	//Check if the node range and the rule range are the same.
	if (rNode->getLR() == rCurrent->getLR(dim) && rNode->getRR() == rCurrent->getRR(dim)){
		rNode->setRule(*rCurrent);
		return;
	}

	// Check if the node range is a subset of the rule range
	else if (rNode->getLR() >= rCurrent->getLR(dim) && rNode->getRR() <= rCurrent->getRR(dim)){
		rNode->setRule(*rCurrent);
		return;
	}
		

	// Check if the node range is partially overlapping with the rule range
	else if (isPartiallyOverlapping(rCurrent, rNode, dim) ){
		populateSubTree(rNode->getLC(), rCurrent, dim);
		populateSubTree(rNode->getRC(), rCurrent, dim);
	}
	else return;
}

bool splitNodeCondition(rule* r, node* rootNode, int dim){
	if (r->getLR(dim) >= rootNode->getLR() && 
		r->getRR(dim) <= rootNode->getRR()){
			return true;
		}
	else return false;
}
node* determineSplitNode(rule* r, node* rootNode, int dim){
	/*
	* Check if the range of the rule is fully inside the
	* range of the current node. 
	*/
	if (splitNodeCondition(r, rootNode, dim)){
			node* splitNode = rootNode;
			if( rootNode->getLC() != (node*)NULL && 
				splitNodeCondition(r, rootNode->getLC(), dim))
				splitNode = determineSplitNode(r, rootNode->getLC(), dim);
			else if (rootNode->getLC() != (node*)NULL && 
				splitNodeCondition(r, rootNode->getRC(), dim))
				splitNode = determineSplitNode(r, rootNode->getRC(), dim);

			else return splitNode;
	}
	else {
		return (node*)NULL;
	}
}
void preOrder(node* n)
{
	if (n == (node*)NULL)
		return;
	n->printNode(); // print node
	preOrder(n->getLC()); // call left child
	preOrder(n->getMC()); // call mid child
	preOrder(n->getRC()); // call right child
}
bool isLeafNode(node* n){
	if(n->getLC() == (node*)NULL && n->getRC() == (node*)NULL)
		return true;
	else return false;
}

// bool isCompleteMatch()
// populate the tree with the rootNode with the rules
void preOrderForNextDim(node* rNode, int nextDim);

//Populate rules to a tree in a given dimension.
void populateRules(node* rNode, vector<rule*> allRules, int dim)
{
    for(int i = 0; i < allRules.size(); i++)
    {
		node* splitNode = determineSplitNode(allRules[i], rNode, dim);

		// Probably don't need this clause
		if(isLeafNode(splitNode)){
			splitNode->setRule(*allRules[i]); // check if this works.
			continue;
		}
		// check if the split node exactly matches the rule range.
		else if(splitNode->getLR() == allRules[i]->getLR(dim) && 
				splitNode->getRR() == allRules[i]->getRR(dim)) {
					splitNode->setRule(*allRules[i]); // check if this works.
					continue;
		}
		else {
			populateSubTree(splitNode->getLC(), allRules[i], dim);
        	populateSubTree(splitNode ->getRC(), allRules[i], dim);
		}
        
    }
	if (allRules[0]->getDimensionality() - dim >= 1) {

		//creating the further dimension trees for the nodes with multiple rules
		preOrderForNextDim(rNode, (dim + 1));
	}
	
}

// Preorder traversal through the populated tree in a given dimension for creating 
// next dimension trees.
void preOrderForNextDim(node* rNode, int nextDim){
	if (rNode == (node*)NULL)
		return;
	else if (rNode->getRulesAssigned().size() > 1) {
		//determine the endpoints for the rules at the node
		vector<rule*> nodeRules = rNode->getRulesAssigned();
		set<int> ruleNums;
		for (int i = 0; i < nodeRules.size(); i++) {
			ruleNums.insert(nodeRules[i]->getLR(nextDim));
			ruleNums.insert(nodeRules[i]->getRR(nextDim));
		}
		// creating the leaf nodes
		// n0: 10:11, n1: 11:12, n2: ... 
		vector<node*> allNodes;
		for(int i = 0; i < ruleNums.size() - 1; i++)
		{
			auto it = next(ruleNums.begin(), i);
			auto it2 = next(ruleNums.begin(), i+1);
			node* tempNode = new node(*it, *it2);
			allNodes.push_back(tempNode);
		}
		node* nextDimNode = createTree(allNodes, 0, allNodes.size() - 1);
		rNode->setMidChild(nextDimNode);
		cout << endl << "Printing dimension " << nextDim << "tree at node: "<< endl;
		rNode->printNode();

		preOrder(nextDimNode);

		// populate the rules in the tree
		populateRules(nextDimNode, nodeRules, nextDim);
		cout << endl << "Printing nodes in dimension " << nextDim <<  "after rule insertion:"<< endl;
		preOrder(nextDimNode);
	}

	preOrderForNextDim(rNode->getLC(), nextDim);
	preOrderForNextDim(rNode->getRC(), nextDim);
}

// function to print set
void printSet(set<int>& vec)
{
	cout << "Printing: " << endl;
	for(auto a : vec)
		cout << a << " ";
	cout << endl;
}

// function to print rules
void printRules(vector<rule*>& vec)
{
	cout << "Printing Rules: " << endl;
	for(auto a : vec)
		a->printRule();
}

// function to print nodes
void printNodes(vector<node*>& vec)
{
	cout << "Printing Leaf Nodes: " << endl;
	for(auto a : vec)
		a->printNode();
}

unsigned int inline atoui(const string& in) {
	std::istringstream reader(in);
	unsigned int val;
	reader >> val;
	return val;
}

std::vector<std::string>& split(const std::string &s, char delim, std::vector<std::string> &elems) {
	std::stringstream ss(s);
	std::string item;
	while (std::getline(ss, item, delim)) {
		elems.push_back(item);
	}
	return elems;
}

std::vector<std::string> split(const std::string &s, char delim) {
	std::vector<std::string> elems;
	split(s, delim, elems);
	return elems;
}

void ReadIPRange(std::array<unsigned int,2>& IPrange,  unsigned int& prefix_length, const string& token)
{
	//cout << token << endl;
	//split slash
	vector<string> split_slash = split(token, '/');
	vector<string> split_ip = split(split_slash[0], '.');
	/*asindmemacces IPv4 prefixes*/
	/*temporary variables to store IP range */
	unsigned int mask;
	int masklit1;
	unsigned int masklit2, masklit3;
	unsigned int ptrange[4];
	for (int i = 0; i < 4; i++)
		ptrange[i] = atoui(split_ip[i]);
	mask = atoui(split_slash[1]);
	
	prefix_length = mask;

	mask = 32 - mask;
	masklit1 = mask / 8;
	masklit2 = mask % 8;

	/*count the start IP */
	for (int i = 3; i>3 - masklit1; i--)
		ptrange[i] = 0;
	if (masklit2 != 0){
		masklit3 = 1;
		masklit3 <<= masklit2;
		masklit3 -= 1;
		masklit3 = ~masklit3;
		ptrange[3 - masklit1] &= masklit3;
	}
	/*store start IP */
	IPrange[0] = ptrange[0];
	IPrange[0] <<= 8;
	IPrange[0] += ptrange[1];
	IPrange[0] <<= 8;
	IPrange[0] += ptrange[2];
	IPrange[0] <<= 8;
	IPrange[0] += ptrange[3];

	//key += std::bitset<32>(IPrange[0] >> prefix_length).to_string().substr(32 - prefix_length);
	/*count the end IP*/
	for (int i = 3; i>3 - masklit1; i--)
		ptrange[i] = 255;
	if (masklit2 != 0){
		masklit3 = 1;
		masklit3 <<= masklit2;
		masklit3 -= 1;
		ptrange[3 - masklit1] |= masklit3;
	}
	/*store end IP*/
	IPrange[1] = ptrange[0];
	IPrange[1] <<= 8;
	IPrange[1] += ptrange[1];
	IPrange[1] <<= 8;
	IPrange[1] += ptrange[2];
	IPrange[1] <<= 8;
	IPrange[1] += ptrange[3];
}

int ReadFilter(vector<string>& tokens, vector<rule*>& ruleset, unsigned int cost)
{
	// 5 fields: sip, dip, sport, dport, proto = 0 (with@), 1, 2 : 4, 5 : 7, 8

	/*allocate a few more bytes just to be on the safe side to avoid overflow etc*/
    int dim = 5;
    int reps = 1;
	rule temp_rule =  rule(dim);
	string key;
	if (tokens[0].at(0) != '@')  {
		/* each rule should begin with an '@' */
		printf("ERROR: NOT A VALID RULE FORMAT\n");
		exit(1);
	}

	int index_token = 0;
	int i = 0;
	for (int rep = 0; rep < reps; rep++)
	{
		/* reading SIP range */
		if (i == 0) {
            
			ReadIPRange(temp_rule.getDimRange(i), temp_rule.prefix_length[i], tokens[index_token++].substr(1));
			i++;
		} else {
			ReadIPRange(temp_rule.range[i], temp_rule.prefix_length[i], tokens[index_token++]);
			i++;
		}
		/* reading DIP range */
		ReadIPRange(temp_rule.range[i], temp_rule.prefix_length[i], tokens[index_token++]);
		i++;
		ReadPort(temp_rule.range[i++], tokens[index_token], tokens[index_token + 2]);
		index_token += 3;
		ReadPort(temp_rule.range[i++], tokens[index_token], tokens[index_token + 2]);
		index_token += 3;
		ReadProtocol(temp_rule.range[i++], tokens[index_token++]);
	}
    // Not working with priority right now.
	//temp_rule.priority = cost;

	ruleset.push_back(temp_rule);

	return 0;
}

//taken from github.com/sorrachai/PartitionSort
void LoadFilters(ifstream& fp, vector<rule*>& ruleset)
{
	int line_number = 0;
	string content;
	while (getline(fp, content)) {
		istringstream iss(content);
		vector<string> tokens{ istream_iterator < string > {iss}, istream_iterator < string > {} };
		ReadFilter(tokens, ruleset, line_number++);
	}
}

//taken from github.com/sorrachai/PartitionSort
vector<rule*> readFilterFileClassBench(const string&  filename)
{
	//assume 5*rep fields

	vector<rule*> rules;
	std::ifstream column_counter(filename);
	ifstream input_file(filename);
	if (!input_file.is_open() || !column_counter.is_open())
	{
		printf("Couldnt open filter set file \n");
		exit(1);
	}


	LoadFilters(input_file, rules);
	input_file.close();
	column_counter.close();

	return	rules;
}
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main function
int main()
{
    // Reading of a file will be completely different due to classbench format.
    // Using code from PartitionSort to read files.
	// read in the number of rules
	int numRules = 0;
	int numDimensions = 0;
	int dim = 1;
	cin >> numRules;
	cin >> numDimensions;

	
	set<int> allNums; // set to store all the numbers in 1 dim - sorted without duplicates
	vector<rule*> allRules; // array of rule pointers
	vector<node*> allNodes; // arra of leaf node pointers for 1st dim

	// read in all the rules
	int leftRange, rightRange;
	vector<pair<int,int>>tempRanges;
	for(int i = 0; i < numRules; i++)
	{
		for (int j = 0; j < numDimensions; j++){
			cin >> leftRange >> rightRange;

			//inserting range endpoints for the first dimension only.
			if (j == 0){
				allNums.insert(leftRange);
				allNums.insert(rightRange);
			}
			
			tempRanges.push_back(make_pair(leftRange, rightRange));
		}
		// temp rule to insert into allRules
		rule* tempRule = new rule(tempRanges, i, numDimensions); // calling the constructor

		tempRanges.clear();
		allRules.push_back(tempRule);
	}

	// creating the leaf nodes
	// n0: 10:11, n1: 11:12, n2: ... 
	for(int i = 0; i < allNums.size() - 1; i++)
	{
		auto it = next(allNums.begin(), i);
		auto it2 = next(allNums.begin(), i+1);
		node* tempNode = new node(*it, *it2);
		allNodes.push_back(tempNode);
	}

	// printNodes(allNodes1);
	// printSet(allNums1);
	// printSet(allNums2);
	// printRules(allRules);

	// create the tree
	node* rootNode = createTree(allNodes, 0, allNodes.size() - 1);

	cout << endl << "Printing nodes before rule insertion in main function:"<< endl;

	preOrder(rootNode);

    // populate the rules in the tree
    populateRules(rootNode, allRules, dim);
	cout << endl << "Printing nodes after rule insertion in main function:"<< endl;
	preOrder(rootNode);

    // destructors for all the local vector of pointers

	return 0;
}