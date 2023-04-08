// create a tree with given ranges

#include <iostream>
#include <vector>
#include <algorithm>
#include <set>

using namespace std;

// rule class to keep track of the input rules
class rule
{
	int ruleNumber;
	int lRange;
	int rRange;

	public:
		// constructors
		rule();
		rule(int lR, int rR, int num);

		// display
		void printRule();
};

// default constructor
rule::rule()
{
	ruleNumber = -1;
	lRange = -1;
	rRange = -1;
}

// non-default constructor
rule::rule(int lR, int rR, int num)
{
	ruleNumber = num;
	lRange = lR;
	rRange = rR;
}

// display method for rule
void rule::printRule()
{
	cout << "Rule " << ruleNumber << ": " << lRange << " - " << rRange << endl;
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

		//getters
		int getLR();
		int getRR();
		node* getLC();
		node* getRC();
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
		leftChild = NULL;
		rightChild = NULL;
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

// display method for node
void node::printNode()
{
	cout << "(" << lRange << " - " << rRange << ")" << endl;
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// FUNCTIONS not part of class

// create the tree and return the root node of the tree
node* createTree(vector<node> allNodes, int l, int r)
{
	// Base Case - needs to change??
	if ((l-r) <= 0)
		return NULL;

	// left recursion
	node* left = createTree(allNodes, l, ((l-r) / 2));

	// create the root node
	// sending the left and the right range of this subset
	node *root = new node(allNodes[l].getLR(), allNodes[r].getRR()); 
	root->setLeftChild(left); // set left child

	// right recursion - set right child
	root->setRightChild(createTree(allNodes, ((l - r) / 2) + 1, r));

	return root; // return root
}

// tree traversal - preorder 
void preOrder(node* n)
{
	if (n == NULL)
		return;
	n->printNode(); // print node
	preOrder(n->getLC()); // call left child
	preOrder(n->getRC()); // call right child
}

// function to print set
void printSet(set<int> vec)
{
	cout << "Printing: " << endl;
	for(auto a : vec)
		cout << a << " ";
	cout << endl;
}

// function to print rules
void printRules(vector<rule> vec)
{
	cout << "Printing Rules: " << endl;
	for(auto a : vec)
		a.printRule();
}

// function to print nodes
void printNodes(vector<node> vec)
{
	cout << "Printing Nodes: " << endl;
	for(auto a : vec)
		a.printNode();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
// Main function
int main()
{
	// read in the number of rules
	int numRanges = 0;
	cin >> numRanges;

	set<int> allNums; // set to store all the numbers - sorted without duplicates
	vector<rule> allRules; // array of rules
	vector<node> allNodes; // arra of leaf nodes

	// read in all the rules
	int leftRange, rightRange;
	for(int i = 0; i < numRanges; i++)
	{
		cin >> leftRange >> rightRange;
		allNums.insert(leftRange);
		allNums.insert(rightRange);

		// temp rule to insert into allRules
		rule tempRule = rule(leftRange, rightRange, i); // calling the constructor
		allRules.push_back(tempRule);
	}

	// creating the leaf nodes
	// n0: 10:11, n1: 11:12, n2: ... 
	for(int i = 0; i < allNums.size() - 1; i++)
	{
		auto it = next(allNums.begin(), i);
		auto it2 = next(allNums.begin(), i+1);
		node tempNode = node(*it, *it2);
		allNodes.push_back(tempNode);
	}

	printNodes(allNodes);
	printSet(allNums);
	printRules(allRules);

	// not working yet
	node* rootNode = createTree(allNodes, 0, allNodes.size() - 1);
	preOrder(rootNode);

	return 0;
}
