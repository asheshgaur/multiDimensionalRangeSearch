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
		rule();
		rule(int lR, int rR, int num);
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
		node();
		node(int lR, int rR);
		void printNode(); 
		void setRule(rule &newRule);
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

// display method for node
void node::printNode()
{
		cout << "Node " << ": " << lRange << " - " << rRange << endl;
}

void node::setRule(rule& newRule)
{
	rulesAssigned.push_back(&newRule); // push pointer to the rule
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void printSet(set<int> vec)
{
	cout << "Printing: " << endl;
	for(auto a : vec)
		cout << a << " ";
	cout << endl;
}

void printRules(vector<rule> vec)
{
	cout << "Printing Rules: " << endl;
	for(auto a : vec)
		a.printRule();
}

void printNodes(vector<node> vec)
{
	cout << "Printing Nodes: " << endl;
	for(auto a : vec)
		a.printNode();
}

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main()
{
	int numRanges = 0;
	cin >> numRanges;
	set<int> allNums;
	vector<rule> allRules;
	vector<node> allNodes;

	int leftRange, rightRange;
	for(int i = 0; i < numRanges; i++)
	{
		cin >> leftRange >> rightRange;
		allNums.insert(leftRange);
		allNums.insert(rightRange);
		rule tempRule = rule(leftRange, rightRange, i);
		allRules.push_back(tempRule);
	}

	// n0: 10:11, n1: 11:12, n2:
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

	return 0;
}
