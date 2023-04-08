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
    node* rihgtChild;

    public: 
        node();
        node(int lR, int rR);
        void setRule(rule& newRule);
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

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
int main()
{
    int numRanges = 0;
    cin >> numRanges;
    set<int> allNums;
    vector<rule> allRules;

    int leftRange, rightRange;
    for(int i = 0; i < numRanges; i++)
    {
        cin >> leftRange >> rightRange;
        allNums.insert(leftRange);
        allNums.insert(rightRange);
        rule tempRule = rule(leftRange, rightRange, i);
        allRules.push_back(tempRule);
    }

    printSet(allNums);
    printRules(allRules);

    return 0;
}
