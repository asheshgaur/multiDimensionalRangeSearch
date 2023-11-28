#ifndef RULES_PACKETS
#define RULES_PACKETS
#include <vector>
#include <utility>
#include <iostream>
typedef std::vector<uint32_t> Packet;
#define MAX_DIMENSIONS 5

const std::vector<int> fieldOrder = {2, 4, 0, 3, 1};

// rule class to keep track of the input rules
class rule
{
    public:
        unsigned int ruleNumber;
        std::vector<unsigned int> prefix_length; //The CIDR prefix for SIP and DIP.
        std::vector<std::pair<unsigned int, unsigned int>> fields;
        int numDimensions;
		// constructors
		rule(int dimensions, int ruleCounter);
		rule(std::vector<std::pair<unsigned int,unsigned int>> fields, int num, int dimensions);
		~rule();
		// display
		void printRule();
		int getLR(int dimension);
		int getRR(int dimension);
		unsigned int getRuleNumber();
		int getDimensionality();
        unsigned int getPrefixLength(int dim);
        unsigned int getDimRangeMin(int dimension);
		unsigned int getDimRangeMax(int dimension);
		void setPrefixLength(int dim, unsigned int value);
		void setDimRangeMin(int dimension, unsigned int value);
		void setDimRangeMax(int dimension, unsigned int value);
		bool satisfiesPacket(const Packet& p);
		bool satisfiesPacketDim(const Packet& p, int dim);
};
#endif