tree:						treeConstruct.cpp
							g++ -g -o tree treeConstruct.cpp
multitree:					treeConstruct_multiDim.cpp
							g++ -g -o multiTree treeConstruct_multiDim.cpp

multiTreeClassbench:		treeConstruct_classbench.cpp
							g++ -g -o multiTreeClassbench treeConstruct_classbench.cpp

multiTreeSridhar:			treeConstruct_Sridhar.cpp
							g++ -g -o multiTreeSridhar treeConstruct_Sridhar.cpp

multiTreePrunedSridhar:		treeConstructPruned_Sridhar.cpp
							g++ -g -o multiTreePrunedSridhar treeConstructPruned_Sridhar.cpp

rangeTree:					range-v3.cpp
							g++ -g -o rangeTree range-v3.cpp

fusionClassify:				fusion_classify.cpp
							g++ -g -o fusionClassify fusion_classify.cpp

simulation:					Simulation_classbench-files.cpp
							g++ -o simulation Simulation_classbench-files.cpp -std=c++14 -O3

simulation_copy:			Simulation_classbench-files-copy.cpp
							g++ -o simulation_copy Simulation_classbench-files-copy.cpp -std=c++14 -O3 -Wall 

simulation_copy-debug:		Simulation_classbench-files-copy.cpp
							g++ -o simulation_copy-debug Simulation_classbench-files-copy.cpp -std=c++14 -g

analytics:					ruleAnalytics.cpp
							g++ -o analytics ruleAnalytics.cpp

percentage:					percentageRulesWrong.cpp
							g++ -o percentage percentageRulesWrong.cpp -O3

clean:
		rm -f tree
		rm -f multiTree
		rm -f multiTreeClassbench
		rm -f multiTreeSridhar
		rm -f multiTreePrunedSridhar
		rm -f rangeTree
		rm -f fusionClassify
		rm -f simulation
		rm -f percentage
		rm -f simulation_copy
		rm -f simulation_copy-debug