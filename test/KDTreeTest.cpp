#include "gtest/gtest.h"

#include <chrono>
#include <random>
#include "algorithm"
#include "../src/kdtree/kdtree.h"
#include "mykdpoint.h"

using namespace std;
typedef kdtree<MyKDPoint<3>, 3> MyKDTree;

class KDTreeTest : public testing::Test {
protected:
	static MyKDTree* mykdtree;

	static void SetUpTestSuite() {
	}

	static void TearDownTestSuite() {
	}

	void SetUp() override
	{
	}

	void TearDown() override
	{
		delete mykdtree;
		mykdtree = nullptr;
	}


private:
	static int getHeightAux(MyKDTree::node* n)
	{
		if (n == nullptr)
			return 0;
		else
			return max(getHeightAux(n->left), getHeightAux(n->right)) + 1;
	}

public:
	static int getHeight(MyKDTree* kd)
	{
		return getHeightAux(kd->get_root());
	}

	template<class vector_t>
	static void shuffle(vector_t& v)
	{
		unsigned seed = (unsigned)std::chrono::system_clock::now().time_since_epoch().count();
		std::shuffle(v.begin(), v.end(), std::default_random_engine(seed));
	}
};

MyKDTree* KDTreeTest::mykdtree = nullptr;

TEST_F(KDTreeTest, BalancedTree) {

	vector<const MyKDPoint<3>*> v1;
	v1.push_back(new MyKDPoint<3>({ -1.857, 1.423, 0.039 }));
	v1.push_back(new MyKDPoint<3>({ 0.209, -1.748, -0.613 }));
	mykdtree = new MyKDTree(v1);
	auto a = mykdtree->range_search(MyKDPoint<3>({ -3.870, 1.255, 0.155 }), 3.5);

	vector<const MyKDPoint<3>*> vec;
	vec.push_back(new MyKDPoint<3>({ 0, 0, 0 }));
	vec.push_back(new MyKDPoint<3>({ 0, 0, 0 }));
	vec.push_back(new MyKDPoint<3>({ 0, 1, 0 }));
	vec.push_back(new MyKDPoint<3>({ 0, 1, 1 }));
	vec.push_back(new MyKDPoint<3>({ 1, 0, 0 }));
	vec.push_back(new MyKDPoint<3>({ 1, 0, 1 }));
	vec.push_back(new MyKDPoint<3>({ 1, 1, 0 }));
	vec.push_back(new MyKDPoint<3>({ 1, 1, 1 }));

	KDTreeTest::shuffle(vec);

	mykdtree = new MyKDTree(vec);
	EXPECT_EQ(KDTreeTest::getHeight(mykdtree), floor(log2(vec.size())) + 1);
}

TEST_F(KDTreeTest, DistanceOnPositiveQuotas) {
	vector<const MyKDPoint<3>*> vec;
	// 3 Node a distance 1 from (0,0,0)
	vec.push_back(new MyKDPoint<3>({ 0, 0, 1 }));
	vec.push_back(new MyKDPoint<3>({ 0, 1, 0 }));
	vec.push_back(new MyKDPoint<3>({ 1, 0, 0 }));

	// 3 Node a distance sqrt(2) from (0,0,0)
	vec.push_back(new MyKDPoint<3>({ 0, 1, 1 }));
	vec.push_back(new MyKDPoint<3>({ 1, 1, 0 }));
	vec.push_back(new MyKDPoint<3>({ 1, 0, 1 }));

	// 1 Node a distance sqrt(3) from (0,0,0)
	vec.push_back(new MyKDPoint<3>({ 1, 1, 1 }));

	KDTreeTest::shuffle(vec);

	mykdtree = new MyKDTree(vec);
	EXPECT_EQ(mykdtree->range_search(MyKDPoint<3>({ 0, 0, 0 }), 1).size(), 3);
	EXPECT_EQ(mykdtree->range_search(MyKDPoint<3>({ 0, 0, 0 }), sqrt(2)).size(), 6);
	EXPECT_EQ(mykdtree->range_search(MyKDPoint<3>({ 0, 0, 0 }), sqrt(3)).size(), 7);
}

TEST_F(KDTreeTest, DistanceOnAllQuotas) {
	vector<const MyKDPoint<3>*> vec;
	// 6 Node a distance 1 from (0,0,0)
	vec.push_back(new MyKDPoint<3>({ 0,  0,  1 }));
	vec.push_back(new MyKDPoint<3>({ 0,  1,  0 }));
	vec.push_back(new MyKDPoint<3>({ 1,  0,  0 }));

	vec.push_back(new MyKDPoint<3>({ 0,  0, -1 }));
	vec.push_back(new MyKDPoint<3>({ 0, -1,  0 }));
	vec.push_back(new MyKDPoint<3>({ -1,  0,  0 }));

	// 12 Node a distance sqrt(2) from (0,0,0)
	vec.push_back(new MyKDPoint<3>({ 0,  1,  1 }));
	vec.push_back(new MyKDPoint<3>({ 1,  1,  0 }));
	vec.push_back(new MyKDPoint<3>({ 1,  0,  1 }));

	vec.push_back(new MyKDPoint<3>({ 0, -1,  1 }));
	vec.push_back(new MyKDPoint<3>({ 0,  1, -1 }));
	vec.push_back(new MyKDPoint<3>({ -1,  1,  0 }));
	vec.push_back(new MyKDPoint<3>({ 1, -1,  0 }));
	vec.push_back(new MyKDPoint<3>({ -1,  0,  1 }));
	vec.push_back(new MyKDPoint<3>({ 1,  0, -1 }));

	vec.push_back(new MyKDPoint<3>({ 0, -1, -1 }));
	vec.push_back(new MyKDPoint<3>({ -1, -1,  0 }));
	vec.push_back(new MyKDPoint<3>({ -1,  0, -1 }));

	// 8 Node a distance sqrt(3) from (0,0,0)
	vec.push_back(new MyKDPoint<3>({ 1,  1,  1 }));

	vec.push_back(new MyKDPoint<3>({ 1,  1, -1 }));
	vec.push_back(new MyKDPoint<3>({ 1, -1,  1 }));
	vec.push_back(new MyKDPoint<3>({ -1,  1,  1 }));

	vec.push_back(new MyKDPoint<3>({ 1, -1, -1 }));
	vec.push_back(new MyKDPoint<3>({ -1, -1,  1 }));
	vec.push_back(new MyKDPoint<3>({ -1,  1, -1 }));

	vec.push_back(new MyKDPoint<3>({ -1, -1, -1 }));

	KDTreeTest::shuffle(vec);

	mykdtree = new MyKDTree(vec);
	EXPECT_EQ(mykdtree->range_search(MyKDPoint<3>({ 0, 0, 0 }), 1).size(), 6);
	EXPECT_EQ(mykdtree->range_search(MyKDPoint<3>({ 0, 0, 0 }), sqrt(2)).size(), 18);
	EXPECT_EQ(mykdtree->range_search(MyKDPoint<3>({ 0, 0, 0 }), sqrt(3)).size(), 26);
}