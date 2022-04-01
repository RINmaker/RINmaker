#include "blackbox_test.h"

#define MY_DBL_EPSILON 0.00005

using namespace std;
using namespace rin;
namespace fs = std::filesystem;

static fs::path runningFolder;
static fs::path runningPath;
void setRunningPath(fs::path exePath) {
    runningPath = exePath;
    runningFolder = exePath.remove_filename();
}

bool compare(const double& a, const double& b) { return abs(a - b) < MY_DBL_EPSILON; }

string interactionName(const edge& e)
{
    string const& interaction = e.interaction();
    return interaction.substr(0, interaction.find_first_of(':'));
}
string source(const edge& e)
{
    string const& source = e.source_id();
    return source.substr(source.find_last_of(':') + 1, string::npos);
}
string target(const edge& e)
{
    string const& target = e.target_id();
    return target.substr(target.find_last_of(':') + 1, string::npos);
}
string source_atom(const edge& e) { return e.source_atom(); }
string target_atom(const edge& e) { return e.target_atom(); }

bool compareDistance(const edge& e, const double& expected) { return compare(stod(e.distance()), expected); }
bool compareAngle(const edge& e, const double& expected) { return compare(stod(e.angle()), expected); }
bool compareEnergy(const edge& e, const double& expected) { return compare(stod(e.energy()), expected); }


class Result
{
public:
    graph rin_graph;
    std::vector<edge> edges;
    std::unordered_map<std::string, node> nodes;

public:
    explicit Result(rin::maker const& rm, rin::parameters const& params) : rin_graph(rm(params))
    {
        edges = rin_graph.get_edges();
        nodes = rin_graph.get_nodes();
    }

private:
    template<typename T>
    static vector<T> find(vector<T> vec, const function<bool(const T&)>& pred)
    {
        vector<T> output;
        for (const T& e : vec)
        {
            if (pred(e))
                output.push_back(e);
        }
        return output;
    }
public:
    vector<edge> findEdges(const function<bool(const edge&)>& pred) { return find(edges, pred); }
    int countEdges(const function<bool(const edge&)>& pred) { return findEdges(pred).size(); }
    bool containEdge(const function<bool(const edge&)>& pred) { return countEdges(pred) > 0; }
};

class BlackBoxTest : public testing::Test {

private:
    inline static const fs::path testCaseFolder = "test_case/";

protected:
    static void SetUpTestSuite() { }

    static void TearDownTestSuite() { }

    Result SetUp(const string& filename)
    {
        string exePath = runningPath.string();
        string pdbPath = (runningFolder / testCaseFolder / filename).string();
        const char* args[2] = { exePath.c_str(), pdbPath.c_str() };

        std::optional<arguments> maybe_args;
        read_args(2, args, maybe_args);

        // won't throw std::bad_optional because args are hand crafted above
        arguments parsed = maybe_args.value();
        return Result(rin::maker(parsed.pdb_path), parsed.params);
    }

    void TearDown() override { }
};

#pragma region IonIon


TEST_F(BlackBoxTest, IonIon1) {
    Result r = SetUp("ionion/1.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 0);
}
TEST_F(BlackBoxTest, IonIon2) {
    Result r = SetUp("ionion/2.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 2);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "IONIC" && source(e) == "HIS" && target(e) == "ASP" && compareDistance(e, 1.21934) &&
               source_atom(e) == "CD2:CE1:CG:ND1:NE2" && target_atom(e) == "CG:OD1:OD2";
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" && compareDistance(e, 1.68740) &&
               source_atom(e) == "NZ" && target_atom(e) == "CD:OE1:OE2";
    }));
}
TEST_F(BlackBoxTest, IonIon3) {
    Result r = SetUp("ionion/3.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 3);
}
TEST_F(BlackBoxTest, IonIon4) {
    Result r = SetUp("ionion/4.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" && compareDistance(e, 2.56018) &&
               source_atom(e) == "NZ" && target_atom(e) == "OE1";
    }));
}
TEST_F(BlackBoxTest, IonIon5) {
    Result r = SetUp("ionion/5.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 2);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" && compareDistance(e, 1.95108) &&
               source_atom(e) == "NZ" && target_atom(e) == "OE1";
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "IONIC" && source(e) == "ARG" && target(e) == "GLU" && compareDistance(e, 2.03282) &&
               source_atom(e) == "NH2" && target_atom(e) == "OE1";
    }));
}
TEST_F(BlackBoxTest, IonIon6) {
    Result r = SetUp("ionion/6.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 0);
}
TEST_F(BlackBoxTest, IonIon7) {
    Result r = SetUp("ionion/7.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "IONIC"; }), 0);
}

#pragma endregion

#pragma region HBond

TEST_F(BlackBoxTest, HBond1) {
    Result r = SetUp("hbond/1.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 0);
}
TEST_F(BlackBoxTest, HBond2) {
    Result r = SetUp("hbond/2.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 0);
}
TEST_F(BlackBoxTest, HBond3) {
    Result r = SetUp("hbond/3.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "ASP" && target(e) == "ARG" && source_atom(e) == "OD1" && target_atom(e) == "NH2" &&
               compareDistance(e, 2.11193) && compareAngle(e, 177.61115);
    }));
}
TEST_F(BlackBoxTest, HBond4) {
    Result r = SetUp("hbond/4.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "GLN" && target(e) == "LYS" && source_atom(e) == "OE1" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.53981) && compareAngle(e, 179.91626);
    }));
}
TEST_F(BlackBoxTest, HBond5) {
    Result r = SetUp("hbond/5.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 2);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "ASN" && target(e) == "ASN" && source_atom(e) == "OD1" && target_atom(e) == "ND2" &&
               compareDistance(e, 2.01970) && compareAngle(e, 179.96223);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "ASN" && target(e) == "ASN" && source_atom(e) == "OD1" && target_atom(e) == "ND2" &&
               compareDistance(e, 2.02333) && compareAngle(e, 179.96160);
    }));
}
TEST_F(BlackBoxTest, HBond6) {
    Result r = SetUp("hbond/6.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 3);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "TYR" && target(e) == "LYS" && source_atom(e) == "OH" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.02752) && compareAngle(e, 179.97440);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "TYR" && target(e) == "LYS" && source_atom(e) == "OH" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.00981) && compareAngle(e, 179.94896);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "TYR" && target(e) == "LYS" && source_atom(e) == "OH" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.02536) && compareAngle(e, 179.95765);
    }));
}
TEST_F(BlackBoxTest, HBond7) {
    Result r = SetUp("hbond/7.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "TYR" && target(e) == "LYS" && source_atom(e) == "OH" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.02752) && compareAngle(e, 179.97440);
    }));
}
TEST_F(BlackBoxTest, HBond8) {
    Result r = SetUp("hbond/8.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "HBOND" && source(e) == "GLN" && target(e) == "LYS" && source_atom(e) == "OE1" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.56130) && compareAngle(e, 67.88154);
    }));
}

#pragma endregion

#pragma region PiPi

TEST_F(BlackBoxTest, PiPi1) {
    Result r = SetUp("pipi/1.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi2) {
    Result r = SetUp("pipi/2.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi3) {
    Result r = SetUp("pipi/3.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi4) {
    Result r = SetUp("pipi/4.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi5) {
    Result r = SetUp("pipi/5.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PIPISTACK"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PIPISTACK" && source(e) == "PHE" && target(e) == "PHE" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" &&
               compareDistance(e, 1) && compareAngle(e, 0);
    }));
}

TEST_F(BlackBoxTest, PiPi6) {
    Result r = SetUp("pipi/6.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PIPISTACK"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PIPISTACK" && source(e) == "PHE" && target(e) == "PHE" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" &&
               compareDistance(e, 1.99982) && compareAngle(e, 0);
    }));
}

#pragma endregion

#pragma region VDW

TEST_F(BlackBoxTest, vdw1) {
    Result r = SetUp("vdw/1.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 1);
}
TEST_F(BlackBoxTest, vdw2) {
    Result r = SetUp("vdw/2.pdb");
    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw3) {
    Result r = SetUp("vdw/3.pdb");
    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw4) {
    Result r = SetUp("vdw/4.pdb");
    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw5) {
    Result r = SetUp("vdw/5.pdb");
    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw6) {
    Result r = SetUp("vdw/6.pdb");
    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw7) {
    Result r = SetUp("vdw/7.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "VDW" && source(e) == "ASN" && target(e) == "GLN" &&
               source_atom(e) == "CB" && target_atom(e) == "NE2" &&
               compareDistance(e, 3.85866) && compareEnergy(e, -0.03763);
    }));
}
TEST_F(BlackBoxTest, vdw8) {
    Result r = SetUp("vdw/8.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "VDW"; }), 4);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "VDW" && source(e) == "ASN" && target(e) == "ASN" &&
               source_atom(e) == "ND2" && target_atom(e) == "CB" &&
               compareDistance(e, 3.74499) && compareEnergy(e, -0.03269);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "VDW" && source(e) == "GLN" && target(e) == "ASN" &&
               source_atom(e) == "OE1" && target_atom(e) == "CB" &&
               compareDistance(e, 1.73218) && compareEnergy(e, 1523.88608);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "VDW" && source(e) == "GLN" && target(e) == "GLN" &&
               source_atom(e) == "OE1" && target_atom(e) == "NE2" &&
               compareDistance(e, 3.49130) && compareEnergy(e, -0.04878);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "VDW" && source(e) == "ASN" && target(e) == "GLN" &&
               source_atom(e) == "CB" && target_atom(e) == "NE2" &&
               compareDistance(e, 1.75928) && compareEnergy(e, 2036.21611);
    }));
}

#pragma endregion

#pragma region Pication

TEST_F(BlackBoxTest, picat1) {
    Result r = SetUp("picat/1.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PICATION"; }), 0);
}
TEST_F(BlackBoxTest, picat2) {
    Result r = SetUp("picat/2.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PICATION"; }), 4);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compareDistance(e, 3.99983) && compareAngle(e, 60.00134);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compareDistance(e, 4) && compareAngle(e, 89.99761);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compareDistance(e, 4.00013) && compareAngle(e, 70.00445);
    }));
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compareDistance(e, 4) && compareAngle(e, 89.99761);
    }));
}
TEST_F(BlackBoxTest, picat3) {
    Result r = SetUp("picat/3.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PICATION"; }), 0);
}
TEST_F(BlackBoxTest, picat4) {
    Result r = SetUp("picat/4.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PICATION"; }), 0);
}
TEST_F(BlackBoxTest, picat5) {
    Result r = SetUp("picat/5.pdb");

    EXPECT_EQ(r.countEdges([](const edge& e) { return interactionName(e) == "PICATION"; }), 1);
    EXPECT_TRUE(r.containEdge([](const edge& e) {
        return interactionName(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compareDistance(e, 2.99985) && compareAngle(e, 70.00357);
    }));
}

#pragma endregion