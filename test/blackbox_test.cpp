#include "blackbox_test.h"

#include "../app/sources/private/rin_graph_impl.h"

#include "../app/sources/private/chemical_entity_impl.h"

#define MY_DBL_EPSILON 0.00005

using namespace std;
using namespace rin;
namespace fs = std::filesystem;

static fs::path running_folder;
static fs::path running_path;
void set_running_path(fs::path exe_path) {
    running_path = exe_path;
    running_folder = exe_path.remove_filename();
}

bool compare(const double& a, const double& b) { return abs(a - b) < MY_DBL_EPSILON; }

string interaction_name(const edge& e)
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

bool compare_distance(const edge& e, const double& expected) { return compare(stod(e.distance()), expected); }
bool compare_angle(const edge& e, const double& expected) { return compare(stod(e.angle()), expected); }
bool compare_energy(const edge& e, const double& expected) { return compare(stod(e.energy()), expected); }


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
    vector<edge> find_edges(const function<bool(const edge&)>& pred) const { return find(edges, pred); }
    unsigned int count_edges(const function<bool(const edge&)>& pred) const { return (unsigned int)find_edges(pred).size(); }
    bool contain_edge(const function<bool(const edge&)>& pred) const { return count_edges(pred) > 0; }
};

class BlackBoxTest : public testing::Test {

private:
    inline static const fs::path test_case_folder = "test_case/";

protected:
    static void SetUpTestSuite() { }

    static void TearDownTestSuite() { }

    Result SetUp(const string& filename)
    {
        string exePath = running_path.string();
        string pdbPath = (running_folder / test_case_folder / filename).string();
        const char* args[2] = { exePath.c_str(), pdbPath.c_str() };

        std::optional<arguments> maybe_args;
        maybe_args = read_args(2, args);
        // won't throw std::bad_optional because args are hand crafted above
        auto const parsed_args = maybe_args.value();
        auto const models = rin::maker::parse_models(parsed_args.pdb_path);
        // again, won't throw because tests are hand crafted to have 1 model each
        return Result(*models.at(0), parsed_args.params);
    }

    void TearDown() override { }
};

#pragma region IonIon


TEST_F(BlackBoxTest, IonIon1) {
    Result r = SetUp("ionion/1.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 0);
}
TEST_F(BlackBoxTest, IonIon2) {
    Result r = SetUp("ionion/2.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 2);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "HIS" && target(e) == "ASP" &&
               compare_distance(e, 1.21934) &&
               source_atom(e) == "CD2:CE1:CG:ND1:NE2" && target_atom(e) == "CG:OD1:OD2";
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               compare_distance(e, 1.68740) &&
               source_atom(e) == "NZ" && target_atom(e) == "CD:OE1:OE2";
    }));
}
TEST_F(BlackBoxTest, IonIon3) {
    Result r = SetUp("ionion/3.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 3);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               compare_distance(e, 3.22121) &&
               source_atom(e) == "NZ" && target_atom(e) == "CD:OE1:OE2";
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "HIS" && target(e) == "GLU" &&
               compare_distance(e, 1.60990) &&
               source_atom(e) == "CD2:CE1:CG:ND1:NE2" && target_atom(e) == "CD:OE1:OE2";
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "ASP" &&
               compare_distance(e, 2.55251) &&
               source_atom(e) == "NZ" && target_atom(e) == "CG:OD1:OD2";
    }));
}
TEST_F(BlackBoxTest, IonIon4) {
    Result r = SetUp("ionion/4.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               compare_distance(e, 2.56018) &&
               source_atom(e) == "NZ" && target_atom(e) == "OE1";
    }));
}
TEST_F(BlackBoxTest, IonIon5) {
    Result r = SetUp("ionion/5.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 2);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               compare_distance(e, 1.95108) &&
               source_atom(e) == "NZ" && target_atom(e) == "OE1";
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "ARG" && target(e) == "GLU" &&
               compare_distance(e, 2.03282) &&
               source_atom(e) == "NH2" && target_atom(e) == "OE1";
    }));
}
TEST_F(BlackBoxTest, IonIon6) {
    Result r = SetUp("ionion/6.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 0);
}
TEST_F(BlackBoxTest, IonIon7) {
    Result r = SetUp("ionion/7.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 0);
}
TEST_F(BlackBoxTest, IonIon8) {
    Result r = SetUp("ionion/8.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "ASP" &&
               compare_distance(e, 3.19828) &&
               source_atom(e) == "NZ" && target_atom(e) == "CG:OD1:OD2";
    }));
}
TEST_F(BlackBoxTest, IonIon9) {
    Result r = SetUp("ionion/9.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "IONIC"; }), 0);
}

#pragma endregion

#pragma region HBond

TEST_F(BlackBoxTest, HBond1) {
    Result r = SetUp("hbond/1.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 0);
}
TEST_F(BlackBoxTest, HBond2) {
    Result r = SetUp("hbond/2.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 0);
}
TEST_F(BlackBoxTest, HBond3) {
    Result r = SetUp("hbond/3.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "ASP" && target(e) == "ARG" && source_atom(e) == "OD1" &&
               target_atom(e) == "NH2" &&
               compare_distance(e, 2.11193) && compare_angle(e, 177.61115);
    }));
}
TEST_F(BlackBoxTest, HBond4) {
    Result r = SetUp("hbond/4.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "GLN" && target(e) == "LYS" && source_atom(e) == "OE1" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.53981) && compare_angle(e, 179.91626);
    }));
}
TEST_F(BlackBoxTest, HBond5) {
    Result r = SetUp("hbond/5.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 2);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "ASN" && target(e) == "ASN" && source_atom(e) == "OD1" &&
               target_atom(e) == "ND2" &&
               compare_distance(e, 2.01970) && compare_angle(e, 179.96223);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "ASN" && target(e) == "ASN" && source_atom(e) == "OD1" &&
               target_atom(e) == "ND2" &&
               compare_distance(e, 2.02333) && compare_angle(e, 179.96160);
    }));
}
TEST_F(BlackBoxTest, HBond6) {
    Result r = SetUp("hbond/6.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 3);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" && source_atom(e) == "OG" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.02752) && compare_angle(e, 179.97440);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" && source_atom(e) == "OG" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.00981) && compare_angle(e, 179.94896);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" && source_atom(e) == "OG" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.02536) && compare_angle(e, 179.95765);
    }));
}
TEST_F(BlackBoxTest, HBond7) {
    Result r = SetUp("hbond/7.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" && source_atom(e) == "OG" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.02752) && compare_angle(e, 179.97440);
    }));
}
TEST_F(BlackBoxTest, HBond8) {
    Result r = SetUp("hbond/8.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "HBOND"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "GLN" && target(e) == "LYS" && source_atom(e) == "OE1" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.56130) && compare_angle(e, 67.88154);
    }));
}

#pragma endregion

#pragma region PiPi

TEST_F(BlackBoxTest, PiPi1) {
    Result r = SetUp("pipi/1.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi2) {
    Result r = SetUp("pipi/2.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi3) {
    Result r = SetUp("pipi/3.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi4) {
    Result r = SetUp("pipi/4.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PIPISTACK"; }), 0);
}
TEST_F(BlackBoxTest, PiPi5) {
    Result r = SetUp("pipi/5.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PIPISTACK"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PIPISTACK" && source(e) == "PHE" && target(e) == "PHE" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" &&
               compare_distance(e, 1) && compare_angle(e, 0);
    }));
}

TEST_F(BlackBoxTest, PiPi6) {
    Result r = SetUp("pipi/6.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PIPISTACK"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PIPISTACK" && source(e) == "PHE" && target(e) == "PHE" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" &&
               compare_distance(e, 1.99982) && compare_angle(e, 0);
    }));
}

#pragma endregion

#pragma region VDW

TEST_F(BlackBoxTest, vdw1) {
    Result r = SetUp("vdw/1.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 1);
}
TEST_F(BlackBoxTest, vdw2) {
    Result r = SetUp("vdw/2.pdb");
    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw3) {
    Result r = SetUp("vdw/3.pdb");
    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw4) {
    Result r = SetUp("vdw/4.pdb");
    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw5) {
    Result r = SetUp("vdw/5.pdb");
    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw6) {
    Result r = SetUp("vdw/6.pdb");
    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 0);
}
TEST_F(BlackBoxTest, vdw7) {
    Result r = SetUp("vdw/7.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "VDW" && target(e) == "ASN" && source(e) == "GLN" &&
                target_atom(e) == "CB" && source_atom(e) == "NE2" &&
               compare_distance(e, 3.85866) && compare_energy(e, -0.13357);
    }));
}
TEST_F(BlackBoxTest, vdw8) {
    Result r = SetUp("vdw/8.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "VDW"; }), 4);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "ASN" && target(e) == "ASN" &&
               source_atom(e) == "ND2" && target_atom(e) == "CB" &&
               compare_distance(e, 3.74499) && compare_energy(e, -0.10873);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "GLN" && target(e) == "ASN" &&
               source_atom(e) == "OE1" && target_atom(e) == "CB" &&
               compare_distance(e, 1.73218) && compare_energy(e, 2022.14052);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "GLN" && target(e) == "GLN" &&
               source_atom(e) == "OE1" && target_atom(e) == "NE2" &&
               compare_distance(e, 3.49130) && compare_energy(e, -0.18889);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "ASN" && target(e) == "GLN" &&
               source_atom(e) == "CB" && target_atom(e) == "NE2" &&
               compare_distance(e, 1.75928) && compare_energy(e, 2653.90001);
    }));
}

#pragma endregion

#pragma region Pication

TEST_F(BlackBoxTest, picat1) {
    Result r = SetUp("picat/1.pdb");


    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PICATION"; }), 0);
}
TEST_F(BlackBoxTest, picat2) {
    Result r = SetUp("picat/2.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PICATION"; }), 4);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 3.99983) && compare_angle(e, 60.00134);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 4) && compare_angle(e, 89.99761);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 4.00013) && compare_angle(e, 70.00445);
    }));
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 4) && compare_angle(e, 89.99761);
    }));
}
TEST_F(BlackBoxTest, picat3) {
    Result r = SetUp("picat/3.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PICATION"; }), 0);
}
TEST_F(BlackBoxTest, picat4) {
    Result r = SetUp("picat/4.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PICATION"; }), 0);
}
TEST_F(BlackBoxTest, picat5) {
    Result r = SetUp("picat/5.pdb");

    EXPECT_EQ(r.count_edges([](const edge &e) { return interaction_name(e) == "PICATION"; }), 1);
    EXPECT_TRUE(r.contain_edge([](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 2.99985) && compare_angle(e, 70.00357);
    }));
}

#pragma endregion