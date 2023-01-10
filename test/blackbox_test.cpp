#include "blackbox_test.h"

#include "../app/sources/private/impl_rin_graph.h"

#include "../app/sources/private/impl_chemical_entity.h"

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
    string const& interaction = e.get_interaction();
    return interaction.substr(0, interaction.find_first_of(':'));
}
string source(const edge& e)
{
    string const& source = e.get_source_id();
    return source.substr(source.find_last_of(':') + 1, string::npos);
}
string target(const edge& e)
{
    string const& target = e.get_target_id();
    return target.substr(target.find_last_of(':') + 1, string::npos);
}
string source_atom(const edge& e) { return e.get_source_atom(); }
string target_atom(const edge& e) { return e.get_target_atom(); }

bool compare_distance(const edge& e, const double& expected) { return compare(stod(e.get_distance()), expected); }
bool compare_angle(const edge& e, const double& expected) { return compare(stod(e.get_angle()), expected); }
bool compare_energy(const edge& e, const double& expected) { return compare(stod(e.get_energy()), expected); }

auto isIonicFunc = [](const edge &e) { return interaction_name(e) == "IONIC"; };
auto isHbondFunc = [](const edge &e) { return interaction_name(e) == "HBOND"; };
auto isPipiFunc = [](const edge &e) { return interaction_name(e) == "PIPISTACK"; };
auto isVdwFunc = [](const edge &e) { return interaction_name(e) == "VDW"; };
auto isPicatFunc = [](const edge &e) { return interaction_name(e) == "PICATION"; };

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

    Result SetUp(const string& filename, const vector<const char*>& additionalParameters = {})
    {
        string exePath = running_path.string();
        string pdbPath = (running_folder / test_case_folder / filename).string();

        // RINmaker -i filename -o dummy --illformed=kall rin
        vector<const char*> mandatoryParameters = { exePath.c_str(), "-i", pdbPath.c_str(), "-o dummy", "--illformed=kall", "rin"};

        vector<const char*> parameters;
        parameters.reserve(mandatoryParameters.size() + additionalParameters.size());
        parameters.insert(parameters.end(), mandatoryParameters.begin(), mandatoryParameters.end());
        parameters.insert(parameters.end(), additionalParameters.begin(), additionalParameters.end());

        auto maybe_args = read_args(static_cast<int>(parameters.size()), parameters.data());
        // won't throw std::bad_optional because args are handcrafted above
        auto const parsed_args = maybe_args.value();
        auto protein_structure = gemmi::read_pdb_file(parsed_args.input().string());
        // again, won't throw because tests are handcrafted to have 1 model each
        return Result(rin::maker{protein_structure.first_model(), protein_structure, parsed_args}, parsed_args);
    }

    void TearDown() override { }
};

#pragma region IonIon


TEST_F(BlackBoxTest, IonIon1) {
    Result r = SetUp("ionion/ionion1.pdb");

    EXPECT_EQ(r.count_edges(isIonicFunc), 0);
}
TEST_F(BlackBoxTest, IonIon2) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "HIS" && target(e) == "ASP" &&
               source_atom(e) == "CD2:CE1:CG:ND1:NE2" && target_atom(e) == "CG:OD1:OD2" &&
               compare_distance(e, 1.21934) && compare_energy(e, 3.96073);
    };
    auto e2 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               source_atom(e) == "NZ" && target_atom(e) == "CD:OE1:OE2" &&
               compare_distance(e, 1.68740) && compare_energy(e, 8.05503);
    };

    {
        Result r = SetUp("ionion/ionion2.pdb");

        EXPECT_EQ(r.count_edges(isIonicFunc), 2);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
    }
    {
        Result r = SetUp("ionion/ionion2.pdb", {"--ionic-bond", "1.23"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("ionion/ionion2.pdb", {"--ionic-bond", "1.2"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 0);
    }
}
TEST_F(BlackBoxTest, IonIon3) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "HIS" && target(e) == "GLU" &&
               source_atom(e) == "CD2:CE1:CG:ND1:NE2" && target_atom(e) == "CD:OE1:OE2" &&
               compare_distance(e, 1.60990) && compare_energy(e, 5.01294);
    };
    auto e2 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "ASP" &&
               source_atom(e) == "NZ" && target_atom(e) == "CG:OD1:OD2" &&
               compare_distance(e, 2.55251) && compare_energy(e, 3.18660);
    };
    auto e3 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               source_atom(e) == "NZ" && target_atom(e) == "CD:OE1:OE2" &&
               compare_distance(e, 3.22121) && compare_energy(e, 4.21956);
    };

    {
        Result r = SetUp("ionion/ionion3.pdb");

        EXPECT_EQ(r.count_edges(isIonicFunc), 3);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
        EXPECT_TRUE(r.contain_edge(e3));
    }
    {
        Result r = SetUp("ionion/ionion3.pdb", {"--ionic-bond", "3.2"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 2);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
    }
    {
        Result r = SetUp("ionion/ionion3.pdb", {"--ionic-bond", "1.8"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("ionion/ionion3.pdb", {"--ionic-bond", "1.59"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 0);
    }
}
TEST_F(BlackBoxTest, IonIon4) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               source_atom(e) == "NZ" && target_atom(e) == "OE1" &&
               compare_distance(e, 2.56018) && compare_energy(e, 5.30904);
    };

    {
        Result r = SetUp("ionion/ionion4.pdb");

        EXPECT_EQ(r.count_edges(isIonicFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("ionion/ionion4.pdb", {"--ionic-bond", "2.54"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 0);
    }
}
TEST_F(BlackBoxTest, IonIon5) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "GLU" &&
               source_atom(e) == "NZ" && target_atom(e) == "OE1" &&
               compare_distance(e, 1.95108) && compare_energy(e, 6.96644);
    };
    auto e2 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "ARG" && target(e) == "GLU" &&
               source_atom(e) == "NH2" && target_atom(e) == "OE1" &&
               compare_distance(e, 2.03282) && compare_energy(e, 2.71632);
    };

    {
        Result r = SetUp("ionion/ionion5.pdb");

        EXPECT_EQ(r.count_edges(isIonicFunc), 2);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
    }
    {
        Result r = SetUp("ionion/ionion5.pdb", {"--ionic-bond", "1.98"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("ionion/ionion5.pdb", {"--ionic-bond", "1.93"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 0);
    }
}
TEST_F(BlackBoxTest, IonIon6) {
    Result r = SetUp("ionion/ionion6.pdb");

    EXPECT_EQ(r.count_edges(isIonicFunc), 0);
}
TEST_F(BlackBoxTest, IonIon7) {
    Result r = SetUp("ionion/ionion7.pdb");

    EXPECT_EQ(r.count_edges(isIonicFunc), 0);
}
TEST_F(BlackBoxTest, IonIon8) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "IONIC" && source(e) == "LYS" && target(e) == "ASP" &&
               source_atom(e) == "NZ" && target_atom(e) == "CG:OD1:OD2" &&
               compare_distance(e, 3.19828) && compare_energy(e, 2.54320);
    };

    {
        Result r = SetUp("ionion/ionion8.pdb");

        EXPECT_EQ(r.count_edges(isIonicFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("ionion/ionion8.pdb", {"--ionic-bond", "3.1"});

        EXPECT_EQ(r.count_edges(isIonicFunc), 0);
    }
}
TEST_F(BlackBoxTest, IonIon9) {
    Result r = SetUp("ionion/ionion9.pdb");

    EXPECT_EQ(r.count_edges(isIonicFunc), 0);
}

#pragma endregion

#pragma region HBond

TEST_F(BlackBoxTest, HBond1) {
    Result r = SetUp("hbond/hbond1.pdb");

    EXPECT_EQ(r.count_edges(isHbondFunc), 0);
}
TEST_F(BlackBoxTest, HBond2) {
    Result r = SetUp("hbond/hbond2.pdb");

    EXPECT_EQ(r.count_edges(isHbondFunc), 0);
}
TEST_F(BlackBoxTest, HBond3) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "ASP" && target(e) == "ARG" &&
               source_atom(e) == "OD1" && target_atom(e) == "NH2" &&
               compare_distance(e, 2.11193) && compare_angle(e, 177.61115) && compare_energy(e, -6864.73748);
    };//angle_adh 1.23

    {
        Result r = SetUp("hbond/hbond3.pdb");

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond3.pdb", {"--hydrogen-bond", "2"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }

    {
        Result r = SetUp("hbond/hbond3.pdb", {"--h-bond-angle", "1.2"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }
}
TEST_F(BlackBoxTest, HBond4) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "GLN" && target(e) == "LYS" &&
               source_atom(e) == "OE1" && target_atom(e) == "NZ" &&
               compare_distance(e, 2.53981) && compare_angle(e, 179.91626) && compare_energy(e, -21052.13189);
    };//angle_adh 0.032959

    {
        Result r = SetUp("hbond/hbond4.pdb");

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond4.pdb", {"--hydrogen-bond", "2.4"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }

    {
        Result r = SetUp("hbond/hbond4.pdb", {"--h-bond-angle", "0.03"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }
}
TEST_F(BlackBoxTest, HBond5) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "ASN" && target(e) == "ASN" &&
               source_atom(e) == "OD1" && target_atom(e) == "ND2" &&
               compare_distance(e, 2.01970) && compare_angle(e, 179.96223) && compare_energy(e, -20882.11645);
    };//angle_adh 0.018706
    auto e2 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "ASN" && target(e) == "ASN" &&
               source_atom(e) == "OD1" && target_atom(e) == "ND2" &&
               compare_distance(e, 2.02333) && compare_angle(e, 179.96160) && compare_energy(e, -20920.96430);
    };//angle_adh 0.018979

    {
        Result r = SetUp("hbond/hbond5.pdb");

        EXPECT_EQ(r.count_edges(isHbondFunc), 2);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
    }
    {
        Result r = SetUp("hbond/hbond5.pdb", {"--hydrogen-bond", "2.02"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond5.pdb", {"--hydrogen-bond", "2"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }

    {
        Result r = SetUp("hbond/hbond5.pdb", {"--h-bond-angle", "0.0189"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond5.pdb", {"--h-bond-angle", "0.0186"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }
}
TEST_F(BlackBoxTest, HBond6) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" &&
               source_atom(e) == "OG" && target_atom(e) == "NZ" &&
               compare_distance(e, 2.00981) && compare_angle(e, 179.94896) && compare_energy(e, -31329.57651);
    };//angle_adh 0.025399
    auto e2 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" &&
               source_atom(e) == "OG" && target_atom(e) == "NZ" &&
               compare_distance(e, 2.02536) && compare_angle(e, 179.95765) && compare_energy(e, -31379.24148);
    };//angle_adh 0.020911
    auto e3 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" &&
               source_atom(e) == "OG" && target_atom(e) == "NZ" &&
               compare_distance(e, 2.02752) && compare_angle(e, 179.97440) && compare_energy(e, -31293.78545);
    };//angle_adh 0.012631

    {
        Result r = SetUp("hbond/hbond6.pdb");

        EXPECT_EQ(r.count_edges(isHbondFunc), 3);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
        EXPECT_TRUE(r.contain_edge(e3));
    }
    {
        Result r = SetUp("hbond/hbond6.pdb", {"--hydrogen-bond", "2.026"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 2);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
    }
    {
        Result r = SetUp("hbond/hbond6.pdb", {"--hydrogen-bond", "2.01"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond6.pdb", {"--hydrogen-bond", "2"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }

    {
        Result r = SetUp("hbond/hbond6.pdb", {"--h-bond-angle", "0.025"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 2);
        EXPECT_TRUE(r.contain_edge(e2));
        EXPECT_TRUE(r.contain_edge(e3));
    }
    {
        Result r = SetUp("hbond/hbond6.pdb", {"--h-bond-angle", "0.020"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e3));
    }
    {
        Result r = SetUp("hbond/hbond6.pdb", {"--h-bond-angle", "0.012"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }
}
TEST_F(BlackBoxTest, HBond7) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "SER" && target(e) == "LYS" && source_atom(e) == "OG" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.02752) && compare_angle(e, 179.97440) && compare_energy(e, -31293.78545);
    };//angle_adh 0.012631

    {
            Result r = SetUp("hbond/hbond7.pdb");

            EXPECT_EQ(r.count_edges(isHbondFunc), 1);
            EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond7.pdb", {"--hydrogen-bond", "2.01"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }

    {
        Result r = SetUp("hbond/hbond7.pdb", {"--h-bond-angle", "0.012"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }
}
TEST_F(BlackBoxTest, HBond8) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "HBOND" && source(e) == "GLN" && target(e) == "LYS" && source_atom(e) == "OE1" &&
               target_atom(e) == "NZ" &&
               compare_distance(e, 2.56130) && compare_angle(e, 67.88154) && compare_energy(e, 0.64979);
    };//angle_adh 0.649794

    {
        Result r = SetUp("hbond/hbond8.pdb");

        EXPECT_EQ(r.count_edges(isHbondFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("hbond/hbond8.pdb", {"--hydrogen-bond", "2.55"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }

    {
        Result r = SetUp("hbond/hbond8.pdb", {"--h-bond-angle", "0.64"});

        EXPECT_EQ(r.count_edges(isHbondFunc), 0);
    }
}

#pragma endregion

#pragma region PiPi

TEST_F(BlackBoxTest, PiPi1) {
    Result r = SetUp("pipi/pipi1.pdb");

    EXPECT_EQ(r.count_edges(isPipiFunc), 0);
}
TEST_F(BlackBoxTest, PiPi2) {
    Result r = SetUp("pipi/pipi2.pdb");

    EXPECT_EQ(r.count_edges(isPipiFunc), 0);
}
TEST_F(BlackBoxTest, PiPi3) {
    Result r = SetUp("pipi/pipi3.pdb");

    EXPECT_EQ(r.count_edges(isPipiFunc), 0);
}
TEST_F(BlackBoxTest, PiPi4) {
    Result r = SetUp("pipi/pipi4.pdb");

    EXPECT_EQ(r.count_edges(isPipiFunc), 0);
}
TEST_F(BlackBoxTest, PiPi5) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "PIPISTACK" && source(e) == "PHE" && target(e) == "PHE" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" &&
               compare_distance(e, 1) && compare_angle(e, 0) && compare_energy(e, -0.52740);
    };//normal center angles = 0 - 0

    {
        Result r = SetUp("pipi/pipi5.pdb");

        EXPECT_EQ(r.count_edges(isPipiFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
}

TEST_F(BlackBoxTest, PiPi6) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "PIPISTACK" && source(e) == "PHE" && target(e) == "PHE" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" &&
               compare_distance(e, 1.99982) && compare_angle(e, 0) && compare_energy(e, -0.52740);
    };//normal center angles = 45.003376 - 45.003376

    {
        Result r = SetUp("pipi/pipi6.pdb");

        EXPECT_EQ(r.count_edges(isPipiFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }

    {
        Result r = SetUp("pipi/pipi6.pdb", {"--pipistack-bond", "1.9"});

        EXPECT_EQ(r.count_edges(isPipiFunc), 0);
    }
    {
        Result r = SetUp("pipi/pipi6.pdb", {"--pipistack-normal-centre", "45"});

        EXPECT_EQ(r.count_edges(isPipiFunc), 0);
    }
}

#pragma endregion

#pragma region VDW

TEST_F(BlackBoxTest, Vdw1) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "ASN" && target(e) == "GLN" &&
               source_atom(e) == "ND2" && target_atom(e) == "OE1" &&
               compare_distance(e, 3.30000) && compare_energy(e, -0.16185);
    };//surface distance 0.150000

    {
        Result r = SetUp("vdw/vdw1.pdb");

        EXPECT_EQ(r.count_edges(isVdwFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("vdw/vdw1.pdb", {"--vdw-bond", "3.2"});

        EXPECT_EQ(r.count_edges(isVdwFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
}
TEST_F(BlackBoxTest, Vdw2) {
    Result r = SetUp("vdw/vdw2.pdb");
    EXPECT_EQ(r.count_edges(isVdwFunc), 0);
}
TEST_F(BlackBoxTest, Vdw3) {
    Result r = SetUp("vdw/vdw3.pdb");
    EXPECT_EQ(r.count_edges(isVdwFunc), 0);
}
TEST_F(BlackBoxTest, Vdw4) {
    Result r = SetUp("vdw/vdw4.pdb");
    EXPECT_EQ(r.count_edges(isVdwFunc), 0);
}
TEST_F(BlackBoxTest, Vdw5) {
    Result r = SetUp("vdw/vdw5.pdb");
    EXPECT_EQ(r.count_edges(isVdwFunc), 0);
}
TEST_F(BlackBoxTest, Vdw6) {
    Result r = SetUp("vdw/vdw6.pdb");
    EXPECT_EQ(r.count_edges(isVdwFunc), 0);
}
TEST_F(BlackBoxTest, Vdw7) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "VDW" && target(e) == "ASN" && source(e) == "GLN" &&
               target_atom(e) == "CB" && source_atom(e) == "NE2" &&
               compare_distance(e, 3.85866) && compare_energy(e, -0.13357);
    };//surface dist 0.488660

    {
        Result r = SetUp("vdw/vdw7.pdb");

        EXPECT_EQ(r.count_edges(isVdwFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("vdw/vdw7.pdb", {"--vdw-bond", "0.48"});

        EXPECT_EQ(r.count_edges(isVdwFunc), 0);
    }
}
TEST_F(BlackBoxTest, Vdw8) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "ASN" && target(e) == "ASN" &&
               source_atom(e) == "ND2" && target_atom(e) == "CB" &&
               compare_distance(e, 3.74499) && compare_energy(e, -0.10873);
    };//surface dist 0.374989
    auto e2 = [](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "GLN" && target(e) == "ASN" &&
               source_atom(e) == "OE1" && target_atom(e) == "CB" &&
               compare_distance(e, 1.73218) && compare_energy(e, 2022.14052);
    };//surface dist -1.587821
    auto e3 = [](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "GLN" && target(e) == "GLN" &&
               source_atom(e) == "OE1" && target_atom(e) == "NE2" &&
               compare_distance(e, 3.49130) && compare_energy(e, -0.18889);
    };//surface dist 0.341303
    auto e4 = [](const edge &e) {
        return interaction_name(e) == "VDW" && source(e) == "ASN" && target(e) == "GLN" &&
               source_atom(e) == "CB" && target_atom(e) == "NE2" &&
               compare_distance(e, 1.75928) && compare_energy(e, 2653.90001);
    };//surface dist -1.610719

    {
        Result r = SetUp("vdw/vdw8.pdb");

        EXPECT_EQ(r.count_edges(isVdwFunc), 4);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2));
        EXPECT_TRUE(r.contain_edge(e3));
        EXPECT_TRUE(r.contain_edge(e4));
    }
    {
        Result r = SetUp("vdw/vdw8.pdb", {"--vdw-bond", "0.37"});

        EXPECT_EQ(r.count_edges(isVdwFunc), 3);
        EXPECT_TRUE(r.contain_edge(e2));
        EXPECT_TRUE(r.contain_edge(e3));
        EXPECT_TRUE(r.contain_edge(e4));
    }
    {
        Result r = SetUp("vdw/vdw8.pdb", {"--vdw-bond", "0.34"});

        EXPECT_EQ(r.count_edges(isVdwFunc), 2);
        EXPECT_TRUE(r.contain_edge(e2));
        EXPECT_TRUE(r.contain_edge(e4));
    }
    {
        Result r = SetUp("vdw/vdw8.pdb", {"--vdw-bond", "-1.59"});

        EXPECT_EQ(r.count_edges(isVdwFunc), 1);
        EXPECT_TRUE(r.contain_edge(e4));
    }
    {
        Result r = SetUp("vdw/vdw8.pdb", {"--vdw-bond", "-1.62"});

        EXPECT_EQ(r.count_edges(isVdwFunc), 0);
    }
}

#pragma endregion

#pragma region Pication

TEST_F(BlackBoxTest, Picat1) {
    Result r = SetUp("picat/picat1.pdb");


    EXPECT_EQ(r.count_edges(isPicatFunc), 0);
}
TEST_F(BlackBoxTest, Picat2) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 3.99983) && compare_angle(e, 60.00134) && compare_energy(e, -0.74231);
    };
    auto e2_3 = [](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 4) && compare_angle(e, 89.99761) && compare_energy(e, -0.74219);
    };
    auto e4 = [](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 4.00013) && compare_angle(e, 70.00445) && compare_energy(e, -0.74209);
    };

    {
        Result r = SetUp("picat/picat2.pdb");

        EXPECT_EQ(r.count_edges(isPicatFunc), 4);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_EQ(r.count_edges(e2_3), 2);
        EXPECT_TRUE(r.contain_edge(e4));
    }
    {
        Result r = SetUp("picat/picat2.pdb", {"--pication-bond", "4.0001"});

        EXPECT_EQ(r.count_edges(isPicatFunc), 3);
        EXPECT_TRUE(r.contain_edge(e1));
        EXPECT_TRUE(r.contain_edge(e2_3));
    }
    {
        Result r = SetUp("picat/picat2.pdb", {"--pication-bond", "3.9999"});

        EXPECT_EQ(r.count_edges(isPicatFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("picat/picat2.pdb", {"--pication-bond", "3.9"});

        EXPECT_EQ(r.count_edges(isPicatFunc), 0);
    }
}
TEST_F(BlackBoxTest, Picat3) {
    Result r = SetUp("picat/picat3.pdb");

    EXPECT_EQ(r.count_edges(isPicatFunc), 0);
}
TEST_F(BlackBoxTest, Picat4) {
    Result r = SetUp("picat/picat4.pdb");

    EXPECT_EQ(r.count_edges(isPicatFunc), 0);
}
TEST_F(BlackBoxTest, Picat5) {
    auto e1 = [](const edge &e) {
        return interaction_name(e) == "PICATION" && source(e) == "TYR" && target(e) == "LYS" &&
               source_atom(e) == "CD1:CD2:CE1:CE2:CG:CZ" && target_atom(e) == "NZ" &&
               compare_distance(e, 2.99985) && compare_angle(e, 70.00357) && compare_energy(e, -2.34615);
    };

    {
        Result r = SetUp("picat/picat5.pdb");

        EXPECT_EQ(r.count_edges(isPicatFunc), 1);
        EXPECT_TRUE(r.contain_edge(e1));
    }
    {
        Result r = SetUp("picat/picat5.pdb", {"--pication-bond", "2.98"});

        EXPECT_EQ(r.count_edges(isPicatFunc), 0);
    }
}

#pragma endregion