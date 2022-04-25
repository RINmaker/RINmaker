#include "rin_graph.h"

#include <filesystem>
#include <string>
#include <utility>

#include "bonds.h"
#include "chemical_entity.h"
#include "config.h"

namespace fs = std::filesystem;

using std::vector, std::string, std::queue, std::to_string, std::unordered_map;
using chemical_entity::aminoacid;

using namespace rin;

struct edge::impl final
{
public:
    string _source, _target, _source_atom, _target_atom;
    string _distance, _energy, _angle, _interaction, _orientation;
    string _donor, _cation, _positive;
};

edge::edge(edge const& other) : pimpl{new impl(*other.pimpl)}
{}

edge::~edge()
{ delete pimpl; }

string const& edge::source_id() const
{ return pimpl->_source; }

string const& edge::target_id() const
{ return pimpl->_target; }

string const& edge::distance() const
{ return pimpl->_distance; }

string const& edge::energy() const
{ return pimpl->_energy; }

string const& edge::interaction() const
{ return pimpl->_interaction; }

string const& edge::source_atom() const
{ return pimpl->_source_atom; }

string const& edge::target_atom() const
{ return pimpl->_target_atom; }

string const& edge::angle() const
{ return pimpl->_angle; }

string const& edge::donor() const
{ return pimpl->_donor; }

string const& edge::cation() const
{ return pimpl->_cation; }

string const& edge::positive() const
{ return pimpl->_positive; }

string const& edge::orientation() const
{ return pimpl->_orientation; }

graph::graph(string name, rin::parameters const& params, vector<aminoacid const*> const& aminoacids,
             vector<std::shared_ptr<bond::base const>> const& bonds) : _name(std::move(name)), _params(params)
{
    for (auto a: aminoacids)
    {
        auto n = (rin::node) *a;
        _nodes.insert({n.get_id(), n});
    }

    // adjust _nodes degree at edge insertion
    for (auto b: bonds)
    {
        auto edge = (rin::edge) *b;

        auto it = _nodes.find(edge.source_id());
        if (it != _nodes.end())
            it->second.inc_degree();

        it = _nodes.find(edge.target_id());
        if (it != _nodes.end())
            it->second.inc_degree();

        _edges.push_back(edge);
    }
}

void graph::write_to_file(std::filesystem::path const& out_path) const
{
    pugi::xml_document doc;

    // <graphml>
    pugi::xml_node graphml = doc.append_child("graphml");
    graphml.append_attribute("xmlns") = "http://graphml.graphdrawing.org/xmlns";
    graphml.append_attribute("xmlns:xsi") = "http://www.w3.org/2001/XMLSchema-instance";
    graphml.append_attribute("xsi:schemaLocation") =
            "http://graphml.graphdrawing.org/xmlns http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd";

    // parameters summary
    auto comment = graphml.parent().insert_child_before(pugi::node_comment, graphml);
    comment.set_value(_params.pretty().c_str());

    // <graph>
    pugi::xml_node graph_node = graphml.append_child("graph");
    graph_node.append_attribute("id") = name().c_str();
    graph_node.append_attribute("edgedefault") = "undirected";

    // graphml requires all key attributes to be listed before the actual node/edges
    bool with_metadata = true;
    for (auto e: _edges)
    {
        e.append_to(graph_node, with_metadata);
        if (with_metadata) with_metadata = false;
    }

    with_metadata = true;
    for (auto kv: _nodes)
    {
        if (kv.second.degree() > 0)
        {
            kv.second.append_to(graph_node, with_metadata);
            if (with_metadata) with_metadata = false;
        }
    }

    doc.save_file(out_path.c_str());
}

edge::edge(bond::ss const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source_id();
    pimpl->_target = bond.target_id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = std::to_string(bond.get_energy());
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = "SG"; // TODO config
    pimpl->_target_atom = "SG"; // TODO config
    pimpl->_donor = cfg::graphml::none;
    pimpl->_cation = cfg::graphml::none;
    pimpl->_positive = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
    pimpl->_angle = cfg::graphml::null;
}

edge::edge(bond::vdw const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source_atom().res().id();
    pimpl->_target = bond.target_atom().res().id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = std::to_string(bond.get_energy());
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = bond.source_atom().name();
    pimpl->_target_atom = bond.target_atom().name();
    pimpl->_donor = cfg::graphml::none;
    pimpl->_angle = cfg::graphml::null;
    pimpl->_cation = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
    pimpl->_positive = cfg::graphml::none;
}

edge::edge(bond::ionic const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source_positive().res().id();
    pimpl->_target = bond.target_negative().res().id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = std::to_string(bond.get_energy());
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = bond.source_positive().name();
    pimpl->_target_atom = bond.target_negative().name();
    pimpl->_positive = bond.source_positive().res().id();
    pimpl->_angle = cfg::graphml::null;
    pimpl->_donor = cfg::graphml::none;
    pimpl->_cation = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
}

edge::edge(bond::hydrogen const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source_atom().res().id();
    pimpl->_target = bond.target_atom().res().id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = std::to_string(bond.get_energy());
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = bond.source_atom().name();
    pimpl->_target_atom = bond.target_atom().name();
    pimpl->_angle = std::to_string(bond.get_angle());
    pimpl->_donor = bond.donor().res().id();
    pimpl->_cation = cfg::graphml::none;
    pimpl->_positive = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
}

edge::edge(bond::pipistack const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source_ring().res().id();
    pimpl->_target = bond.target_ring().res().id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = std::to_string(bond.get_energy());
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = bond.source_ring().name();
    pimpl->_target_atom = bond.target_ring().name();
    pimpl->_angle = std::to_string(bond.angle());
    pimpl->_donor = cfg::graphml::none;
    pimpl->_cation = cfg::graphml::none;
    pimpl->_positive = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
}

edge::edge(bond::pication const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source_ring().res().id();
    pimpl->_target = bond.target_cation().res().id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = std::to_string(bond.get_energy());
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = bond.source_ring().name();
    pimpl->_target_atom = bond.target_cation().name();
    pimpl->_cation = bond.target_cation().res().id();
    pimpl->_angle = std::to_string(bond.angle());
    pimpl->_donor = cfg::graphml::none;
    pimpl->_positive = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
}

edge::edge(bond::generic_bond const& bond) : pimpl{new impl()}
{
    pimpl->_source = bond.source().id();
    pimpl->_target = bond.target().id();
    pimpl->_distance = std::to_string(bond.get_length());
    pimpl->_energy = cfg::graphml::none;
    pimpl->_interaction = bond.get_interaction();
    pimpl->_source_atom = bond.source().name();
    pimpl->_target_atom = bond.target().name();
    pimpl->_angle = cfg::graphml::null;
    pimpl->_donor = cfg::graphml::none;
    pimpl->_cation = cfg::graphml::none;
    pimpl->_positive = cfg::graphml::none;
    pimpl->_orientation = cfg::graphml::none;
}

void add_data(
        pugi::xml_node& node, string const& prefix, string const& type, string const& key_name, string const& key_value,
        string const& key_type, bool with_metadata)
{
    pugi::xml_node data = node.append_child("data");
    data.append_attribute("key") = (prefix + key_name).c_str();
    data.append_child(pugi::node_pcdata).set_value(key_value.c_str());

    if (with_metadata)
    {
        // add key info in (node -> graph -> graphml) before (node -> graph)
        pugi::xml_node key = node.parent().parent().insert_child_before("key", node.parent());
        key.append_attribute("id") = (prefix + key_name).c_str();
        key.append_attribute("for") = type.c_str();
        key.append_attribute("attr.name") = key_name.c_str();
        key.append_attribute("attr.type") = key_type.c_str();
    }
}

void edge::append_to(pugi::xml_node& rin, bool metadata)
{
    // the xml node representing a rin edge
    pugi::xml_node edge = rin.append_child("edge");
    edge.append_attribute("source") = pimpl->_source.c_str();
    edge.append_attribute("target") = pimpl->_target.c_str();

    add_data(edge, "e_", "edge", "NodeId1", pimpl->_source, "string", metadata);
    add_data(edge, "e_", "edge", "NodeId2", pimpl->_target, "string", metadata);

    add_data(edge, "e_", "edge", "Energy", pimpl->_energy, "double", metadata);
    add_data(edge, "e_", "edge", "Distance", pimpl->_distance, "double", metadata);

    add_data(edge, "e_", "edge", "Interaction", pimpl->_interaction, "string", metadata);
    add_data(edge, "e_", "edge", "Atom1", pimpl->_source_atom, "string", metadata);
    add_data(edge, "e_", "edge", "Atom2", pimpl->_target_atom, "string", metadata);

    add_data(edge, "e_", "edge", "Angle", pimpl->_angle, "double", metadata);
    add_data(edge, "e_", "edge", "Donor", pimpl->_donor, "string", metadata);
    add_data(edge, "e_", "edge", "Cation", pimpl->_cation, "string", metadata);
    add_data(edge, "e_", "edge", "Positive", pimpl->_positive, "string", metadata);
    add_data(edge, "e_", "edge", "Orientation", pimpl->_orientation, "string", metadata);
}

struct node::impl final
{
public:
    string _id;
    string _pdb_name;
    string _chain;
    string _seq;
    string _name;
    string _x, _y, _z;
    string _bfactor;
    string _secondary;

    int _degree = 0;
};


node::node(chemical_entity::aminoacid const& res) : pimpl{new impl()}
{
    pimpl->_id = res.id();
    pimpl->_chain = res.chain_id();
    pimpl->_seq = std::to_string(res.sequence_number());
    pimpl->_name = res.name();
    pimpl->_x = std::to_string(res[0]);
    pimpl->_y = std::to_string(res[1]);
    pimpl->_z = std::to_string(res[2]);
    pimpl->_bfactor = res.ca() == nullptr ? "NULL" : std::to_string(res.ca()->temp_factor());
    pimpl->_secondary = res.secondary_structure_id();
    pimpl->_pdb_name = res.pdb_name();
    pimpl->_degree = 0;
}

node::node(node const& other) : pimpl{new impl(*other.pimpl)}
{}

node::~node()
{ delete pimpl; }

void node::inc_degree()
{ ++pimpl->_degree; }

int node::degree() const
{ return pimpl->_degree; }

[[nodiscard]]
string const& node::get_id() const
{ return pimpl->_id; }


void node::append_to(pugi::xml_node& graph, bool with_metadata) const
{
    pugi::xml_node node;

    node = graph.prepend_child("node");
    node.append_attribute("id") = pimpl->_id.c_str();

    add_data(node, "v_", "node", "Degree", std::to_string(pimpl->_degree), "double", with_metadata);
    add_data(node, "v_", "node", "NodeId", pimpl->_id, "string", with_metadata);

    add_data(node, "v_", "node", "Residue", pimpl->_id, "string", with_metadata);
    add_data(node, "v_", "node", "Chain", pimpl->_chain, "string", with_metadata);
    add_data(node, "v_", "node", "Position", pimpl->_seq, "double", with_metadata);
    add_data(node, "v_", "node", "Name", pimpl->_name, "string", with_metadata);

    add_data(node, "v_", "node", "x", pimpl->_x, "double", with_metadata);
    add_data(node, "v_", "node", "y", pimpl->_y, "double", with_metadata);
    add_data(node, "v_", "node", "z", pimpl->_z, "double", with_metadata);

    add_data(node, "v_", "node", "Bfactor_CA", pimpl->_bfactor, "double", with_metadata);
    add_data(node, "v_", "node", "Secondary_Structure", pimpl->_secondary, "string", with_metadata);

    add_data(node, "v_", "node", "PdbName", pimpl->_pdb_name, "string", with_metadata);
}
